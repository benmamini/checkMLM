checkMLM <- function(uDataframe, groupVar, maxCat){
  
  dataFrame <- uDataframe
  outcomeVars <- setdiff(names(dataFrame), groupVar) 
 
  
  dataFrame[outcomeVars] <- lapply(
    dataFrame[outcomeVars],
    makeNumer
  )
  dataFrame[groupVar] <- lapply(
    dataFrame[groupVar],
    makeFactor
  )
  
group <- dataFrame[[groupVar]] 
  
    plots <- lapply(outcomeVars, function(varName) {
    var <- dataFrame[[varName]]
    
    ok  <- checkUsablevar(column = var, group = group, varName = varName)
    if (is.null(ok)) return(NULL)
    df <- dataFrame[ok, , drop = FALSE]
    var <- df[[varName]]
   
    
    
    isCat <- checkCat(var, maxCat) 
    

    modelBinary  <- isBinary(var)
    fit <- fitMLM(varName, df, groupVar, modelBinary)
    icc <- iccFit(fit)
    deseff <- designEfffit(fit, df, groupVar)
    label <- mlmLabel(icc, deseff)
    propdf <- makedfProp(var, isCat)
    makeVarplot(isCat, propdf, var, varName, label)
    # call your prewritten helpers here, using var/ok/varName/group
    # return a single ggplot object
  })
   
   
   
   plots
 }

#---------------------------------------------------------------
isBinary <- function(var) {
  length(unique(var[!is.na(var)])) == 2
}


makeVarplot <- function(type, dfProp, var, varName, iccLabel) {
  
  if (type) {
    
    ggplot2::ggplot(dfProp, ggplot2::aes(x = factor(x), y = n)) +
      ggplot2::geom_col() +
      ggplot2::annotate(
        "text",
        x = -Inf, y = Inf,
        label = iccLabel,
        hjust = -0.05, vjust = 1.1
      ) +
      ggplot2::labs(
        x = varName,
        y = "Proportion",
        title = "Bar chart"
      )
    
  } else {
    
    descrip <- computeDesc(var)
    descripLabel <- paste0(
      "mean: ", sprintf("%.3f", descrip["mean"]), "\n",
      "sd: ", sprintf("%.3f", descrip["sd"]), "\n",
      "skewness: ", sprintf("%.3f", descrip["skewness"]), "\n",
      "kurtosis: ", sprintf("%.3f", descrip["kurtosis"]), "\n",
      "min: ", sprintf("%.3f", descrip["min"]), "\n",
      "max: ", sprintf("%.3f", descrip["max"])
    )
    
    dfCont <- data.frame(x = var)
    
    ggplot2::ggplot(dfCont, ggplot2::aes(x = x)) +
      ggplot2::geom_histogram() +
      ggplot2::geom_vline(xintercept = descrip["mean"], linetype = "dashed") +
      ggplot2::geom_vline(xintercept = descrip["min"],  linetype = "dotted") +
      ggplot2::geom_vline(xintercept = descrip["max"],  linetype = "dotted") +
      ggplot2::annotate(
        "text",
        x = -Inf, y = Inf,
        label = iccLabel,
        hjust = -0.05, vjust = 1.1
      ) +
      ggplot2::annotate(
        "text",
        x = Inf, y = Inf,
        label = descripLabel,
        hjust = 1.05, vjust = 1.1
      ) +
      ggplot2::labs(
        x = varName,
        y = "Count",
        title = "Histogram"
      )
  }
}



computeDesc <- function(var) {
  
  x <- var[!is.na(var)]
  
  c(
    mean     = mean(x),
    sd       = stats::sd(x),
    skewness = moments::skewness(x),
    kurtosis = moments::kurtosis(x),  
    min      = min(x),
    max      = max(x)
  )
}


makedfProp <- function(var, isCat) {
if(isCat == TRUE){
  uni <- unique(var[!is.na(var)])
  tot <- sum(!is.na(var))
  
  props <- sapply(uni, function(u) {
    sum(var == u, na.rm = TRUE) / tot
  })
  
  data.frame(
    x = uni,
    n = props
  )
  } else return(NULL)
}


extractVarcom <- function(fit) {
  if (!inherits(fit, c("lmerMod", "glmerMod")))
    stop("fit must be lme4 model")
  
  
  vc <- as.data.frame(lme4::VarCorr(fit))
  sigmaU2 <- vc$vcov[vc$grp != "Residual" & vc$var1 == "(Intercept)"]  # between - takes group name as parameter
  
  
  
  
  if (inherits(fit, "lmerMod")) {
    sigmaE2 <- vc$vcov[vc$grp == "Residual"]                # WITHIN
  } else {
    link <- stats::family(fit)$link
    sigmaE2 <- switch(
      link,
      logit   = (pi^2) / 3,
      probit  = 1,
      cloglog = (pi^2) / 6,
      stop("Unsupported link for latent ICC: ", link)
    )
  }
  
  list(sigmaU2 = sigmaU2, sigmaE2 = sigmaE2)
}


mlmLabel <- function(icc, deff){
  paste0(
    "ICC: ", sprintf("%.3f", icc), "\n",
    "DEFF: ", sprintf("%.3f", deff)
  )
}

iccFit <- function(fit) {
  comp <- extractVarcom(fit)
  comp$sigmaU2 / (comp$sigmaU2 + comp$sigmaE2)
}


designEfffit <- function(fit, data, group) {
  icc <- iccFit(fit)
  mbar <- mean(as.numeric(table(data[[group]])))
  1 + (mbar - 1) * icc
}



fitMLM <- function(outcome,
                  data,
                  group,
                  isBinary, 
                  link = "logit",
                  REML = FALSE,
                  ...) {
  
  
  
  if (!is.character(group) || length(group) != 1L) { #ensure group (name of column) is character vector of length 1
    stop("`group` must be a single column name (character scalar).")
  }
  if (!group %in% names(data)) {
    stop("`group` is not a column in `data`.")
  }
  
  if (!is.factor(data[[group]])) {
    stop("`data[[group]]` must be a factor.") # we will convert to factor in the main function 
  }
  
  full_formula <- as.formula(
    paste0("`", outcome, "` ~ 1 + (1 | `", group, "`)")
  )
  
  if (isTRUE(isBinary)) {
    data[[outcome]] <- factor(data[[outcome]])
    return(lme4::glmer(full_formula, data = data, family = stats::binomial(link = link)))
  }
  
  lme4::lmer(full_formula, data = data, REML = REML)
}


checkCat <- function(var, maxCat){
  l <- length(unique(var[!is.na(var)]))
  l <= maxCat
}



checkUsablevar <- function(column, group, varName, min_n = 100) {
  ok <- !is.na(column) & !is.na(group)
  
  if (sum(ok) < min_n) {
    warning(
      "Skipping ", varName,
      ": too few usable cases (n = ", sum(ok), ")"
    )
    return(NULL)
  }
  
  if (length(unique(group[ok])) < 2) {
    warning(
      "Skipping ", varName,
      ": < 2 groups with data"
    )
    return(NULL)
  }
  
  ok
}


makeFactor <- function(x) {
  if (is.numeric(x)) factor(x)
  else if (inherits(x, "haven_labelled")) factor(as.numeric(x))
  else x
}

makeNumer <- function(x){
  
  if (is.factor(x)) as.numeric(as.character(x))
  else if (inherits(x, "haven_labelled")) as.numeric(x)
  else as.numeric(x)
  
}

