#' Run Multilevel Modeling (MLM) Diagnostics
#'
#' @description This function performs descriptive analysis, calculates ICC and Design Effects, 
#' and generates diagnostic plots for multilevel data structures.
#'
#' @param dataFrame A data.frame containing the variables to be analyzed.
#' @param maxCat Integer. The maximum number of unique values a variable can have to be treated as categorical.
#' @param displayVar Character. The name of the grouping variable to be used for plot labels. Defaults to NULL.
#' @param groupVars Character vector. The names of the grouping/clustering variables. Defaults to NULL.
#' @param multiLvl Logical. If TRUE, calculates multilevel statistics (ICC/DEFF). Defaults to TRUE.
#' @param minCase Integer. Minimum number of usable cases required for a variable pair to be analyzed. Defaults to 100.
#'
#' @return A list of class 'mlm_diag' containing:
#' \itemize{
#'   \item Plots: A list of ggplot2 objects for each outcome variable.
#'   \item Correlation_Matrix: A data.frame of pairwise correlations.
#'   \item ICC_Table: A data.frame of ICC and Design Effect values (if multiLvl is TRUE).
#' }
#' @export
#'
#' @examples
#' # mlmDiagnostics(my_data, maxCat = 5, groupVars = "school_id")


mlmDiagnostics <- function(dataFrame, maxCat, displayVar = NULL, groupVars = NULL, multiLvl = TRUE, minCase = 100){ # if display var is null just don't put anything on the graph
  
  missingVars <- unique(c(
      setdiff(groupVars, names(dataFrame)),
      setdiff(displayVar, names(dataFrame))
      ))
  
  if(length(missingVars) > 0){
    stop(
      "Critical error: ", 
       paste0(missingVars, collapse = ", "),
      " not found in data.")
  }

  varNames <- names(dataFrame)
  outcomeVars <- setdiff(names(dataFrame), groupVars)


  dataFrame[outcomeVars] <- lapply(
  dataFrame[outcomeVars],
  makeNumer
  )

  dataFrame[groupVars] <- lapply(
  dataFrame[groupVars],
  makeFactor
  )  

  corrMat <- buildcorrMat(df = dataFrame, vars = outcomeVars, maxCat = maxCat)

  iccTable <- NULL 

  if(multiLvl == TRUE){
  
    varPairs<- expand.grid(groupVar = groupVars, outVar = outcomeVars, stringsAsFactors = FALSE)
  
    iccTable <- buildiccTable(df = dataFrame, pairs = varPairs, minCase = minCase)
  
  }

  plots <- lapply(outcomeVars, function(varName) {  
  
    var <- dataFrame[[varName]]
  
    label <- NULL 
  
    if(multiLvl == TRUE){
    
    label <- mlmLabel(iccTable = iccTable, displayVar = displayVar, varName = varName) 
    
  }
 
 isCat <- checkCat(var = var, maxCat = maxCat)
 
 propdf <- makedfProp(var = var, isCat = isCat)  
 
 makevarPlot(type = isCat, dfProp = propdf, var = var, varName = varName, iccLabel = label) 
})


  return(list(
    Plots = plots,
    Correlation_Matrix = corrMat,
    ICC_Table = iccTable
))
}



buildcorrMat <- function(df, vars, maxCat){
  
  pairs <- utils::combn(vars, 2)
  
  corList <- apply(pairs, 2, function(p){
    v1 <- p[1]
    v2 <- p[2]
   
    method <- if (!checkCat(var = df[[v1]], maxCat = maxCat) && !checkCat(var = df[[v2]], maxCat = maxCat)) "pearson" else "spearman"
    
    r <- stats::cor(df[[v1]], df[[v2]], use = "pairwise.complete.obs", method = method)
    data.frame(Var1 = v1, Var2 = v2, Correlation = r)
  } )
  
  corDf <- do.call(rbind, corList)
  return(corDf)
  
}


buildiccTable <- function(df, pairs, minCase){
  
  iccList <- Map(function(outVar, groupVar, df, minCase){
    
    outCol <- df[[outVar]]
    groupCol <- df[[groupVar]]
    
    ok  <- checkUsablevar(column = outCol, group = groupCol, varName = outVar, groupName = groupVar, minCase = minCase) 
    if (is.null(ok)) return(NULL)
    
    subDf <- df[ok, , drop = FALSE]
    
    suboutCol <- subDf[[outVar]]
    subgroupCol <- subDf[[groupVar]]
    
    modelBinary  <- isBinary(var = suboutCol) 
    fit <- fitMLM(outcome = outVar, data = subDf, group = groupVar, isBinary = modelBinary)
    icc <- iccFit(fit = fit)
    deseff <- designEfffit(fit = fit, data = subDf,group = groupVar)
    
    return(data.frame(
      Outcome = outVar,
      Group = groupVar,
      ICC = icc,
      Design_Effect = deseff,
      stringsAsFactors = FALSE
    ))
    
  }, outVar = pairs$outVar, groupVar = pairs$groupVar, MoreArgs = list(df = df, minCase = minCase))
  
  iccDf <- do.call(rbind, iccList)
  
  return(iccDf)
}


isBinary <- function(var) {
  length(unique(var[!is.na(var)])) == 2
}


makevarPlot <- function(type, dfProp, var, varName, iccLabel = NULL) {
  
  if (type) {
   
    p <- ggplot2::ggplot(dfProp, ggplot2::aes(x = as.factor(category), y = proportion)) +
      ggplot2::geom_col() +
      ggplot2::labs(
        x = varName,
        y = "Proportion",
        title = "Bar chart"
      )
    
  } else {
    
    descrip <- computeDesc(var = var)
    
    descripLabel <- paste0(
      "mean: ", sprintf("%.3f", descrip["mean"]), "\n",
      "sd: ", sprintf("%.3f", descrip["sd"]), "\n",
      "skewness: ", sprintf("%.3f", descrip["skewness"]), "\n",
      "kurtosis: ", sprintf("%.3f", descrip["kurtosis"]), "\n",
      "min: ", sprintf("%.3f", descrip["min"]), "\n",
      "max: ", sprintf("%.3f", descrip["max"])
    )
    
    dfCont <- data.frame(x = var)
    
    p <- ggplot2::ggplot(dfCont, ggplot2::aes(x = x)) +
      ggplot2::geom_histogram() +
      ggplot2::geom_vline(xintercept = descrip["mean"], linetype = "dashed") +
      ggplot2::geom_vline(xintercept = descrip["min"],  linetype = "dotted") +
      ggplot2::geom_vline(xintercept = descrip["max"],  linetype = "dotted") +
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
  
  
  if (!is.null(iccLabel)) {
    p <- p + ggplot2::annotate(
      "text",
      x = -Inf, y = Inf,
      label = iccLabel,
      hjust = -0.05, vjust = 1.1
    )
  }
  
  return(p)
  
}



computeDesc <- function(var) {

  x <- var[!is.na(var)]

  c(
  mean     = stats::mean(x),
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
    category = as.character(uni),
    proportion = as.numeric(props),
    stringsAsFactors = FALSE
    )
} else return(NULL)
}


extractVarcom <- function(fit) {
  
  if (!inherits(fit, c("lmerMod", "glmerMod")))
    
    stop("fit must be lme4 model")
  
    vc <- as.data.frame(lme4::VarCorr(fit))
    
    sigmaU2 <- vc$vcov[vc$grp != "Residual" & vc$var1 == "(Intercept)"]  


  if (inherits(fit, "lmerMod")) {
    sigmaE2 <- vc$vcov[vc$grp == "Residual"]              
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


mlmLabel <- function(iccTable, displayVar, varName){
  
  if(is.null(displayVar)){ return(NULL) }
  
    display <- iccTable[iccTable$Group == displayVar, ]

    idx <- display$Outcome == varName 
    
  if(!any(idx, na.rm = TRUE)){
    
    return(paste0(displayVar, "\n No ICC data found"))
    
  }
  
  paste0(
    displayVar,
    "\n",
    "ICC: ", sprintf("%.3f", display[idx, "ICC"]), "\n",
    "DEFF: ", sprintf("%.3f", display[idx, "Design_Effect"])
)
}


iccFit <- function(fit) {
  comp <- extractVarcom(fit = fit)
  comp$sigmaU2 / (comp$sigmaU2 + comp$sigmaE2)
}


designEfffit <- function(fit, data, group) {
  icc <- iccFit(fit = fit)
  mbar <- mean(as.numeric(table(data[[group]])))
  1 + (mbar - 1) * icc
}


fitMLM <- function(outcome, data, group, isBinary,  link = "logit", REML = FALSE) {
  
  if (!is.character(group) || length(group) != 1L) { 
    
    stop("`group` must be a single column name (character scalar).")
  }
  
  if (!group %in% names(data)) {
    
    stop("`group` is not a column in `data`.")
    
  }

  if (!is.factor(data[[group]])) {
    
    stop("`data[[group]]` must be a factor.") 
    
  }

  full_formula <- stats::as.formula(
    
  paste0("`", outcome, "` ~ 1 + (1 | `", group, "`)")

  )

  if (isTRUE(isBinary)) {
    
    data[[outcome]] <- factor(data[[outcome]])
    
    return(lme4::glmer(formula = full_formula, data = data, family = stats::binomial(link = link)))
  }

  return(lme4::lmer(formula = full_formula, data = data, REML = REML))
}


checkCat <- function(var, maxCat){
  
  l <- length(unique(var[!is.na(var)]))
  
  l <= maxCat
}



checkUsablevar <- function(column, group, varName, groupName, minCase) {
  
  ok <- !is.na(column) & !is.na(group)
  
  if (sum(ok) < minCase) {
    
    warning(
      "Skipping pair: ", varName, groupName,
      ": too few usable cases (n = ", sum(ok), ")"
    )
    
    return(NULL)
    
  }
  
  if (length(unique(group[ok])) < 2) {
    
    warning(
      "Skipping pair: ", varName, groupName,
      ": < 2 groups with data"
    )
    
    return(NULL)
  }
  
  ok
}


makeFactor <- function(x) {

  if (inherits(x, "haven_labelled")) {
    
    x <- haven::zap_labels(x)
    
  }

  factor(x)
  
}


makeNumer <- function(x) {

  if (inherits(x, "haven_labelled")) {
    
    x <- haven::zap_labels(x)
    
  }

  if (is.factor(x)) {
    
    x <- as.character(x)
    
  }

  as.numeric(x)
}

