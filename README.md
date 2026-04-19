# checkMLM
Lightweight utilities for **diagnostic multilevel modeling** in R.

mlmDiagnostics() is the main entry point for the package, designed to rapidly evaluate the nesting structure of your data before fitting complex mixed-effects models. It automatically runs random-intercept models across multiple outcomes to compute Intraclass Correlation Coefficients (ICC) and Design Effects. It also features optional Likelihood Ratio Testing (LRT) to evaluate random slopes. The function returns a clean, fault-tolerant list containing annotated diagnostic plots, spaghetti plots for slope variance, a pairwise correlation matrix, and comprehensive summary tables.

This is intended for **exploratory / screening use**, not full model specification.



## Data expectations

- `Dataframe` must be a data.frame
- `groupVars` A character vector of your cluster/nesting variables. All other variables in the dataframe are automatically treated as outcomes for ICC computations.
- `depVars`(Optional) For evaluating random slopes. All other variables (except groupVars) will be tested as predictors for likelihood ratio tests.
- Binary outcomes: Automatically detected and fit using logistic mixed-effects models.
- Missing Data: Allowed. The function performs listwise deletion per model. Variables with insufficient cases (defined by `minCase`) are skipped.
- Plots: Bar charts (using proportions) are generated for variables with ≤ maxCat unique values; histograms are generated for continuous variables.

## Dependencies

- ggplot2  
- lme4  
- moments  

## Usage
```R
library(checkMLM)

# Example 1: Basic ICC Screening
results <- mlmDiagnostics(
  dataFrame = my_data, 
  maxCat = 5, 
  groupVars = "school_id",
  displayVar = "school_id"
)

# View the ICC Table
results$ICC_Table

# View the diagnostic plot for the first variable
results$Plots[[1]]


# Example 2: Testing for Random Slopes
results_slopes <- mlmDiagnostics(
  dataFrame = my_data, 
  maxCat = 5, 
  groupVars = "school_id",
  slopes = TRUE,
  depVars = c("math_score", "reading_score")
)

# View the Likelihood Ratio Test results
results_slopes$Slope_Table

# View a spaghetti plot of slope variance
results_slopes$Slope_Plots[[1]]
```

## Project Structure


```
|   .gitignore
|   .Rbuildignore
|   .Rprofile
|   DESCRIPTION
|   LICENSE
|   NAMESPACE
|   README.md
|   renv.lock
|
+---man
|       mlmDiagnostics.Rd
|
+---R
|       checkMLM.R
|
\

