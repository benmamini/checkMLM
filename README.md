# checkMLM
Lightweight utilities for **diagnostic multilevel modeling** in R.

The main entry point, `checkMLM()`, fits empty multilevel models (random intercept only) for multiple outcome variables, computes **ICC** and **design effects**, and produces basic diagnostic plots (histograms or bar charts) annotated with these quantities.

This is intended for **exploratory / screening use**, not full model specification.



## Data expectations

- `uDataframe` must be a data.frame
- `groupVar` must be the name of a grouping variable (factor or coercible to factor)
- All other columns are treated as outcomes
- Missing values are allowed
- Categorical variables are defined as having ≤ `maxCat` unique values

## Dependencies

- ggplot2  
- lme4  
- moments  
