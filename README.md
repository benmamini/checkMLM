# checkMLM
Lightweight utilities for **diagnostic multilevel modeling** in R.

The main entry point, `checkMLM()`, fits empty multilevel models (random intercept only) for multiple outcome variables, computes **ICC** and **design effects**, and produces basic diagnostic plots (histograms or bar charts) annotated with these quantities. Outputs pairwise correlation table and ICC/Design Effect table.

This is intended for **exploratory / screening use**, not full model specification.



## Data expectations

- `Dataframe` must be a data.frame
- `groupVars` are the names of your cluster variables
- All other columns are treated as outcomes
- Missing values are allowed
- Categorical variables are defined as having ≤ `maxCat` unique values

## Dependencies

- ggplot2  
- lme4  
- moments  


```R
#Usage
library(checkMLM)

# Run the diagnostic pipeline on your nested data
results <- checkMLM(
  Dataframe = data, 
  groupVar = "school_id", 
  maxCat = 5 
)

# The function returns a list containing your plots and tables
print(results$ICC_Table)

