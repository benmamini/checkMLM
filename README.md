# checkMLM
Lightweight utilities for **diagnostic multilevel modeling** in R.

The main entry point, `mlmDiagnostics`, fits empty multilevel models (random intercept only) for multiple outcome variables, computes **ICC** and **design effects**, and produces basic diagnostic plots (histograms or bar charts) annotated with these quantities. Outputs pairwise correlation table and ICC/Design Effect table.

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

##Usage
```R
library(checkMLM)

# Run the diagnostic pipeline on your nested data
results <- mlmDiagnostics(
Dataframe = data, 
maxCat = 5,
displayVar = "school_id",
groupVar = "school_id", 
multiLvl = TRUE,
maxCat = 5 
)
```
## The function returns a list containing your plots and tables

```
print(results$ICC_Table)
```
##Project Structure

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
```
