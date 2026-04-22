test_that("checkMod correctly identifies insufficient data", {
  # 1. Setup mock data
  mock_df <- data.frame(y = rnorm(10), g = rep(1:2, each = 5))
  
  # 2. Test for Failure (Skip Logic)
  # We have 10 cases, but we require 100. Should return NULL and a warning.
  expect_warning(
    res_null <- checkMod(df = mock_df, allCols = c("y", "g"), 
                         groupVars = "g", minCase = 100, outVar = "y"),
    "skipped: only 10 usable cases"
  )
  expect_null(res_null)
  
  # 3. Test for Success (Pass Logic)
  # We have 10 cases, and we only require 5. Should return a dataframe.
  res_df <- checkMod(df = mock_df, allCols = c("y", "g"), 
                     groupVars = "g", minCase = 5, outVar = "y")
  
  expect_s3_class(res_df, "data.frame")
  expect_equal(nrow(res_df), 10)
})