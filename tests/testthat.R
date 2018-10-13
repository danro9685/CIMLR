Sys.setenv("R_TESTS" = "")

library("testthat")
library("CIMLR")

test_check("CIMLR")
