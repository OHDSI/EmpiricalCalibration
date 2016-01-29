This is the first submission of this package

---

## Test environments
* Ubuntu 14.04.3 LTS (Travis), R 3.2.3
* Windows 7, R 3.2.3

## R CMD check results

There were no ERRORs or WARNINGs. I see 3 NOTES (on Travis):

* CRAN New submission

* No repository set, so cyclic dependency check skipped 
  (Not sure how to fix this on Travis. Does not appear on Windows where I've set the repository)
  
* Possible R code problems: These are all due to non-standard evaluations not detected by devtools

## Downstream dependencies

This is a new submission, there are no downstream dependencies.