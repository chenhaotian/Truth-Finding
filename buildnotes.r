startingPATH <- "."

Rcpp::Rcpp.package.skeleton(name = "truthr",path = startingPATH)
packagePATH <- paste0(startingPATH,"/truthr")
## copy .r and .cpp files into /truthr/R/ and truthr/src/ respectively
Rcpp::compileAttributes(packagePATH,verbose = TRUE)
devtools::document(packagePATH)
devtools::check(packagePATH)
## fix warning and errors according to check() results

devtools::build(packagePATH)

## remove Read-and-delete-me
