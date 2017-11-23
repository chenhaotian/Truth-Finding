res$s_aa_n_claims$a_truth <- attributesmapper$a[match(res$s_aa_n_claims$a_truth,attributesmapper$aid)]
res$s_aa_n_claims$a <- attributesmapper$a[match(res$s_aa_n_claims$a,attributesmapper$aid)]

load("~/Downloads/tmpdata")
rawdb <- rawdb[rawdb$a!=4,]

res <- TF(rawdb,"ss",beta  = 1,alpha1 = 1,burnin = 500,maxit = 10000,sample_step = 50,considerpi = TRUE)

xxx <- plyr::ddply(rawdb,"e",function(d){
    tmp <- table(d$a)
    d$mj <- as.integer(names(tmp)[which.max(tmp)])
    return(d)
},.progress = "text")

xxx2 <- res$rawdb_original

xxx2 <- xxx2[order(xxx2$e),]
xxx <- xxx[order(xxx$e),]

##' @title \%IN\%
##' @description data.fram version of %in%.
##' @return a logical vector of length nrow(df1).


##' @title checkprior
##' @description a wrapper function in checking parameters.
##' @return a logical value.
