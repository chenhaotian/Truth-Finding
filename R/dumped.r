load("~/Downloads/tmpdata")
rawdb <- rawdb[rawdb$a!=4,]

library(Rcpp)
source("R/TFR.r")
sourceCpp("src/TFR.cpp")

res <- TF(rawdb,"ss",beta  = 1,alpha1 = 1,burnin = 500,maxit = 10000,sample_step = 50)

xxx <- plyr::ddply(rawdb,"e",function(d){
    tmp <- table(d$a)
    d$mj <- as.integer(names(tmp)[which.max(tmp)])
    return(d)
},.progress = "text")

xxx2 <- res$rawdb_original

xxx2 <- xxx2[order(xxx2$e),]
xxx <- xxx[order(xxx$e),]
tmp <- data.frame(t=xxx2$e_truth,m=xxx$mj)


## precision
tmp2 <- plyr::ddply(res$s_aa_n_claims,"s",function(d){
    data.frame(s=d$s[1],rate=sum(d$count[d$a_truth==d$a]) / sum(d$count),cnt=sum(d$count))
})

tmp2 <- tmp2[order(tmp2$cnt),]

plot(tmp2$rate,type = "l")
