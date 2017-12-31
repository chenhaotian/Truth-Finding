load("~/Downloads/tmpdata")
rawdb <- rawdb[rawdb$a!=4,]

library(Rcpp)
source("R/TFR.r")
sourceCpp("src/TFR.cpp")

res <- TF(rawdb,"ss",beta  = 1,alpha1 = 1,burnin = 500,maxit = 5000,sample_step = 50)

xxx <- plyr::ddply(rawdb,"e",function(d){
    tmp <- table(d$a)
    d$mj <- as.integer(names(tmp)[which.max(tmp)])
    return(d)
},.progress = "text")

table(unique(xxx[,c("e","mj")])$mj)
res$pi

head(res$s_aa_n_claims,30)

xxx2 <- res$rawdb_original

xxx2 <- xxx2[order(xxx2$e),]
xxx <- xxx[order(xxx$e),]
tmp <- data.frame(t=xxx2$e_truth,m=xxx$mj)


## precision
tmp2 <- plyr::ddply(res$s_aa_n_claims,"s",function(d){
    data.frame(s=d$s[1],rate=sum(d$count[d$a_truth==d$a]) / sum(d$count),cnt=sum(d$count))
})

tmp2 <- tmp2[order(tmp2$cnt),]

plot(tmp2$cnt,tmp2$rate,type = "l")


table(res$rawdb_original$e_truth[res$rawdb_original$e_truth_prob<0.9])/length(res$rawdb_original$e_truth[res$rawdb_original$e_truth_prob<0.9])

table(res$rawdb_original$e_truth[res$rawdb_original$e_truth_prob>=0.9])/length(res$rawdb_original$e_truth[res$rawdb_original$e_truth_prob>=0.9])

93.7
97.5
75.0
78.0
