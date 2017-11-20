## data.fram version of %in%
## returns a logical vector indicating if there is a match or not for each ROW of the left operand.
`%IN%` <- function(df1,df2){
    if( !is(df1,"data.frame") | !is(df2,"data.frame")){
        stop("inputs must be data.frames")
    }
    if(!identical(vapply(df1,class,FUN.VALUE = "a",USE.NAMES = FALSE),vapply(df2,class,FUN.VALUE = "a",USE.NAMES = FALSE))){
        stop("columns of inputs not aligned")
    }
    if(ncol(df1)==1){
        return(df1[[1]]%in%df2[[1]])
    }else{
        names(df1) <- names(df2)
        return(duplicated(rbind(df2,df1))[-(1:nrow(df2))])
    }

}


checkprior <- function(prior,CLASS=character(),NROW=NULL,NCOL=NULL,LENGTH=NULL,ELEMENTCLASSES=NULL,VALUES=NULL){
    ## print(class(prior))
    is(prior,CLASS) &
    ifelse(!is.null(NROW),nrow(prior)%in%NROW,TRUE) &
    ifelse(!is.null(NCOL),ncol(prior)%in%NCOL,TRUE) &
    ifelse(!is.null(LENGTH),length(prior)%in%LENGTH,TRUE) &
    ifelse(!is.null(ELEMENTCLASSES),all(sapply(prior,class)==ELEMENTCLASSES),TRUE) &
    ifelse(!is.null(VALUES),all(prior %in% VALUES), TRUE)
}

library(Rcpp)
sourceCpp("src/TFR.cpp")

## note that in LTM, "recall" and "specificity" of an source are calculated with regard to all the entities that contain claims from this source, this is only meaning full when there're multiple right attributes associated to each entitiy. When there's only one right attribute in each entity, there's no FN and TN in the confusion matrix, only "precision" can be calculated.
## for example:
## Entity e has right attributes (a1,a2,a3), source s1 has one claim a1, source s2 has three claims (a2,a3,a4), then the recall for s1 and s2 are both 1/3, while the specificity for s1 and s2 are 1 and 0, one must use both recall and specificity to measure s1 and s2's qualities. What if e has only one right attribute, which is a1, s1 has claim a1 and s2 has claim a2, then the precision for s1 and s2 are 1 and 0, and only precision alone is able to measure qualities of this two sources.

## beta: Fx2 numeric matrix, each row is the beta prior the corresponding fact
## alpha0: Sx2 numeric matrix, each row is the beta prior for the corresponding source FPR
## alpha1: Sx2 numeric matrix, each row is the beta prior for the corresponding source sensitivity

## model: specify model in use.
##     mn: each entity has MULTIPLE true attributes, attributes NOT sharing among entities.
##     sn: each entity has SINGLE true attribute, attributes NOT sharing among entities.
##     ss: each entity has SINGLE true attribute, attributes SHARING among entities.

## alpha0: , hyper parameter for specificity when model="mn". NULL for model="sn" and "ss"
## alpha1: hyper parameter for recall when model="mn", hyper parameter for precision when model="sn" and "ss"
TF <- function(rawdb,model=c("mn","sn","ss"),beta,alpha0=NULL,alpha1,burnin,maxit,sample_step){

    cat("Checking integrity...")
    burnin <- as.integer(burnin)
    maxit <- as.integer(maxit)
    sample_step<- as.integer(sample_step)
    ## integrity check part1
    if(is.matrix(rawdb)) rawdb <- as.data.frame(rawdb)
    if(!checkprior(rawdb,CLASS = "data.frame",NCOL = 3) | any(duplicated(rawdb))){
        stop("rawdb must be a data.frame or matrix with 3 columns representing entity, attribute and source respectively")
    }else{
        names(rawdb) <- c("e","a","s")
    }
    model <- match.arg(model)
    nentities <- length(unique(rawdb$e))
    nsourcecs <- length(unique(rawdb$s))
    nattributes <- ifelse(model=="ss",length(unique(rawdb$a)),1)
    if((!checkprior(alpha0,CLASS = "matrix",NCOL = 2,NROW=nsourcecs)& model=="mn") |
       !checkprior(alpha1,CLASS = "matrix",NCOL = 2,NROW=nsourcecs))
        stop("alpha1/0 should be a numeric matrix of dimension number_of_sources x 2")
    if((!checkprior(beta,CLASS = "numeric",LENGTH = 2) & model=="mn") |
       (!checkprior(beta,CLASS = "numeric",LENGTH = 1) & model=="sn") |
       (!checkprior(beta,CLASS = "numeric",LENGTH = nattributes) & model=="ss"))
        stop("beta should be a length 2 numberic when model = mn, a length 1 numeric when model = sn, a lengh number_of_attributes vector when model = ss")
    cat("done.\n")



    sourcesmapper <- unique(rawdb[,"s",drop=FALSE])
    sourcesmapper$sid <- 0L:(nsourcecs-1L)
    if(model=="mn"){
        cat("Preparing mappers...")
        factsmapper <- unique(rawdb[,c("e","a")])
        factsmapper$fid <- 0L:(nrow(factsmapper)-1L)
        cat("done.\n")

        cat("Preparing claims...")
        claims <- do.call(rbind.data.frame,lapply(split(rawdb,rawdb$e),function(d){
            tmp <- expand.grid(unique(d$a),unique(d$s),stringsAsFactors = FALSE)
            names(tmp) <- c("a","s")
            tmp$o <- as.integer(tmp %IN% d[,c("a","s")])
            tmp$e <- d$e[1]
            tmp
        }))
        row.names(claims) <- NULL
        claims <- merge(merge(claims,factsmapper,all.x = TRUE,sort = FALSE)[c("fid","s","o")],sourcesmapper,all.x = TRUE,sort = FALSE)[c("fid","sid","o")]
        claims <- claims[order(claims$fid,claims$sid,claims$o),]
        cat("done.\n")

        cat("Initializing truth & truth-source-claim count...")
        facts <- data.frame(fid=factsmapper$fid,t=sample(c(0L,1L),nrow(factsmapper),replace = TRUE),stringsAsFactors = FALSE)
        facts <- facts[order(facts$fid),]
        ## count, truth-source-claim
        ## raw sensitivity and FPR data
        ctsc <- expand.grid(c(0L,1L),sourcesmapper$sid,c(0L,1L),stringsAsFactors = FALSE)
        names(ctsc) <- c("t","sid","o")
        cwit <- merge(claims,facts,all.x = TRUE,sort = FALSE) # claim with initial truth, a tmp variable in calculating ctsc
        ctsc <- merge(ctsc,aggregate(fid~sid+t+o,data = cwit,FUN = length),all.x = TRUE,sort = FALSE)
        rm(cwit)
        names(ctsc)[4] <- "count"
        ctsc$count[is.na(ctsc$count)] <- 0L
        ctsc <- ctsc[order(ctsc$t,ctsc$sid,ctsc$o),]
        ## fact-claim index, claim start posiitons of each fact
        ## to simplify C++ code, here the index starts with zero
        fcidx <- match(facts$fid,claims$fid)-1L             #facts and claims must be sorted beforehand
        fcidx <- c(fcidx,nrow(claims))                      #append end position for easier representation
        cat("done.\n")

        cat("Sampling...\n")
        res <- truthfinding_mn(facts=facts$t,fcidx=fcidx, claims=as.matrix(claims), ctsc=as.matrix(ctsc), beta=beta, alpha0=alpha0, alpha1=alpha1,burnin = burnin, maxit = maxit,sample_step = sample_step)
        cat("\n all done \n")
        
    }
    else if(model=="sn"){
        cat("Preparing mappers...")
        entities <- unique(rawdb[,"e",drop=FALSE])
        entities$eid <- 0L:(nentities-1L)
        rawdb <- rawdb[order(rawdb$e),] #group by e, in order to split
        attributesmapper <- rawdb[,c("e","a")]
        attributesmapper$aid <- unlist(sapply(split(rawdb$a,rawdb$e),function(a){
            match(a,unique(a))-1L
        },simplify = FALSE,USE.NAMES = FALSE),use.names = FALSE)
        ## tokenize rawdb
        rawdb_original <- rawdb         #save a backup
        rawdb$a <- attributesmapper$aid
        rawdb$e <- entities$eid[match(rawdb$e,entities$e)]
        rawdb$s <- sourcesmapper$sid[match(rawdb$s,sourcesmapper$s)]
        rawdb <- rawdb[order(rawdb$e,rawdb$a,rawdb$s),]  # !!! important!!!
        attributesmapper <- unique(attributesmapper)
        cat("done.\n")
        
        cat("Initializing counts & turths ...")
        e_n_attr <- aggregate(a~e,data = rawdb,function(l){length(unique(l))})
        e_n_attr <- e_n_attr$a[order(e_n_attr$e)] #claimed attributes count for each entity
        s_n_claims <- aggregate(e~s,data = rawdb,length)
        s_n_claims <- s_n_claims$e[order(s_n_claims$s)] #claims of each source(some claims may be too few)
        max_nattributes <- max(e_n_attr)
        e_truths <- sapply(e_n_attr,function(l){ #initialize truth for each entity
            one_cat_zero_begin(rep(1,l+1)) #?????? l+1 or l ???????
        },simplify = TRUE,USE.NAMES = FALSE)
        s_n_right <- sapply(split(rawdb$a==e_truths[match(rawdb$e,entities$eid)],rawdb$s),sum,simplify = TRUE,USE.NAMES = TRUE)
        s_n_right <- s_n_right[order(as.integer(names(s_n_right)))] #number of right cliams for each source
        ecidx <- c(match(entities$eid,rawdb$e)-1L,nrow(rawdb))
        cat("done.\n")

        cat("Sampling...\n")
        res <- truthfinding_sn(ecidx = ecidx,e_truths = e_truths,e_n_attr = e_n_attr,s_n_right = s_n_right, s_n_claims = s_n_claims,rawdb = as.matrix(rawdb),beta = beta, alpha1 = alpha1, max_nattributes = max_nattributes, burnin = burnin, maxit = maxit,sample_step = sample_step)
        cat("\n all done \n")
    }


    ## res$ctsc <- ctsc
    ## res$claims <- claims
    ## res$factsmapper <- factsmapper
    ## res$sourcesmapper <- sourcesmapper
    return(res)
}


beta=1;alpha0=1;alpha1=1;burnin=100;maxit=1000;sample_step=30

testdata <- rawdb[rawdb$e %in% sample(rawdb$e,10000),]

## res <- TF_binary(rawdb = testdata,binary = TRUE,beta = 1,alpha0 = matrix(c(8,2),nrow = length(unique(testdata$s)),ncol = 2,byrow = TRUE),alpha1 = 0.1,burnin = 500,maxit = 2000,sample_step = 10)

load("~/Downloads/crowd_task_am_history")
rawdb <- crowd[,c("tracking_id","decision","user_id")]
names(rawdb) <- c("e","a","s")
rawdb <- rawdb[rawdb$e %in% sample(rawdb$e,20000),]

xxx <- aggregate(e~s,data = rawdb,length)
rawdb <- rawdb[rawdb$s%in%xxx$s[xxx$e > quantile(xxx$e,0.2)],]



res <- TF(rawdb = rawdb,model = "sn",beta = 1,alpha1 = matrix(c(8,2),nrow = length(unique(rawdb$s)),ncol = 2,byrow = TRUE),burnin = 1000,maxit = 3000,sample_step = 100)

res <- TF(rawdb = rawdb,model = "sn",beta = 1,alpha1 = matrix(1,nrow = length(unique(rawdb$s)),ncol = 2),burnin = 1000,maxit = 10000,sample_step = 100)

library(reshape2)
spa <- read.table("~/Downloads/consolidation_spa.csv",header = TRUE,sep = ",",quote = "\"",na.strings = "NA",stringsAsFactors = FALSE)
spa <- spa[,-1]
spa <- na.omit(melt(spa,id.vars="item_id"))

