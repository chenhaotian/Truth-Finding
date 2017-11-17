## data.fram version of %in%
## returns a logical vector indicating if there is a match or not for EACH ROW of its left operand.
`%IN%` <- function(df1,df2){
    if( !is(df1,"data.frame") | !is(df2,"data.frame")){
        stop("inputs must be data.frames")
    }
    if(!identical(vapply(df1,class,FUN.VALUE = "a",USE.NAMES = FALSE),vapply(df2,class,FUN.VALUE = "a",USE.NAMES = FALSE))){
        stop("columns of inputs not aligned")
    }
    if(ncol(df1)==1){
        return(df1[[1]]%in%df2[[1]])
    }
    res <- vapply(1:ncol(df1),function(i){
        match(df1[[i]],df2[[i]], nomatch = 0L)
    },FUN.VALUE = rep(0L,nrow(df1)),USE.NAMES = FALSE)

    return(apply(res,MARGIN = 1,function(x){
        if(x[1]!=0L){
            all(x[-1]==x[1])
        }else{
            FALSE
        }
    }))
}

## label:
## e: entity
## a: attribute
## s: source
## o: claim
## t: truth

## raw data
rawdb <- unique(data.frame(e=sample(letters,3000,replace = TRUE),a=sample(c("spa","breakfast","iron"),size = 3000,replace = TRUE),s=sample(toupper(letters[1:3]),300,replace = TRUE),stringsAsFactors = FALSE))

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

## beta: Fx2 numeric matrix, each row is the beta prior the corresponding fact
## alpha0: Sx2 numeric matrix, each row is the beta prior for the corresponding source FPR
## alpha1: Sx2 numeric matrix, each row is the beta prior for the corresponding source sensitivity
TF_binary <- function(rawdb,binary=FALSE,beta,alpha0,alpha1,burnin,maxit,sample_step){

    cat("Preparing mappers...")
    burnin <- as.integer(burnin)
    maxit <- as.integer(maxit)
    sample_step<- as.integer(sample_step)
    ## integrity check part1
    if(is.matrix(rawdb)) rawdb <- as.data.frame(rawdb)
    if(!checkprior(rawdb,CLASS = "data.frame",NCOL=c(3,4)) | any(duplicated(rawdb))) stop("rawdb must be a data.frame or matrix with 3 or 4 columns and no duplicate entries")
    if(ncol(rawdb)==3){
        names(rawdb) <- c("e","a","s")
    }else{                                 #ncol(rawdb)==4)
        names(rawdb) <- c("e","a","s","o") #binary claim included
        if(!checkprior(rawdb$o,CLASS = "integer",VALUES=c(0L,1L))) stop("the forth column of rawdb must be integer of 0s and 1s")
    }

    ## 1. mappers
    factsmapper <- unique(rawdb[,c("e","a")])
    factsmapper$fid <- 0L:(nrow(factsmapper)-1L)
    sourcesmapper <- unique(rawdb[,c("s")])
    sourcesmapper <- data.frame(s=sourcesmapper,sid=0L:(length(sourcesmapper)-1L),stringsAsFactors = FALSE)


    ## integrity check part2
    if(is.numeric(beta)) beta <- matrix(beta,nrow=nrow(factsmapper),ncol = 2)
    if(is.numeric(alpha0)) alpha0 <- matrix(alpha0,nrow=nrow(sourcesmapper),ncol = 2)
    if(is.numeric(alpha1)) alpha1 <- matrix(alpha1,nrow=nrow(sourcesmapper),ncol = 2)
    if(!checkprior(beta,"matrix",nrow(factsmapper),2)) stop("beta should be of length 1 or a numeric matrix of dimension n.facts x 2")
    if(!checkprior(alpha0,"matrix",nrow(sourcesmapper),2) |
       !checkprior(alpha1,"matrix",nrow(sourcesmapper),2)) stop("alpha0 or alpha1 should both be of length 1 or a numeric matrix of dimension n.sources x 2")
    cat("done.\n")

    cat("Preparing claims...")
    ## 2. claims
    if(ncol(rawdb)==3){
        claims <- do.call(rbind.data.frame,lapply(split(rawdb,rawdb$e),function(d){
            tmp <- expand.grid(unique(d$a),unique(d$s),stringsAsFactors = FALSE)
            names(tmp) <- c("a","s")
            tmp$o <- as.integer(tmp %IN% d[,c("a","s")])
            tmp$e <- d$e[1]
            tmp
        }))
    }else{                              #ncol(rawdb)==4
        claims <- rawdb
    }
    row.names(claims) <- NULL
    claims <- merge(merge(claims,factsmapper,all.x = TRUE,sort = FALSE)[c("fid","s","o")],sourcesmapper,all.x = TRUE,sort = FALSE)[c("fid","sid","o")]
    claims <- claims[order(claims$fid,claims$sid,claims$o),]
    cat("done.\n")

    cat("Initializing truth & truth-source-claim count...")
    ## 3. initialize truth
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
    res <- truthfinding_binary(facts=facts$t,fcidx=fcidx, claims=as.matrix(claims), ctsc=as.matrix(ctsc), beta=beta, alpha0=alpha0, alpha1=alpha1,burnin = burnin, maxit = maxit,sample_step = sample_step)
    cat("\n all done \n")

    res$ctsc <- ctsc
    res$claims <- claims
    res$factsmapper <- factsmapper
    res$sourcesmapper <- sourcesmapper
    return(res)
}

## match <- crowd[crowd$decision==1,]

load("~/Downloads/crowd_task_am_history")

ones <- crowd[,c("tracking_id","decision","user_id","decision")]; names(ones) <- c("e","a","s","o")
ones$o <- 1L
zeros <- crowd[,c("tracking_id","decision","user_id","qa_decision")]; names(zeros) <- c("e","a","s","o")

rawdb <- rbind(ones,zeros)
rawdb <- unique(rawdb)

beta=1;alpha0=1;alpha1=1;burnin=100;maxit=1000;sample_step=30

res <- TF_binary(rawdb = rawdb,binary = TRUE,beta = 1,alpha0 = matrix(c(8,2),nrow = 262,ncol = 2,byrow = TRUE),alpha1 = 0.1,burnin = 500,maxit = 2000,sample_step = 10)






## crowd <- merge(crowd,x,by.x = "external_task_id",by.y = "task_id",all.x = FALSE, all.y = FALSE)

