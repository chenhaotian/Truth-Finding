##' @title %IN%
##' @description data.fram version of %in%.
##' @return a logical vector of length nrow(df1).
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

##' @title checkprior
##' @description a wrapper function in checking parameters.
##' @return a logical value.
checkprior <- function(prior,CLASS=character(),NROW=NULL,NCOL=NULL,LENGTH=NULL,ELEMENTCLASSES=NULL,VALUES=NULL){
    ## print(class(prior))
    is(prior,CLASS) &
    ifelse(!is.null(NROW),nrow(prior)%in%NROW,TRUE) &
    ifelse(!is.null(NCOL),ncol(prior)%in%NCOL,TRUE) &
    ifelse(!is.null(LENGTH),length(prior)%in%LENGTH,TRUE) &
    ifelse(!is.null(ELEMENTCLASSES),all(sapply(prior,class)==ELEMENTCLASSES),TRUE) &
    ifelse(!is.null(VALUES),all(prior %in% VALUES), TRUE)
}

##' @title TF
##' @description Truth finding algorithms with Bayes nets.
##' @details return a data.frame of current queuing orders, each row of
##' the data.frame representing an order, queryorder will return all of the
##' queuing orders if orderid is NULL. when there is no queuing orders,
##' queryorder will return a data.frame with 0 rows.
##' @param rawdb data.frame with three columns, representing "entity", "attribute",
##' "information source" respectively.
##' @param model character, specifying the model to be used
##' "ss", each entity can only have one true attribute, attribute labels are sharing
##' among entities.
##' "sn", each entity can only have one true attribute, attribute label in each entity can be unique.
##' "mn": each entity can have multiple attributes.
##' @param beta numeric, hyper parameter for attribute distribution, default 1
##' @param alpha0 numeric, hyper parameter for false positive rate in model "mn"
##' @param alpha1 numeric, hyper parameter for source recall rate.
##' @param burnin integer, number of burn-in samples.
##' @param maxit integer, number of iterations.
##' @param sample_step integer, sample step length in posterior estimation.
##' @param considerpi logical, 
##' @return a list
##' @examples
##' \dontrun{
##'   #nothing here yet
##' }
##' @export
TF <- function(rawdb,model=c("ss","sn","mn"),beta=1,alpha0=NULL,alpha1=1,burnin=500L,maxit=3000L,sample_step=100L,considerpi=TRUE){

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
    nsources <- length(unique(rawdb$s))
    nattributes <- length(unique(rawdb$a))
    if(nattributes<2 & model=="ss") stop("there need to be at least two unique attributes in rawdb")
    ## if(!checkprior(alpha0,CLASS = "matrix",NCOL = 2L,NROW=nsources))
    ##     stop("alpha1/0 should be a numeric matrix of dimension number_of_sources x 2")
    ## if((!checkprior(beta,CLASS = "numeric",LENGTH = 2L) & model=="mn") |
    ##    (!checkprior(beta,CLASS = "numeric",LENGTH = 1L) & model=="sn") |
    ##    (!checkprior(beta,CLASS = "numeric",LENGTH = nattributes) & model=="ss"))
    ##     stop("beta should be a length 2 numberic when model = mn, a length 1 numeric when model = sn, a lengh number_of_attributes vector when model = ss")
    ## cat("done.\n")

    sourcesmapper <- unique(rawdb[,"s",drop=FALSE])
    sourcesmapper$sid <- 0L:(nsources-1L)
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
        gc(verbose=FALSE)
        names(ctsc)[4] <- "count"
        ctsc$count[is.na(ctsc$count)] <- 0L
        ctsc <- ctsc[order(ctsc$t,ctsc$sid,ctsc$o),]
        ## fact-claim index, claim start posiitons of each fact
        ## to simplify C++ code, here the index starts with zero
        fcidx <- match(facts$fid,claims$fid)-1L             #facts and claims must be sorted beforehand
        fcidx <- c(fcidx,nrow(claims))                      #append end position for easier representation
        cat("done.\n")

        claim <- as.matrix(claims)
        ctsc <- as.matrix(ctsc)
        gc(verbose = FALSE)
        
        cat("Sampling...\n")
        res <- truthfinding_mn(facts=facts$t,fcidx=fcidx, claims=claims, ctsc=ctsc, beta=beta, alpha0=alpha0, alpha1=alpha1,burnin = burnin, maxit = maxit,sample_step = sample_step)
        cat("\n all done \n")
        
    }
    else if(model%in%c("sn","ss")){
        cat("Preparing mappers...")
        entities <- unique(rawdb[,"e",drop=FALSE])
        entities$eid <- 0L:(nentities-1L)
        rawdb_original <- rawdb         #save a backup
        ## tokenize attributes
        if(model=="sn"){
            rawdb <- rawdb[order(rawdb$e),] #group by e, in order to split
            attributesmapper <- rawdb[,c("e","a")]
            attributesmapper$aid <- unlist(sapply(split(rawdb$a,rawdb$e),function(a){
                match(a,unique(a))-1L
            },simplify = FALSE,USE.NAMES = FALSE),use.names = FALSE)
            rawdb$a <- attributesmapper$aid
            attributesmapper <- unique(attributesmapper)
        }
        else if(model=="ss"){
            attributesmapper <- unique(rawdb[,"a",drop=FALSE])
            attributesmapper$aid <- 0L:(nattributes-1L)
            rawdb$a <- attributesmapper$aid[match(rawdb$a,attributesmapper$a)]
        }
        
        rawdb$e <- entities$eid[match(rawdb$e,entities$e)]
        rawdb$s <- sourcesmapper$sid[match(rawdb$s,sourcesmapper$s)]
        rawdb <- rawdb[order(rawdb$e,rawdb$a,rawdb$s),]  # !!! important!!!
        cat("done.\n")
        
        cat("Initializing counts & turths ...")
        ecidx <- c(match(entities$eid,rawdb$e)-1L,nrow(rawdb))
        if(model=="sn"){
            e_n_attr <- aggregate(a~e,data = rawdb,FUN = function(l){length(unique(l))})
            e_n_attr <- e_n_attr$a[order(e_n_attr$e)] #claimed attributes count for each entity
            s_n_claims <- aggregate(e~s,data = rawdb,FUN = length)
            s_n_claims <- s_n_claims$e[order(s_n_claims$s)] #claims of each source(some claims may be too few)
            max_nattributes <- max(e_n_attr)
            if(max_nattributes<2) stop("rawdb: there need to be at least one entity with at least two different attributes")
            e_truths <- sapply(e_n_attr,function(l){ #initialize truth for each entity
                one_cat_zero_begin(rep(1,l+1L)) #?????? l+1 or l ???????
            },simplify = TRUE,USE.NAMES = FALSE)
            s_n_right <- sapply(split(rawdb$a==e_truths[match(rawdb$e,entities$eid)],rawdb$s),sum,simplify = TRUE,USE.NAMES = TRUE)
            s_n_right <- s_n_right[order(as.integer(names(s_n_right)))] #number of right cliams for each source
            cat("done.\n")
            cat("Sampling...\n")

            rawdb <- as.matrix(rawdb)
            
            res <- truthfinding_sn(ecidx = ecidx,e_truths = e_truths,e_n_attr = e_n_attr,s_n_right = s_n_right, s_n_claims = s_n_claims,rawdb = rawdb,beta = beta, alpha1 = alpha1, max_nattributes = max_nattributes, burnin = burnin, maxit = maxit,sample_step = sample_step)
        }
        else if(model=="ss"){
            ## random initial value will lead to varies local maxima, here use majority rules instead.
            rawdb$a_truth <- unname(do.call(c,lapply(split(rawdb$a,rawdb$e),function(l){
                tmp <- table(l)
                rep(as.integer(names(tmp)[which.max(tmp)]),length(l))
            })))

            e_truths <- rawdb$a_truth[!duplicated(rawdb$e)]

            
            ## nsouce x nattributes matrix
            s_a_n_claims <- expand.grid(0L:(nsources-1L),0L:(nattributes-1L))
            names(s_a_n_claims) <- c("s","a_truth")
            s_a_n_claims <- merge(s_a_n_claims,aggregate(e~s+a_truth,data = rawdb,FUN = length),all.x = TRUE,sort = FALSE)
            s_a_n_claims$e[is.na(s_a_n_claims$e)] <- 0L
            names(s_a_n_claims)[3] <- "count"
            s_a_n_claims <- s_a_n_claims[order(s_a_n_claims$s,s_a_n_claims$a_truth),] #claims of each source(some claims may be too few)

            ## nsource x nattributes x nattributes array
            s_aa_n_claims <- expand.grid(0L:(nsources-1L),0L:(nattributes-1L),0L:(nattributes-1L))
            names(s_aa_n_claims) <- c("s","a_truth","a")
            s_aa_n_claims <- merge(s_aa_n_claims,aggregate(e~s+a_truth+a,data = rawdb,FUN = length),all.x = TRUE,sort = FALSE)
            s_aa_n_claims$e[is.na(s_aa_n_claims$e)] <- 0L
            names(s_aa_n_claims)[4] <- "count"
            s_aa_n_claims <- s_aa_n_claims[order(s_aa_n_claims$s,s_aa_n_claims$a_truth,s_aa_n_claims$a),] #claims of each source(some claims may be too few)

            pi <- as.integer(table(e_truths)[as.character(0L:(nattributes-1))])
            pi[is.na(pi)] <- 0
            ## print(pi)

            rawdb <- as.matrix(rawdb[,-4])
            
            res <- truthfinding_ss_fullpar(ecidx = ecidx,e_truths = e_truths, s_aa_n_claims = s_aa_n_claims$count, s_a_n_claims = s_a_n_claims$count,rawdb = rawdb,beta = beta, pi = pi,alpha1 = alpha1, nattributes = nattributes, nsources = nsources, nentities = nentities, burnin = burnin, maxit = maxit,sample_step = sample_step,considerpi = considerpi)
        }
        cat("\n all done \n")
    }

    if(model=="ss"){
        rawdb_original$e_truth <- res$e_truths[match(rawdb_original$e,entities$e)]
        rawdb_original$e_truth <- attributesmapper$a[match(rawdb_original$e_truth,attributesmapper$aid)]
        res$rawdb_original <- rawdb_original
        names(res$pi) <- as.character(attributesmapper$a)
        
        s_aa_n_claims$count <- res$s_aa_n_claims
        s_aa_n_claims$a_truth <- attributesmapper$a[match(s_aa_n_claims$a_truth,attributesmapper$aid)]
        s_aa_n_claims$a <- attributesmapper$a[match(s_aa_n_claims$a,attributesmapper$aid)]
        res$s_aa_n_claims <- s_aa_n_claims
    }
    else if(model=="mn"){
        res$ctsc <- ctsc
        res$claims <- claims
        res$factsmapper <- factsmapper
        res$sourcesmapper <- sourcesmapper
    }
    
    return(res)
}


