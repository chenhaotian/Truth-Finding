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
rawdb <- unique(data.frame(e=sample(letters[1:4],30,replace = TRUE),a=sample(c("spa","breakfast","iron"),size = 30,replace = TRUE),s=sample(toupper(letters[1:3]),30,replace = TRUE),stringsAsFactors = FALSE))

if(any(duplicated(rawdb))) stop("duplicated entries in input data")

## mappers
factsmapper <- unique(rawdb[,c("e","a")])
factsmapper$fid <- 1:nrow(factsmapper)

sourcesmapper <- unique(rawdb[,c("s")])
sourcesmapper <- data.frame(s=sourcesmapper,sid=1:length(sourcesmapper),stringsAsFactors = FALSE)

## claims
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

## initialize truth
facts <- data.frame(fid=factsmapper$fid,t=sample(c(0,1),nrow(factsmapper),replace = TRUE),stringsAsFactors = FALSE)
facts <- facts[order(facts$fid),]

## count, truth-source-claim
## raw sensitivity and FPR data
ctsc <- expand.grid(c(0,1),sourcesmapper$sid,c(0,1),stringsAsFactors = FALSE)
names(ctsc) <- c("t","sid","o")
cwit <- merge(claims,facts,all.x = TRUE,sort = FALSE) # claim with initial truth, a tmp variable in calculating ctsc
ctsc <- merge(ctsc,aggregate(fid~sid+t+o,data = cwit,FUN = length),all.x = TRUE,sort = FALSE)
rm(cwit)
names(ctsc)[4] <- "count"
ctsc$count[is.na(ctsc$count)] <- 0L
ctsc <- ctsc[order(ctsc$t,ctsc$sid,ctsc$o),]

## fact-claim index, claim start posiitons of each fact
## to simplify C++ code, here the index starts with zero
fcidx <- match(facts$fid,claims$fid)-1             #facts and claims must be sorted beforehand
fcidx <- c(fcidx,nrow(claims))                     #append end position for easier representation


input: facts claims fcidx ctsc

for( f in fid){
    
    cf <- claims[claims$fid==f,]
    ## get startidx and endidx of f in claim
    
    for(c in cf){
        calculate p with ctsc
    }

    update facts
    for(c in cf){
        update ctsc
    }
}



## claims <- merge(claims,factsmapper,by = c("e","a"),all.x = TRUE)

## claims$t <- facts$t[claims$fid]

## s <- unique(rawdb$s)
## t <- c(TRUE,FALSE)
## o <- c(TRUE,FALSE)




## (fid)-(source,claim) mapper
fsc <- subset(claims,select=c(fid,s,o))



## count of source-attribute-claim
