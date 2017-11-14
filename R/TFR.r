## data
rawdb <- unique(data.frame(e=sample(letters[1:4],30,replace = TRUE),a=sample(c("spa","breakfast","iron"),size = 30,replace = TRUE),s=sample(toupper(letters[1:3]),30,replace = TRUE),stringsAsFactors = FALSE))

if(any(duplicated(rawdb))) stop("duplicated entries in input data")

## facts
facts <- unique(rawdb[,c("e","a")])
facts <- facts[order(facts$e,facts$a),]

## claims

factsource <- unique(rawdb[,c("e","s")])

lapply(split(facts,1:nrow(facts)),function(d){
    
})

lapply(split(rawdb,rawdb$e),function(d){
    claims <- expand.grid(d$a,d$s,stringsAsFactors = FALSE)
    claims$o <- claims %IN% d[,c("a","s")]
    claims
})

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

a br
b sp

a br
a sp
b br
b sp
