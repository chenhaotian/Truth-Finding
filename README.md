# Truth-Finding
Find hidden truth from conflicting information sources.

## Illustrations

=="Website A claims that person B has a tank parking in the garage."==

Website A is an "information source", person B is an "entity" and "has a tank parking in the garage" is an "attribute", they together form a claim:

| Entity        | Attribute     | Source  |
| ------------- |:-------------:| -----:|
| Person B      | Has a tank paring in the garage | Website A |

In a typical truth finding scenarios, there may be multiple information sources providing claims about multiple entities having some kind of attributes, and there may exist conflicts between these claims. The model is to help finding the hidden truths, meanwhile calculate the credibilities of information sources.

Generally there are three types of truth finding scenarios:
+ **S**ingle true attribute and **S**haring attributes(**SS**): `M sources claiming that N people's gender are either male or not male` , in this example all sources' claims are about the same set of attributes, which is {"male", "mot male"}, that means **attribute labels are sharing among entities,** moreover, each person can not be "male" and "not male" at the same time, suggesting that **each entity can only have one true attribute.** 
+ **S**ingle true attribute and attributes **N**ot sharing among entities(**SN**): `M sources claiming the one most favourite things of N people` , here the "one most favourite thing" can be anything, so **attribute labels are not sharing among entities,** meanwhile each person has only one most favourite thing, i.e. **each entity can only have one true attribute.**
+ **M**ultiple true attributes and attributes **N**ot sharing among entities(**MN**): `M sources claiming the casts of N movies` , casts are combinations of names, the number of unique names can be infinitely large, and some actors or actresses may only appear in one movie, so **attribute labels are not sharing among entities,** while each movie can have multiple actors and actresses, indicating **each entity can have multiple attributes.**  [ [2](https://arxiv.org/pdf/1203.0058.pdf) ]


## Usage 
Load the function `TF()` directly from source: 
```R
library(Rcpp) #for loading cpp functions
sourceCpp(code=paste(readLines("https://raw.githubusercontent.com/chenhaotian/Truth-Finding/master/src/TFR.cpp"),collapse = "\n"))
source("https://raw.githubusercontent.com/chenhaotian/Truth-Finding/master/R/TFR.r") 
``` 

Parameters: 
```R 
args(TF) 
``` 
+ **rawdb** : data.frame, must satisfy ncol(rawdb)=3, column names can be arbitary but the content of each column must represent "entity", "attribute" and "information source" respectively.
+ **model** : character, specifying model to be used, should be one of "ss","sn" and "mn".
	+ ss: **each entity can only have one true attribute, attribute labels are sharing among entities.**
	+ sn: **each entity can only have one true attribute, attribute label in each entity can be unique.**
	+ mn: **each entity can have multiple attributes.**

Example:
```R
## read your own rawdb as a three column data.frame, then run the following code to find latent truths with model="ss":
res <- TF(rawdb,model="ss")
```

## Reference
1. [Machine learning: a probabilistic perspective](http://cds.cern.ch/record/1981503)
2. [A Bayesian Approach to Discovering Truth from Conflicting Sources for Data Integration](https://arxiv.org/pdf/1203.0058.pdf)
