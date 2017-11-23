# Truth-Finding
Find hidden truth from conflicting information sources.

## Illustrations
"Website A claims that person B has a tank parking in the garage."
Webside A is an "information source", person B is an "entity" and "has a tank parking in the garage" is an "attribute", they together form a claim:

| Entity        | Attribute     | Source  |
| ------------- |:-------------:| -----:|
| Person B      | Has a tank paring in the garage | Website A |

In a typical truth finding scenarios, there are be multiple information providers claiming that multiple entities have some kind of attributes, and among the claims there may exist conflicts. 

Generally there are three types of truth finding problems.
+ `Some sources claiming that a group of people's gender are either male or not male` , here all soures' claimings are about the same set of attributes {"male", "mot male"}, that means **attribute labels are sharing among entities,** and each person can not have "male" and "not male" attributes at the same time, that means **each entity can only have one true attribute.**
+ `Some sources claiming the one most favourite thing of a group of people` , here the attributes are not fixed, it can be anything, so **attribute labels are not sharing among entities,** meanwhile each person has only one most favourite thing, so **each entity can only have one true attribute.**
+ `Some sources claiming the cast of some of movies` , cast can be any combination of any person's name, the size of the attribute set an unlimited, so **attribute labels are not sharing among entities,** while each movie can have multiple actors and actresses, so **each entity can have multiple attributes.**  (ref[2](https://arxiv.org/pdf/1203.0058.pdf) )

The model provided here is to help find the hidden truth in all three types of scenarios, meanwhile calculate the credibilities of the information sources.

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
+ **model** : character, specifying model to be used, should be one of c("ss","sn","mn") 
	+ ss: each entity can only have one true attribute, attribute labels are sharing among entities.
	+ sn: each entity can only have one true attribute, attribute label in each entity can be unique.
	+ mn: each entity can have multiple attributes.


## Reference
1. [Machine learning: a probabilistic perspective](http://cds.cern.ch/record/1981503)
2. [A Bayesian Approach to Discovering Truth from Conflicting Sources for Data Integration](https://arxiv.org/pdf/1203.0058.pdf)
