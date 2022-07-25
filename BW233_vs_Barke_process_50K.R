# Author Sybille U Mittmann
# 50K iSelect BW233 vs Barke data processing

library(gdata)
library(stringr)

mydata <- read.xls(file.choose(), row.names = 1)
mydata[mydata == ""] <- NA
mydata[mydata == "NA"] <- NA

## don't use : mydata<- data.frame(lapply(mydata, as.character), stringsAsFactors=FALSE) ##

## Transform factors into characters: ##
mydata = data.frame(mydata,stringsAsFactors = F)

## REMOVE MONOMORPHIC MARKERS ##
mydata2 <- mydata[, apply(mydata, 2,function(x) length(unique(na.omit(x)))) > 1]

## Transpose data frame ##
mydata3 <- t(mydata2)

## Remove remaining monomorphic markers ##
mydata4<-mydata3[!mydata3[,1]==mydata3[,2],]

## Transpose back ##
mydata5 <- t(mydata4)

## Turn matrix back into dataframe ##
mydata5 = data.frame(mydata5,stringsAsFactors = F)

##Remove empty columns ##
mydata6 <- mydata5[,!apply(mydata5, 2, function(x) all(gsub(" ", "", x)=="", na.rm=TRUE))]

## Change Names of genotypes ##
row.names(mydata6)<-str_sub(row.names(mydata6), start = -14, end = -1)
View(mydata6)

## CHANGE ALLELE NAMES TO MAPPING FORMAT ##
# Allele names can be changed referred to each of the parental alleles. For example if the allele is the same as the parental "mutant" we can call that allele "m", if it´s like the other parental "Barke", we can call it "b". Finally, the rest of the alleles can be called heterozygous, "h"
mydata6[-c(1:2),]<-lapply(mydata6,function(x)
ifelse(x[-c(1,2)]==x[1],'b',
ifelse(x[-c(1,2)]==x[2],'m','h')))
write.csv(mydata6, file = "ChrX.csv", na="")
