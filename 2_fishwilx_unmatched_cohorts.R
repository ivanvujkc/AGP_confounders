# This script tests for differences in distribution of confounding factors among cohorts such as those that are outputs from the '1_unmatched_cohorts_nperms.R' script. For each condition, median P values will be calculated among all permutations of re-sampled cohorts for that given condition. Benjamini-Hochberg false discovery rate Q values will then be calculated for all median P values.

library(ggplot2); library(ape);library(vegan);library(reshape2);library(RColorBrewer); library(exactRankTests)

all<-list.files("~/pathto_unmatched_cohort_nperms_directory",full.names=F,pattern=".csv")
allpaths<-list.files("~/pathto_unmatched_cohort_nperms_directory",full.names=T,pattern=".csv")

metadat<-read.csv(file="~/pathto_metadata.csv",header=T,row.names=1)


# The following set of calculations are for categorical variables (e.g. alcohol consumption frequency), and thus Fisher's exact tests are performed
resultmat<-data.frame(matrix(nrow=length(all),ncol=3))
colnames(resultmat)<-c("var","pvalue","distribution")

for(j in 1:length(allpaths)){
	ftab<-read.table(file= allpaths[j],header=T,sep=",",row.names=1)
	target1<-rownames(ftab[ftab $target=="case",])
	mapp2<-metadat[rownames(metadat)%in%rownames(ftab),]
	mapp3<-cbind(mapp2,"target")
	mapp3[target1,"target"]<-1
	mapp3[is.na(mapp3[,"target"]),"target"]<-0
	m2<-table(mapp3 $target, mapp3 $alcohol_frequency)
	m2[is.na(m2)]<-0
	m3<-t(m2)

	tryCatch(
	rm(f1)
	,error=function(e){})
	
	tryCatch(
	f1<-fisher.test(m2,workspace = 200000000)
	,error=function(e){})
	
	tryCatch(
	resultmat[j,2]<-f1$p.value
	,error=function(e){})
	
	tryCatch(
	resultmat[j,3]<-log(sum(m3[4:5,2])/sum(m3[1:2,2]),base=10)-log(sum(m3[4:5,1])/sum(m3[1:2,1]),base=10)
	,error=function(e){})
		
resultmat[j,1]<-all[j]
}

# The following block is performed in the case that the Fisher's exact test workspace was not large enough to accommodate all metadata in the prior block. Thus, P values are simulated.
for(j in which(is.na(resultmat$pvalue))){
	ftab<-read.table(file= allpaths[j],header=T,sep=",",row.names=1)
	target1<-rownames(ftab[ftab $target=="case",])
	mapp2<-metadat[rownames(metadat)%in%rownames(ftab),]
	mapp3<-cbind(mapp2,"target")
	mapp3[target1,"target"]<-1
	mapp3[is.na(mapp3[,"target"]),"target"]<-0
	m2<-table(mapp3 $target, mapp3 $alcohol_frequency)
	m2[is.na(m2)]<-0
	m3<-t(m2)

	tryCatch(
	rm(f1)
	,error=function(e){})
	
	tryCatch(
	f1<-fisher.test(m2,simulate.p.value=T)
	,error=function(e){})

	tryCatch(
	resultmat[j,2]<-f1$p.value
	,error=function(e){})
	
	tryCatch(
	resultmat[j,3]<-log(sum(m3[4:5,2])/sum(m3[1:2,2]),base=10)-log(sum(m3[4:5,1])/sum(m3[1:2,1]),base=10)
	,error=function(e){})	
	
resultmat[j,1]<-all[j]
}

resultmat<-resultmat[order(resultmat$pvalue),]
resultmat3<-cbind(resultmat,colsplit(resultmat $var,"_",c("dis","nperm","cond","mat")))
ressum<-tapply(resultmat3$pvalue,resultmat3$dis,summary)
ressum2<-cbind(unlist(ressum),names(unlist(ressum)))
ressum3<-cbind(ressum2,colsplit(ressum2[,2],"\\.",c("dis","stat")))
resmed<-subset(ressum3,stat=="Median")
resmed2<-cbind(resmed,BH_qval=p.adjust(as.numeric(as.vector(resmed [,1])),method="BH"))
colnames(resmed2)[1]<-"pvalue"
write.csv(resmed2,file="~/pathto_outputdirectory_mismatchtest/unmatch_fishers-alcfreq.csv")










# The following set of calculations are for continuous variables (e.g. bmi), and thus Mann-Whitney U tests are performed

resultmat<-data.frame(matrix(nrow=length(all),ncol=3))
colnames(resultmat)<-c("var","pvalue","distribution")

for(j in 1:length(allpaths)){
	ftab<-read.table(file= allpaths[j],header=T,sep=",",row.names=1)
	target1<-rownames(ftab[ftab $target=="case",])
	mapp2<-metadat[rownames(metadat)%in%rownames(ftab),]
	mapp3<-cbind(mapp2,"target")
	mapp3[target1,"target"]<-1
	mapp3[is.na(mapp3[,"target"]),"target"]<-0
	m2<-mapp3[,c("target","bmi")]
	m3<-t(m2)

	tryCatch(
	resultmat[j,2]<-wilcox.test(m2$bmi~ m2$target)$p.value
	,error=function(e){})
	
	tryCatch(
	resultmat[j,3]<-log(mean(m2[m2$target=="case",2])/mean(m2[m2$target==0,2]),base=10)
	,error=function(e){})

resultmat[j,1]<-all[j]
}

resultmat<-resultmat[order(resultmat$pvalue),]
resultmat3<-cbind(resultmat,colsplit(resultmat $var,"_",c("dis","nperm","cond","mat")))
ressum<-tapply(resultmat3$pvalue,resultmat3$dis,summary)
ressum2<-cbind(unlist(ressum),names(unlist(ressum)))
ressum3<-cbind(ressum2,colsplit(ressum2[,2],"\\.",c("dis","stat")))
resmed<-subset(ressum3,stat=="Median")
resmed2<-cbind(resmed,BH_qval=p.adjust(as.numeric(as.vector(resmed [,1])),method="BH"))
colnames(resmed2)[1]<-"pvalue"
write.csv(resmed2,file="~/pathto_outputdirectory_mismatchtest/unmatch_mannwhitney-bmi.csv")
