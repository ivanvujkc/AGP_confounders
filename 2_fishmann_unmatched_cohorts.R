# This script tests for differences in distribution of confounding factors among cohorts such as those that are outputs from the '1_unmatched_cohorts_nperms.R' script. For each condition, median P values will be calculated among all permutations of re-sampled cohorts for that given condition. Benjamini-Hochberg false discovery rate Q values will then be calculated for all median P values.

library(ggplot2); library(ape);library(vegan);library(reshape2);library(RColorBrewer); library(exactRankTests)

outpath<-"~/Downloads/outputdirectory/"
fishmannpath<-"~/Downloads/outputdirectory/fishmann_unmat/"
metadat<-read.csv(file="~/Downloads/metadat.csv",header=T,row.names=1)


all<-list.files(outpath,full.names=F,pattern="unmatched.csv")
allpaths<-list.files(outpath,full.names=T,pattern="unmatched.csv")

# microbiota-associated variables list:
varscat<-c("sex","alcohol_frequency","vegetable_frequency", "milk_cheese_frequency", "meat_eggs_frequency","bowel_movement_quality","whole_grain_frequency","sugary_sweets_frequency","salted_snacks_frequency")

varscontinuous<-c("age_years","bmi")

# The following set of calculations are for categorical variables (e.g. alcohol consumption frequency), and thus Fisher's exact tests are performed


for(k in 1:length(varscat))
{
resultmat<-data.frame(matrix(nrow=length(all),ncol=3))
colnames(resultmat)<-c("var","pvalue","distribution")
for(j in 1:length(allpaths)){
	ftab<-read.table(file= allpaths[j],header=T,sep=",",row.names=1)
	target1<-rownames(ftab)[ftab $target=="case"]
	mapp3<-merge(ftab, metadat,by=0)
	# mapp2<-metadat[rownames(metadat)%in%rownames(ftab),]
	# mapp3<-cbind(mapp2,"target")
	# mapp3[target1,"target"]<-1
	# mapp3[is.na(mapp3[,"target"]),"target"]<-0
	m2<-table(mapp3 $target, mapp3[,colnames(mapp3)%in% varscat[k]])
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
	mapp3<-merge(ftab, metadat,by=0)
	m2<-table(mapp3 $target, mapp3[,colnames(mapp3)%in% varscat[k]])
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
write.csv(resmed2,file=paste(fishmannpath,"/unmatch_fishers-",varscat[k],".csv",sep=""))
}









# The following set of calculations are for continuous variables (e.g. bmi), and thus Mann-Whitney U tests are performed
for(k in 1:length(varscontinuous))
{
resultmat<-data.frame(matrix(nrow=length(all),ncol=3))
colnames(resultmat)<-c("var","pvalue","distribution")

for(j in 1:length(allpaths)){
	ftab<-read.table(file= allpaths[j],header=T,sep=",",row.names=1)
	mapp3<-merge(ftab, metadat,by=0)
	m2<-mapp3[,c("target", varscontinuous[k])]
	m3<-t(m2)

	tryCatch(
	resultmat[j,2]<-wilcox.test(m2[,varscontinuous[k]]~ m2$target)$p.value
	,error=function(e){})
	
	tryCatch(
	resultmat[j,3]<-log(mean(m2[m2$target=="case",2])/mean(m2[m2$target=="control",2]),base=10)
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
write.csv(resmed2,file=paste(fishmannpath,"/unmatch_mann-",varscontinuous[k],".csv",sep=""))
}





library(ggplot2);library(viridis)
library(ggplot2);library(viridis)

fishmannpaths<-grep("unmatch",list.files(path= fishmannpath,full.names=T),value=T)
fishmann_names<-grep("unmatch",list.files(path= fishmannpath,full.names=F),value=T)
fishmannsummary<-read.table(file=fishmannpaths[1],sep=",",header=T)
fishmannsummary<-cbind(fishmannsummary,"varthatdiffers"=fishmann_names[1])
for(i in 2:length(fishmannpaths))
{
temp<-read.table(file=fishmannpaths[i],sep=",",header=T)
temp<-cbind(temp,"varthatdiffers"=fishmann_names[i])
fishmannsummary<-rbind(fishmannsummary,temp)
}

fishmannsummary[,ncol(fishmannsummary)]<-gsub(".csv","",as.vector(fishmannsummary[,ncol(fishmannsummary)]))

fishmannsummary$matchstat[fishmannsummary$BH_qval<0.05]<-"mismatched"
fishmannsummary$matchstat[fishmannsummary$BH_qval>0.05]<-"matched"
fishmannsummary2<-cbind(fishmannsummary,colsplit(fishmannsummary$varthatdiffers,"-",c("X3","variablename")))

ggplot(fishmannsummary2, aes(x= dis, y= variablename)) + geom_point(aes(colour= matchstat),size=3)+coord_flip()+ scale_color_viridis(option="viridis",discrete=T)+theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),panel.grid.major.y = element_blank(), panel.background = element_blank(),axis.text.x = element_text(angle = 90, hjust = 0,vjust=0.3))+scale_y_discrete(position = "right")





write.csv(fishmannsummary2,file=paste(fishmannpath,"/fishmann_summary.csv",sep=""))






