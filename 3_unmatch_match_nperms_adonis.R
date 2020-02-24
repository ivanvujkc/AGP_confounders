# The following script performs Adonis (PERMANOVA) tests and thus calculates P values and F statistics for 25 unmatched and 25 confounder variable-matched cohorts for each condition. This script requires assessment of which confounder variables differ between cases and controls (output of "2_fishwilx_unmatched_cohorts.R" script), as well as lists of cases within metadata (input of "1_unmatch_nperms"). N.b. inclusion/exclusion criteria differed in our study for: autism spectrum disease, type 2 diabetes, and inflammatory bowel disease, as written in the manuscript. Below, inclusion/exclusion criteria for all remaining diseases is shown.

library(ggplot2);library(labdsv);library(reshape2);library(vegan)

fishwilx2 <-read.csv(file="~/pathto_outputdirectory_mismatchtest/fishwilx_summary.csv",header=T)
dis<-unique(fishwilx2$dis)

cohortnames <-list.files("~/pathto_casesdirectory",pattern=".csv",full.names=F)
cohortpaths <-list.files("~/pathto_casesdirectory",pattern=".csv",full.names=T)

# Ensure that vectors "cohortnames" and "dis" include the same list of conditions
cohortnames%in%dis
dis%in%cohortnames

# Define number of cohort permutations to create
nperms<-25

# Define maximum cohort size (total number of cases and controls)
maxcohortcutoff<-1000

# Define number of Adonis (PERMANOVA) permutations to perform for P value calculation (can significantly increase CPU time)
adonperms<-999

# Define output filepath for case-control cohorts
outpath<-"~/pathto_outputdirectory_matunmatcohorts/"

# Define output filepath for visualizations and Adonis results
outpathadon<-"~/pathto_outputdirectory_matunmat_adonis_graphs/"

tab<-read.table("~/pathto_otutable.txt",sep="\t",header=T,row.names=1)


adonmatall<-data.frame(matrix(ncol=10))
colnames(adonmatall)<-c("Df","SumsOfSqs","MeanSqs","Fstat","R2","P_value","sampIDparam","numsamps","distparam","disease")


for(h in 1:length(cohortnames)){

unmatched<-read.csv(file=cohortpaths[h],header=T,row.names=1)

targetu<-rownames(unmatched[unmatched$target==1,])

map1<-read.csv(file="~/pathto_metadata.csv",header=T,row.names=1)

# Impose inclusion/exclusion criteria
map2<-subset(map1,age_years>=20&age_years<=80&antibiotic_history<2&ibd==0&diabetes==0&bmi>=12.5&bmi<=40)
map2<-map2[rowSums(map2[,c("country_Canada","country_USA","country_United.Kingdom")])>0,]

# Remove subjects not in OTU/ASV table
map2<-map2[rownames(map2)%in%colnames(tab),]

# Scale and center all matching variables
map2$latitude<-scale(map2$latitude)
map2$longitude<-scale(map2$longitude)
map2$sex<-scale(map2$sex)
map2$age_years<-scale(map2$age_years)
map2$bmi<-scale(map2$bmi)
map2$alcohol_frequency<-scale(map2$alcohol_frequency)
map2$vegetable_frequency<-scale(map2$vegetable_frequency)
map2$milk_cheese_frequency<-scale(map2$milk_cheese_frequency)
map2$meat_eggs_frequency<-scale(map2$meat_eggs_frequency)
map2$bowel_movement_quality <-scale(map2$bowel_movement_quality)
map2$whole_grain_frequency <-scale(map2$whole_grain_frequency)
map2$sugary_sweets_frequency <-scale(map2$sugary_sweets_frequency)
map2$salted_snacks_frequency <-scale(map2$salted_snacks_frequency)
map3<-map2[,c("latitude","longitude","sex","age_years","bmi", "alcohol_frequency","vegetable_frequency", "milk_cheese_frequency", "meat_eggs_frequency","bowel_movement_quality","whole_grain_frequency","sugary_sweets_frequency","salted_snacks_frequency")]
locvars<-c("longitude","latitude")
sub1<-subset(fishwilx2,dis== cohortnames[h])
varsthatdiffer<-sub1[sub1$BH_qval=="mismatched","variablename"]
map3<-map2[,c(locvars,varsthatdiffer)]
map3<-map3[!is.na(rowSums(map3)),]

vars<-c("matched","unmatched")

justvars<-colnames(map3)[!colnames(map3)%in%c("longitude","latitude")]

allmet2<-allmetadat[rownames(map3),]

# Remove all subjects not reporting information for given condition of interest
allnonnas<-rownames(allmet2[!is.na(allmet2[, dis2[h]]==0),])

# Create mapping file with just location metadata, randomize order of all subjects
tempmap<-map3[sample(allnonnas,size=length(allnonnas)),c(locvars)]

# Create Euclidean distance matrix based on location metadata
eudist<-dist(tempmap[,],method="euclidean")
eu2<-as.matrix(eudist)


# target1<-targetu[targetu%in%rownames(map3)]

adonmat<-data.frame(matrix(nrow=length(vars)*nperms,ncol=10))
colnames(adonmat)<-colnames(adonmatall)
for(j in 1:length(vars))
{
  if(vars[j]=="unmatched")
	  {
	    tempmap<-map3[allnonnas,locvars]	
	  }	else if(vars[j]=="matched")	{
	    tempmap<-map3[allnonnas,c(locvars, justvars)]
	  }
	  
# Create Euclidean distance matrix based on given metadata
eudist<-dist(tempmap[,],method="euclidean")
eu2<-as.matrix(eudist)

for (l in 1:nperms)
   {
    # Select only case subjects with matching variable metadata in mapping file
	target1b<-targetu[targetu%in%rownames(map3)]
  	# Randomize order of cases, and impose maximum cohort limit (limiting number of cases, hence maxcohortcutoff/2)
  	target1<-sample(target1b,size=min(length(target1b),maxcohortcutoff/2))
  	
  	# Create matrix containing names of all cases, dimensions to accommodate one control per case, and annotated by case/control status
  	ftab<-matrix(nrow=length(target1)*2,ncol=2)
    colnames(ftab)<-c("sampID","target")
    ftab[1:length(target1),1]<-target1
    ftab[1:length(target1),2]<-"case"
    ftab[(length(target1)+1):(length(target1)*2),2]<-"control"
    for (i in 1:length(target1))
    {
      sampname<-target1[i]
      # Rank controls by shortest Euclidean distance to case subject
      contsinord<-colnames(eu2[,order(eu2[sampname,])])
      # Remove from consideration all cases and previously selected controls
      contsinord2<-contsinord[!contsinord%in%ftab[,1]]
      # Select top 5 most similar controls for this given case
      topconts <-contsinord2[1:5]
      # Randomly select one from top 5 most similar controls
      randcont<-sample(topconts,size=1)
      ftab[i+length(target1),1]<-randcont
    }

rownames(ftab)<-ftab[,"sampID"]
write.csv(ftab,file=paste(outpath,as.character(cohortnames[h]),"-", vars[j],"-",l,".csv",sep=""))
tab2<-tab[,colnames(tab)%in%rownames(ftab)]

# Take out samples from mapping file that are not in OTU table
ftab2<-ftab[rownames(ftab)%in%colnames(tab2),]
tab3<-tab2[,colnames(tab2)%in%rownames(ftab2)]

tOTUtab<-t(tab3)

canb<-as.matrix(dist(tOTUtab, method = "canberra", diag = T, upper = T))
canb3<-canb[order(rownames(canb)),order(colnames(canb))]
ftab2 <-as.data.frame(ftab2[order(rownames(ftab2)),])
formu<-formula(paste("canb3 ~ target"))
MDSadonis <-adonis(formu, data= ftab2, permutations=adonperms)

    k<-(j-1)*nperms+l
    adonmat[k,1]<-MDSadonis$aov.tab[1,1]
    adonmat[k,2]<-MDSadonis$aov.tab[1,2]
    adonmat[k,3]<-MDSadonis$aov.tab[1,3]
    adonmat[k,4]<-MDSadonis$aov.tab[1,4]
    adonmat[k,5]<-MDSadonis$aov.tab[1,5]
    adonmat[k,6]<-MDSadonis$aov.tab[1,6]
    adonmat[k,7]<-vars[j]
    adonmat[k,8]<-nrow(canb3)
    adonmat[k,9]<-toString(colnames(tempmap))
	adonmat[k,10]<-as.character(cohortnames[h])
	  } 
}

adonmatall<-rbind(adonmatall,adonmat)


write.table(adonmat,file=paste(outpathadon, cohortnames,"-adonis_results.csv",sep=""),sep=",",col.names=NA)

adonsum<-aggregate(adonmat[,"F.Model"], list(adonmat$sampIDparam), mean)
adonmat $sampIDparam = factor(adonmat $sampIDparam, adonsum[order(adonsum[,2]),1])
ggplot(adonmat,aes(x= sampIDparam,y=F.Model)) + theme_bw() +geom_boxplot(fatten = NULL,outlier.size=-1) + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = 0.75, size = 1, linetype = "solid") +geom_jitter(height=0,width=0.2,alpha=0.4) +theme(axis.text.x = element_text(angle = 90, hjust =1,vjust=0.3))
ggsave(filename=paste(outpathadon, cohortnames,"-adonis.pdf",sep=""),height= 3.5,width=4.7)
}

 write.csv(adonmatall,file=paste(outpathadon,"/alldis.csv",sep=""))
