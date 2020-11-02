# This script creates n permutations of location-matched ('unmatched' for confounding variables) case-control cohorts for a list of 'conditions'. Throughout, 'conditions' refers to the host variable that defines cases and controls, for which tests are being performed to understand whether there is a microbiota differences between cases and controls. For example, our study utilized these scripts to test for differences between cases and controls for several diseases, in this case the diseases are compatible with our usage of the word 'condition'.
# It requires as input a csv file for each condition, with cases being marked as '1' under a column titled 'target'. Such files are the output of python scripts at github.com/jacksklar/AGPMicrobiomeHostPredictions, found in the 'Feature_Cohorts/binary_cohorts_no_matching' output folder. Such files must include all cases that are present in the mapping file. Downstream analysis will be facilitated if files are named such that the condition is one word (without spaces or underscores). Furthermore, to remove subjects with no information regarding a given subject in metadata ('NA' response), filenames denoting cases of a condition should be equivalent to column names referring to that condition in the metatdata file.

library(ggplot2);library(labdsv);library(reshape2);library(vegan)

# Define output filepath
outpath<-"~/Downloads/outputdirectory/"

map1<-read.csv(file="~/Downloads/metadat.csv",header=T,row.names=1)

# Impose inclusion/exclusion criteria, following are, as example, those for diabetes. These differ from special disease cohorts as defined in methods (IBD, ASD) and from remaining diseases. Namely, BMI is limited to 40 for other diseases and diabetics are excluded.
map2<-subset(map1,age_years>=20&age_years<=80&antibiotic_history<2&ibd==0&bmi>=12.5)
map2<-map2[rowSums(map2[,c("country_Canada","country_USA","country_United.Kingdom")])>0,]

# Define number of cohort permutations to create
nperms<-25

# Define maximum cohort size (total number of cases and controls)
maxcohortcutoff<-1000

tab<-read.csv("~/Downloads/merged2_s50_otus_mc01p.csv",header=T,row.names=1)

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

map3<-map3[!is.na(rowSums(map3)),]


# select cases from samples passing filters
allmet2<-map1[rownames(map3),]
target1<-rownames(subset(allmet2, diabetes==1&diabetes_type=="Type II diabetes"))
# write cases to file for later script
write.csv(target1,file="~/Downloads/cases_diabetes.csv",row.names=F)

# Remove all subjects not reporting information for given condition of interest
allnonnas<-rownames(allmet2[!is.na(allmet2[, "diabetes"]==0),])

# Create mapping file with just location metadata, randomize order of all subjects
tempmap<-map3[sample(allnonnas,size=length(allnonnas)),c(locvars)]


# Create Euclidean distance matrix based on location metadata
eudist<-dist(tempmap[,],method="euclidean")
eu2<-as.matrix(eudist)

for (l in 1:nperms)
{
  	# Select only case subjects with matching variable metadata in mapping file
	target1b<-target1[target1%in%rownames(map3)]
  	# Randomize order of cases, and impose maximum cohort limit (limiting number of cases, hence maxcohortcutoff/2)
  	target1<-sample(target1b,size=min(length(target1b),maxcohortcutoff/2))

  	# Create matrix containing names of all cases, dimensions to accommodate one control per case, and annotate by case/control status
  	ftab<-matrix(nrow=length(target1)*2,ncol=2)
    colnames(ftab)<-c("sampID","target")
    ftab[1:length(target1),1]<-target1
    ftab[1:length(target1),2]<-"case"
    ftab[(length(target1)+1):(length(target1)*2),2]<-"control"
    for (i in 1:length(target1))
    {
      sampname<-target1[i]
      # Rank controls by shortest Euclidean distance to case subject and extract sample names
      contsinord<-colnames(eu2[,order(eu2[sampname,])])
      # Remove from consideration all cases and previously selected controls
      contsinord2<-contsinord[!contsinord%in%ftab[,1]]
      # Select top 5 most similar controls for this given case
      topconts <-contsinord2[1:5]
      # Randomly select one from top 5 most similar controls
      randcont<-sample(topconts,size=1)
      ftab[i+length(target1),1]<-randcont
    }

write.csv(ftab,paste(outpath, "/diabetes_",l, "_unmatched.csv",sep=""),row.names=F)
}
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

map3<-map3[!is.na(rowSums(map3)),]

allmet2<-map1[rownames(map3),]

# Remove all subjects not reporting information for given condition of interest
allnonnas<-rownames(allmet2[!is.na(allmet2[, "diabetes"]==0),])

# Create mapping file with just location metadata, randomize order of all subjects
tempmap<-map3[sample(allnonnas,size=length(allnonnas)),c(locvars)]

# Create Euclidean distance matrix based on location metadata
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
write.csv(ftab,paste(outpath, "diabetes_",l, "_unmatched.csv",sep=""))
}
}
