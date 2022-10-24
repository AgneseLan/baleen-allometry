
#===========================================================#
#                                                           #
#       CURVES AND POINTS ANALYSES - MYSTICETI ONLY         #
#                                                           #
#===========================================================#

#CH.1 - Reading in curves and landmarks, resampling, and sliding
#Code adapted from Ryan Felice and Ellen Coombs

#LOAD LIBRARIES ----
#always do this first!!
library(tidyverse)
library(Morpho)
library(geomorph)
library(Rvcg)
library(paleomorph)
library(EMMLi)
library(qgraph)
library(ape)
library(geiger)
library(abind)
library("devtools")
library(SURGE)
library(magick)


#devtools::install_github("rnfelice/SURGE")

#DATA IMPORT ----

#Import table defining curves
curve_table <- read_csv('Data/curves_all.csv')

#curve_table$bone[11:12] <- "premaxilla" to change bone group if needed for plots

#Identify the folder where your pts files are (I have a different one for each suborder)
ptsfolder <- "Data/pts"

#Import the pts file names
ptslist <- dir(ptsfolder, pattern='.pts', recursive=F)

#Create curve info and indicate fixed points
my_curves <- create_curve_info(curve_table, n_fixed = 64)

##Check that your pts and ply match
#Identify ply folder
plyfolder <- "Data/ply"

#Make list of ply file names
plylist <-  dir(plyfolder, pattern='.ply', recursive=F)

#Make sublist to check names
ptslist2<-gsub(pattern="\\.pts$","",ptslist)
plylist2<-gsub(pattern="\\.ply$","",plylist)

#Check that names match up
identical(plylist2,ptslist2) #should be TRUE if order the same
#all.equal() to get explanation why they differ

#Save pts list to file - specimens names associated with classifiers
write.table(plylist2, file = "Output/specimens_list.csv", sep = ",", col.names = NA)

#Set WD where the landmark data is
setwd(ptsfolder)

#RESAMPLE AND FIX MISSING DATA ----

#Import pts data

#Capture the output first if it's too big to read in the console (and if you want to check each specimen is ok)
#Only do once, only checks if curve present in specimen
#Check Output folder path
capture.output(import_chkpt_data(ptslist, my_curves, subsampl = TRUE, verbose=TRUE), 
               file = "D:/OneDrive/my research project/R analyses/evo-devo-whales-2/Output/check_curves.txt", append = F)

#Create object with imported data
subsampled.lm <- import_chkpt_data(ptslist, my_curves, subsampl = TRUE, verbose=F)

#If you have any missing points, Checkpoint will set them to x,y,z=9999
#This makes for bad plots in checkLM, so switch them to NA
subsampled.lm[subsampled.lm == 9999] <- NA

###SET WD to ply from console!! -->

#Make sure you have converted your plys from binary to ASCII
#Check to make sure your curves look okay on each specimen
specs_tofix <- checkLM(subsampled.lm, path="", pt.size = 8, suffix=".ply", render = "s", begin = 35, point = "p")
#Write down specs numbers as you will need to clean environment to import them again

#Use exported ptslist to find names of specimens based on numbers

##Re-import starting from the beginning

#Check again
specs_tofix2 <- checkLM(subsampled.lm, path="", pt.size = 8, suffix=".ply", render = "s", begin = 15, point = "p")

#Create object with new resampled points
newpts <- subsampled.lm

###Fix missing landmarks####

#Create missing list 
misslist <- createMissingList(dim(newpts)[3])

#Get landmark numbers
for (j in 1:dim(newpts)[[3]]){
  misslist[[j]]<-which(is.na(newpts[,1,j]))
} 

#Create missing list table and file
df <- rbind(misslist)
df <- t(df)

misstable <- as.data.frame(df)

#Add names specimens
misstable <- cbind(misstable, specimens = dimnames(newpts)[[3]])

View(misstable)

#Write txt file (list class, no csv)
sink("D:/OneDrive/my research project/R analyses/evo-devo-whales-2/Output/missing_landmarks.txt")
print(misstable)
sink() 
#Convert in csv from excel if needed

##Fix missing landmarks
newpts2 <- fixLMtps(newpts)

#checkLM(newpts2$out, path="", pt.size = 1, suffix=".ply", render="s", begin = 36)

#SLIDE ----

#================================#
#      RUN Slider3D CODE         #
#================================#

###!!SET WD to ply from console!! -->

#Slide points on surface
slided4.all <- slider3d_2(newpts2$out, SMvector= my_curves$Sliding.LMs,
                          outlines = my_curves$Curves, 
                          #copy ply folder path to ensure it works - set as working directory from console to print
                          sur.path = "D:/OneDrive/my research project/R analyses/evo-devo-whales-2/Data/ply",
                          sur.name = NULL, 
                          meshlist = paste("",dimnames(newpts2$out)[[3]],".ply",sep=""), ignore = NULL,
                          sur.type = "ply", tol = 1e-10, deselect = FALSE, inc.check = FALSE,
                          recursive = TRUE, iterations = 3, initproc = TRUE,
                          pairedLM = 0, mc.cores = 1, bending=TRUE,
                          fixRepro = F, stepsize=0.2,
                          missingList = misslist)

#Name specimens slided data
dimnames(slided4.all[["dataslide"]])[3] <- dimnames(newpts2$out)[3]

#Re-estimate missing post sliding
slided4.all$dataslide[which(is.na(newpts))] <- NA

#Fix missing landmarks
slid.lms <- fixLMtps(slided4.all$dataslide)

#Extract landmark array from object
slidedlms <- slid.lms$out

#Check to see how they look (plys must be ASCII format)
specs_tofix_slid <- checkLM(slidedlms, path="", pt.size = 6, 
                            suffix=".ply", render = "s", begin = 171, point = "p")

#Save slided LMs as R data file
save(slidedlms, file = '"D:/OneDrive/my research project/R analyses/evo-devo-whales-2/Output/slidedlms.R')
#Save to Output folder


###### 
#Next - ch. 2 - absent bones 

