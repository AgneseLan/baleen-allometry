
#===========================================================#
#                                                           #
#       CURVES AND POINTS ANALYSES - MYSTICETI ONLY         #
#                                                           #
#===========================================================#


#CH.3 - Prepare final dataset for analysis, run GPA and PCA

#LOAD LIBRARIES ----
#always do this first!!
library(geomorph) 
library(geiger)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggrepel)
library(gginnards)
library(ggphylomorpho)
library(ggfortify)
library(RColorBrewer) 
library(borealis)
library(ggthemes)
library(ggpubr)
library(ggplotify)
library(Morpho)
library(rphylopic)
library(png)
library(gridExtra)
library(phytools)
library(evomap)
library(car)
library(Rvcg)
library(scales)

#devtools::install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")

#DATA PREP ----

###SET WD to root folder from console!! -->

#Import classifiers
classifiers <- read_csv("Data/specimens_myst.csv")
glimpse(classifiers)

#Make sure the specimens are in the same order before proceeding!!!
identical(dimnames(final_dataset)[[3]], classifiers$specimen, attrib.as.set = T)

##Order shape data by category, useful for plot legend
#Check levels/category names
as.factor(classifiers$category)

#Order shape data as category
order_dataframe <- geomorph.data.frame(raw_data = final_dataset, category = classifiers$category)
#Check specimens order
dimnames(order_dataframe$raw_data)[[3]]

#Order dataframe
order_dataframe$category <- factor(order_dataframe$category,
                                   levels = c("1-early", "2-late/new", "3-immature", "4-adult"))

#Create new shape data object ordered by category
shape_array <- order_dataframe$raw_data[,,order(order_dataframe$category)]
#Check it worked
dimnames(shape_array)[[3]]

##Order classifiers data by category, useful for plot legend
#Make factor for variable
classifiers$category <- factor(classifiers$category, 
                               levels = c("1-early", "2-late/new", "3-immature", "4-adult")) #copy from string printed with the code above
#Order
classifiers <- classifiers[order(classifiers$category),]
View(classifiers)

#Check specimens and classifiers are in the same order
identical(dimnames(shape_array)[[3]], classifiers$specimen, attrib.as.set = T) 

##Save mesh with plotted landmarks
#Find mean specimen raw data
findMeanSpec(shape_array)
#if fetus with missing bones look for well scanned adult or juvenile
#Get spec number
match("Bal.acut_CAS23867", dimnames(shape_array)[[3]])

#Import simplified ply
#Less faces, no holes or isolated triangles
refmesh_all <- vcgImport("Data/refmesh_myst.ply")

#Plot on surface
shade3d(refmesh_all, col = "white", alpha = 0.5)
spheres3d(shape_array[fixed_LMs,,52], col =  "firebrick", type = "s",
          radius = 8, aspect = T, main = "mean",axes = F, main = F, fov = 0)
spheres3d(shape_array[-fixed_LMs,,52], col =  "tomato", type = "s",
          radius = 6, aspect = T, main = "mean",axes = F, main = F, fov = 0)
text3d(shape_array_LM[,1,52], shape_array_LM[,2,52], shape_array_LM[,3,52], 
       texts = fixed_LMs, pos = 4, offset = 1, font = 2) #change pos

rgl.snapshot(filename = "Output/landmarks_dorsal.png") 
rgl.snapshot(filename = "Output/landmarks_lateral1.png") 
rgl.snapshot(filename = "Output/landmarks_lateral2.png") 
rgl.snapshot(filename = "Output/landmarks_ventral.png") 
rgl.snapshot(filename = "Output/landmarks_posterior.png") 
rgl.snapshot(filename = "Output/landmarks_anterior.png") 

play3d(spin3d(axis = c(1, 0,0), rpm = 10), duration = 6)
movie3d(spin3d(axis = c(1, 0,0), rpm = 10), duration = 6, movie = "landmarks" ,dir = "Output/")

#Check for outliers in raw data shape array, they would be displayed in red
#Might be due to absent bones, check misstable list
plotOutliers(shape_array)
#Plot outliers landmarks to search for possible problems - check file name to find number in classifiers
checkLM(shape_array, path="Data/ply/", pt.size = 5, suffix=".ply", render = "s", begin = 162, point = "s")

#Save txt file of output
sink("Output/outliers_raw.txt")
print(plotOutliers(shape_array))
sink() 
#Mark where outliers end based on plot

#Save specimens names as object
specimens <- dimnames(shape_array)[[3]]

#Save Id as object, useful for later analysis
Ids <- classifiers$code

#Save growth stages as factor, useful for later analysis
categories <- as.factor(classifiers$category)
#Check how many colors are needed
levels(categories) #4
#To capitalize or change spelling
levels(categories) <- c("Early Fetus", "Late Fetus/Neonate", "Juvenile", "Adult")

#Save genera as factor, useful for later analysis
genera <- as.factor(classifiers$genus2)
#Check how many colors are needed
length(levels(genera)) #9

#Save families as factor, useful for later analysis
families <- as.factor(classifiers$family)
#Check how many colors are needed
length(levels(families)) #6
#To capitalize or change spelling
levels(families) <- c("Balaenidae" ,     "Balaenopteridae" , "Delphinidae"  ,       
                      "Neobalaenidae" , "Phocoenidae",     "Physeteroidea")

#Save groups as factor, useful for later analysis
groups <- as.factor(classifiers$group)
#Check how many colors are needed
groups #2
#To capitalize or change spelling
levels(groups) <- c("Mysticeti","Odontoceti")

#Save feeding modes as factor, useful for later analysis
feeding <- as.factor(classifiers$Feeding_BL20)
#Check how many colors are needed
length(levels(feeding)) #6
#To capitalize or change spelling
levels(feeding) <- c("Biting/Crushing-Suction", "Filter (lunge)",
                     "Filter (skim)", "Filter (suction)", "Suction")

##Create project palettes----
mypalette_paired <- brewer.pal(12,"Paired")
image(1:12, 1, as.matrix(1:12), col = mypalette_paired, xlab = "Paired",
      ylab = "", yaxt = "n")

mypalette_hueyellowblue <- brewer.pal(9,"YlGnBu")
image(1:9, 1, as.matrix(1:9), col = mypalette_hueyellowblue, xlab = "Y-B Hues",
      ylab = "", yaxt = "n")

mypalette_huepurple <- brewer.pal(9,"RdPu")
image(1:9, 1, as.matrix(1:9), col = mypalette_huepurple, xlab = "Purple Hues",
      ylab = "", yaxt = "n")

mypalette_huegreen <- brewer.pal(9,"YlGn")
image(1:9, 1, as.matrix(1:9), col = mypalette_huegreen, xlab = "Green Hues",
      ylab = "", yaxt = "n")

mypalette_blue <- as.matrix(ggthemes_data[["tableau"]][["color-palettes"]][["ordered-sequential"]][["Blue"]][["value"]])
image(1:20, 1, as.matrix(1:20), col = mypalette_blue, xlab = "Blue",
      ylab = "", yaxt = "n")

#Palette for 27 genera - based on genus2, different shades for myst and odont
mypalette_myst <- colorRampPalette(c(mypalette_huepurple[3:6], mypalette_huepurple[8]))
mypalette_myst(6)
plot(rep(1,6),col=mypalette_myst(6),pch=19,cex=3)
mypalette_myst2 <- mypalette_myst(6)

mypalette_odont <- c(mypalette_huegreen[4],
                     mypalette_hueyellowblue[6],
                     mypalette_hueyellowblue[8]) #mix different hues to help
plot(rep(1,3),col=mypalette_odont,pch=19,cex=3)

#Make palette for taxa by listing genera
levels(genera)

mypalette_taxa <- c(mypalette_myst2[4], mypalette_myst2[5], mypalette_myst2[1], mypalette_myst2[2], #"Balaena"   "Balaenoptera_l" "Balaenoptera_s" "Caperea" 
                    mypalette_myst2[6], # "Eschrichtius" 
                    mypalette_odont[1], #   "Kogia" 
                    mypalette_odont[2], mypalette_myst2[3], # "Lagenorhynchus" "Megaptera"     
                    mypalette_odont[3]) #"Phocoena"      

plot(rep(1,9),col=mypalette_taxa ,pch=19,cex=3, main = "Taxa colors", ylab = "", xlab = "" ,cex.main = 2, 
     yaxt = "n")
title(xlab = "1-Balaena 2-Balaenoptera_l 3-Balaenoptera_s 4-Caperea 5-Eschrichtius 
   6-Kogia 7-Lagenorhynchus 8-Megaptera 9-Phocoena", cex.lab = 1.3, font.lab = 3, line = -2)
text(x = seq_along(1:9), y = 1.05, labels = seq_along(1:9))

#Palette for families
levels(families)

mypalette_families <- c(mypalette_myst2[4], mypalette_myst2[1], mypalette_odont[2], #"balaenidae"      "balaenopteridae" "delphinidae"     
                        mypalette_myst2[2], mypalette_odont[3], mypalette_odont[1]) #"neobalaenidae"   "phocoenidae"     "physeteroidea"

plot(rep(1,6),col=mypalette_families ,pch=19,cex=3, main = "Families colors", ylab = "", xlab = "" ,cex.main = 2)
title(xlab = "1-Balaenidae 2-Balaenopteridae 3-Delphinidae 4-Neobalaenidae 
5-Phocoenidae 6-Physeteroidea", cex.lab = 1.3, font.lab = 3, line = -3)
text(x = seq_along(1:6), y = 1.05, labels = seq_along(1:6))

#Palette for categories - early, late/new, immature, adult
mypalette_category <- c(mypalette_blue[3,], mypalette_blue[7,], mypalette_blue[13,], mypalette_blue[18,])
image(1:4, 1, as.matrix(1:4), col = mypalette_category, main = "Categories colors", 
      xlab =  "1-early 2-late/new 3-immature 4-adult", cex.lab = 1.3, cex.main =2,
      ylab = "", yaxt = "n")

#Palette for groups (Mysticeti, Odontoceti)
mypalette_groups <- c(mypalette_myst2[4], mypalette_odont[1])
image(1:2, 1, as.matrix(1:2), col = mypalette_groups, main = "Groups colors", xlab = "1-Mysticeti 2-Odontoceti", 
      cex.lab = 1.2, cex.main =2,ylab = "", yaxt = "n", xaxt = "n")
axis(1, at = c(1, 2))

#Palette for feeding modes (biting/crushing, biting/crushing/suction, filter, filter_skim, filter_suc, suction)
mypalette_feeding <- c(mypalette_odont[2], mypalette_myst2[3], 
                       mypalette_myst2[2], mypalette_myst2[6], mypalette_odont[1] )
plot(rep(1,5),col=mypalette_feeding ,pch=19,cex=3, main = "Feeding modes colors", ylab = "", xlab = "" ,cex.main = 2)
title(xlab = "1-biting/crushing/suction 
      2-filter 3-filter_skim 4-filter_suc 
      5-suction",
      cex.lab = 1.3, line = -4)
text(x = seq_along(1:5), y = 1.05, labels = seq_along(1:5))

#Create shape palette 2 groups (Mysticeti, Odontoceti), 4 categories and 13 families
shapes <- c(21,22)
shapes_cat <- c(23,24,22,21)
shapes_fam <- c(17, 19, 3, 15, 8, 9) #mysticeti solid, odontoceti empty
shapes_feed <- c(2, 17, 19, 15, 5)

##Images for plots
Bals <- readPNG("Data/b.bona.png")
Ball <- readPNG("Data/b.physalus.png")
Bala <- readPNG("Data/balaena.png")
Cap <- readPNG("Data/caperea.png")
Esch <- readPNG("Data/eschrichtius.png")
Kog <- readPNG("Data/kogia.png")
Lage <- readPNG("Data/lagenorhynchus.png")
Meg <- readPNG("Data/megaptera.png")
Phoc <- readPNG("Data/phocoena.png")

myst <- readPNG("Data/megaptera.png")
odont <- readPNG("Data/lagenorhynchus.png")

bala <- readPNG("Data/balaena.png")
balpt <- readPNG("Data/b.physalus.png")
delp <- readPNG("Data/lagenorhynchus.png")
neob <- readPNG("Data/caperea.png")
phoc <- readPNG("Data/phocoena.png")
phys <- readPNG("Data/kogia.png")

#GPA ALIGNMENT ----

#Procrustes alignment, should also show mean config coordinates
gpa <- gpagen(shape_array) 

#Save Centroid size as object
Csize <- gpa$Csize 
#Log-transform Centroid size as object
logCsize <- log10(Csize) 

#Save mean shape to create links
mean_shape <- gpa$consensus 

#Coordinates of all specimens after GPA alignment
coords <- gpa$coords 

#Plot all specimens with mean to check that all landmarks are ok
plotAllSpecimens(coords, mean = TRUE, label = F, plot.param = list(pt.cex = 0.05, mean.cex = 3, mean.bg = "black"))
#Save screenshot of 3D viewer
rgl.snapshot(filename = "Output/plot_gpa_points.png") 
rgl.snapshot(filename = "Output/plot_gpa_points1.png") 

#Check for outliers, they would be displayed in red - most immature ones are normal as outliers
plotOutliers(coords)
#Plot landmarks from outliers in 3D to check how they look
spheres3d(coords[,,28], r = 0.002)

#checkLM(shape_array, path="", pt.size = 2, suffix=".ply", render="s", begin = 65) 
#to run if needed to check plotting of points on mesh

##Make data frame for analyses in geomorph
gdf <- geomorph.data.frame(coords = coords,
                           Id = classifiers$code, genus = classifiers$genus2, category = classifiers$category,
                           family = classifiers$family, group = classifiers$group,
                           TL_100 = classifiers$TL_100, BZW_100 = classifiers$BZW_100, 
                           feeding = classifiers$Feeding_BL20, size = logCsize)
glimpse(gdf)

gdf$stage <- if_else(gdf$category %in% levels(classifiers$category)[1:2], "prenatal", "postnatal")

#============================================================#
#             CHECK GPA OUTLIERS separate script             #
#============================================================#

#PCA COMPLETE DATASET ----

#Run PCA on complete dataset
PCA_all <- gm.prcomp(gdf$coords)

#List of PC components and proportion of variation
PCA_all 

#Save PCA results to file
sink("Output/PCA_all_components.txt")
print("PCA complete dataset")
print(PCA_all)
sink() 

#Change row names to codes to make plot readable
row.names(PCA_all$x) <- gdf$Id

##View plot
pca_plot <- plot(PCA_all, main = "PCA all data - PC1-PC2",  pch = 21, #title and type of point to be used
                 col = "deeppink",   #border of points
                 bg = "deeppink",    #fill of points
                 cex = 1,            #size of points (1=regular)
                 font.main = 2)       #bold font for title
#Add quick labels to plot
text(x = PCA_all$x[,1], y = PCA_all$x[,2], labels = rownames(PCA_all$x), 
     pos = 1,       #position relative to data point
     offset = 0.5,  #distance from data point
     cex = 0.75)    #font size (1=regular)

#Save PC scores as object to use later
pcscores_all <- PCA_all$x 

#Save shapes of extremes for axes used in plot
PC1min_all <- PCA_all[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_all <- PCA_all[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_all <- PCA_all[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_all <- PCA_all[["shapes"]][["shapes.comp2"]][["max"]] 

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - using spheres3D for points
#PC1min colors
#spheres3d(mean_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
PC1min_all_points <- spheres3d(PC1min_all, radius=.001, color = col_modules)

rgl.snapshot(filename = "Output/PC1min_all.png") 
clear3d()

#PC1max colors
PC1max_all_points <- spheres3d(PC1max_all, radius=.001, color = col_modules)

rgl.snapshot(filename = "Output/PC1max_all.png") 
clear3d()

#PC2min colors
PC2min_all_points <- spheres3d(PC2min_all, radius=.001, color = col_modules)

rgl.snapshot(filename = "Output/PC2min_all.png") 
clear3d()

#PC2max colors
PC2max_all_points <- spheres3d(PC2max_all, radius=.001, color = col_modules)

rgl.snapshot(filename = "Output/PC2max_all.png") 
clear3d()

##Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_all <- as_tibble(pcscores_all)
#Add labels and other attributes to tibble as columns
pcscores_all <- pcscores_all %>% mutate(specimens = gdf$Id, group = gdf$group, category = gdf$category,
                                        genus = gdf$genus, family = gdf$family, size = gdf$size, 
                                        feeding = gdf$feeding)
glimpse(pcscores_all)

#Nice PCA plot with stages and groups
PCA_all_ggplot <- ggplot(pcscores_all, aes(x = Comp1, y = Comp2, label = specimens, colour = genus, fill = genus))+
  geom_point(size = 3, aes(shape = category))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 40)+
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = shapes_cat)+            #legend and color adjustments
  scale_colour_manual(name = "Genus", labels = levels(genera), #copy from as.factor(genera)
                      values = mypalette_taxa, aesthetics = c("colour","fill"), 
                      guide = guide_legend(label.theme = element_text(size =10, angle = 0, face = "italic")))+
  theme_bw()+
  xlab("PC 1 (45.41%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (23.85%)")+
  theme(legend.title = element_text(size = 12, face = "bold")) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_ggplot

#Make hulls for PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
#Divide by group first
hulls_all_category_myst <- pcscores_all %>%
  filter(group %in% c("mysticeti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

hulls_all_category_odont <- pcscores_all %>%
  filter(group %in% c("odontoceti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
PCA_all_category_ggplot <- ggplot(pcscores_all, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = family))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_all_category_myst, aes(x = x, y = y, fill = category), 
               alpha = .1, size = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_all_category_odont, aes(x = x, y = y, fill = category), 
               alpha = .1, size = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Family", labels = levels(families),
                     values = shapes_fam)+
  theme_bw()+
  xlab("PC 1 (45.41%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (23.85%)")+
  theme(legend.title = element_text(size = 11, face = "bold"))+
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_category_ggplot

#Add phylopics for groups
PCA_all_category_ggplot <- 
  PCA_all_category_ggplot +
  add_phylopic(myst, alpha = 1, x = -0.05, y = -0.15, ysize = 0.03, color = "gray30")+
  add_phylopic(odont, alpha = 1, x = 0.2, y = 0.2, ysize = 0.05, color = "gray20")
PCA_all_category_ggplot

#Make hulls for PCA plot with hulls around genera
hulls_all_genus <- pcscores_all %>%
  group_by(genus) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around genera 
PCA_all_genus_ggplot1 <- ggplot(pcscores_all, aes(x = Comp1, y = Comp2))+
  geom_polygon(data = hulls_all_genus, aes(x = x, y = y, group = genus, colour = genus, linetype = group), inherit.aes = F,
               size = 0.8, alpha = 0.005, show.legend = FALSE)+ #colored hulls with transparency
  geom_point(size = 4, aes(shape = category, alpha = category, colour = genus, fill = genus))+
  scale_colour_manual(name = "Genus", labels = levels(genera), 
                      values = mypalette_taxa, aesthetics = c("colour","fill"),
                      guide = guide_legend(label.theme = element_text(size =10, angle = 0, face = "italic")))+ #to be ordered as they appear in tibble
  scale_alpha_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = c(0.3, 0.5, 0.7, 1))+            #legend and color adjustments
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), 
                     values = shapes_cat)+
  scale_linetype_manual(values = c(1, 2))+
  theme_bw()+
  xlab("PC 1 (45.41%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (23.85%)")+
  theme(legend.title = element_text(size = 11, face = "bold"))+
  guides(shape = guide_legend(override.aes = list(fill = "grey10", colour = "grey10")))

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_genus_ggplot1

#Nice PCA plot with hulls around genera - myst only
pcscores_all_myst <- pcscores_all %>% filter(group == "mysticeti")
hulls_all_genus_myst <- hulls_all_genus %>% filter(group == "mysticeti")

PCA_all_genus_myst_ggplot <- ggplot(pcscores_all_myst, aes(x = Comp1, y = Comp2))+
  geom_polygon(data = hulls_all_genus_myst, aes(x = x, y = y, group = genus, colour = genus, linetype = family), inherit.aes = F,
               size = 1, alpha = 0.005, show.legend = FALSE)+ #colored hulls with transparency
  geom_point(size = 4, aes(shape = category, alpha = category, colour = genus, fill = genus))+
  scale_colour_manual(name = "Genus", labels = levels(as.factor(pcscores_all_myst$genus)), 
                      values = mypalette_taxa[-c(6,7,9)], aesthetics = c("colour","fill"),
                      guide = guide_legend(label.theme = element_text(size =10, angle = 0, face = "italic")))+ #to be ordered as they appear in tibble
  scale_alpha_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = c(0.3, 0.5, 0.7, 1))+            #legend and color adjustments
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), 
                     values = shapes_cat)+
  scale_linetype_manual(values = c(3,1,2))+
  theme_bw()+
  xlab("PC 1 (45.41%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (23.85%)")+
  theme(legend.title = element_text(size = 11, face = "bold"), 
        legend.position = c(0.85,0.23), legend.key = element_blank(),
        legend.background=element_blank())+
  guides(shape = guide_legend(override.aes = list(fill = "grey10", colour = "grey10")),
         color = "none", fill = "none")

#Visualize plot and save as PDF using menu in bar on the right
PCA_all_genus_myst_ggplot

#Add phylopics for groups
PCA_all_genus_myst_ggplot <- 
  PCA_all_genus_myst_ggplot +
  add_phylopic(Bala, alpha = 1, x = -0.22, y = -0.12, ysize = 0.017, color = mypalette_taxa[1])+
  add_phylopic(Ball, alpha = 1, x = 0, y = 0, ysize = 0.012, color = mypalette_taxa[2])+
  add_phylopic(Bals, alpha = 1, x = -0.12, y = 0.07, ysize = 0.014, color = mypalette_taxa[3])+
  add_phylopic(Cap, alpha = 1, x = -0.23, y = 0, ysize = 0.015, color = mypalette_taxa[4])+
  add_phylopic(Esch, alpha = 1, x = -0.09, y = -0.15, ysize = 0.012, color = mypalette_taxa[5])+
  add_phylopic(Meg, alpha = 1, x = -0.03, y = 0.14, ysize = 0.015, color = mypalette_taxa[8])
PCA_all_genus_myst_ggplot

#Get shapes at extremes Mysticeti
#Specify predicted shapes based on position in morphospace (here x/y coordinates)

PC_myst <- PCA_all$x[rows_mysticeti,1:2]
#User-picked spots can be anything, but it in this case evenly-spaced PCA coordinates
preds <- shape.predictor(gdf$coords[,,rows_mysticeti], x= PC_myst, Intercept = FALSE,
                         predPC1_PC2min = c(-0.26, -0.21),
                         predPC1_PC2max = c(0.07, 0.17))

#PC1_PC2max colors
PC1_PC2max_myst_points <- spheres3d(preds$predPC1_PC2max, radius=.002, color = col_modules)

rgl.snapshot(filename = "Output/PC1_PC2max_myst.png") 
rgl.snapshot(filename = "Output/PC1_PC2max_myst1.png") 
clear3d()

#PC1_PC2min colors
PC1_PC2min_myst_points <- spheres3d(preds$predPC1_PC2min, radius=.002, color = col_modules)

rgl.snapshot(filename = "Output/PC1_PC2min_myst.png") 
rgl.snapshot(filename = "Output/PC1_PC2min_myst1.png") 
clear3d()


###Regression PC1 and PC2 ----

#Calculate regression for each component for size
reg_PC1all_size <- lm(Comp1 ~ size, data = pcscores_all)
reg_PC2all_size <- lm(Comp2 ~ size, data = pcscores_all)

#View results and p-value
summary(reg_PC1all_size)
summary(reg_PC2all_size)
anova(reg_PC1all_size)
anova(reg_PC2all_size)

#Save results of significant regression to file
sink("Output/PC1-2all_size_lm.txt")
print("PC1")
summary(reg_PC1all_size)
anova(reg_PC1all_size)
print("PC2")
summary(reg_PC2all_size)
anova(reg_PC2all_size)
sink() 

#Calculate regression for each component taking family into account
reg_PC1all_family <- lm(Comp1 ~ family, data = pcscores_all)
reg_PC2all_family <- lm(Comp2 ~ family, data = pcscores_all)

#View results and p-value
summary(reg_PC1all_family)
summary(reg_PC2all_family)
anova(reg_PC1all_family)
anova(reg_PC2all_family)

#Save results of significant regression to file
sink("Output/PC1-2all_family_lm.txt")
print("PC1")
summary(reg_PC1all_family)
anova(reg_PC1all_family)
print("PC2")
summary(reg_PC2all_family)
anova(reg_PC2all_family)
sink() 

#Calculate regression for each component taking genus into account
reg_PC1all_genus <- lm(Comp1 ~ genus, data = pcscores_all)
reg_PC2all_genus <- lm(Comp2 ~ genus, data = pcscores_all)

#View results and p-value
summary(reg_PC1all_genus)
summary(reg_PC2all_genus)
anova(reg_PC1all_genus)
anova(reg_PC2all_genus)

#Save results of significant regression to file
sink("Output/PC1-2all_genus_lm.txt")
print("PC1")
summary(reg_PC1all_genus)
anova(reg_PC1all_genus)
print("PC2")
summary(reg_PC2all_genus)
anova(reg_PC2all_genus)
sink() 

#Calculate regression for each component taking category into account
reg_PC1all_category <- lm(Comp1 ~ category, data = pcscores_all)
reg_PC2all_category <- lm(Comp2 ~ category, data = pcscores_all)

#View results and p-value
summary(reg_PC1all_category)
summary(reg_PC2all_category)
anova(reg_PC1all_category)
anova(reg_PC2all_category)

#Save results of significant regression to file
sink("Output/PC1-2all_category_lm.txt")
print("PC1")
summary(reg_PC1all_category)
anova(reg_PC1all_category)
print("PC2")
summary(reg_PC2all_category)
anova(reg_PC2all_category)
sink() 

#Calculate regression for each component taking group into account
reg_PC1all_group <- lm(Comp1 ~ group, data = pcscores_all)
reg_PC2all_group <- lm(Comp2 ~ group, data = pcscores_all)

#View results and p-value
summary(reg_PC1all_group)
summary(reg_PC2all_group)
anova(reg_PC1all_group)
anova(reg_PC2all_group)

#Save results of significant regression to file
sink("Output/PC1-2all_group_lm.txt")
print("PC1")
summary(reg_PC1all_group)
anova(reg_PC1all_group)
print("PC2")
summary(reg_PC2all_group)
anova(reg_PC2all_group)
sink() 

#Calculate regression for each component taking feeding mode into account
reg_PC1all_feeding <- lm(Comp1 ~ feeding, data = pcscores_all)
reg_PC2all_feeding <- lm(Comp2 ~ feeding, data = pcscores_all)

#View results and p-value
summary(reg_PC1all_feeding)
summary(reg_PC2all_feeding)
anova(reg_PC1all_feeding)
anova(reg_PC2all_feeding)

#Save results of significant regression to file
sink("Output/PC1-2all_feeding_lm.txt")
print("PC1")
summary(reg_PC1all_feeding)
anova(reg_PC1all_feeding)
print("PC2")
summary(reg_PC2all_feeding)
anova(reg_PC2all_feeding)
sink() 

#Calculate regression for each component taking feeding mode into account - Mysticeti only
reg_PC1all_feeding_myst <- lm(Comp1[rows_mysticeti] ~ feeding[rows_mysticeti], data = pcscores_all)
reg_PC2all_feeding_myst <- lm(Comp2[rows_mysticeti] ~ feeding[rows_mysticeti], data = pcscores_all)

#View results and p-value
summary(reg_PC1all_feeding_myst)
summary(reg_PC2all_feeding_myst)
anova(reg_PC1all_feeding_myst)
anova(reg_PC2all_feeding_myst)

#Save results of significant regression to file
sink("Output/PC1-2all_feeding_myst_lm.txt")
print("PC1")
summary(reg_PC1all_feeding_myst)
anova(reg_PC1all_feeding_myst)
print("PC2")
summary(reg_PC2all_feeding_myst)
anova(reg_PC2all_feeding_myst)
sink() 

#Calculate regression for each component taking family into account - Mysticeti only
reg_PC1all_family_myst <- lm(Comp1[rows_mysticeti] ~ family[rows_mysticeti], data = pcscores_all)
reg_PC2all_family_myst <- lm(Comp2[rows_mysticeti] ~ family[rows_mysticeti], data = pcscores_all)

#View results and p-value
summary(reg_PC1all_family_myst)
summary(reg_PC2all_family_myst)
anova(reg_PC1all_family_myst)
anova(reg_PC2all_family_myst)

#Save results of significant regression to file
sink("Output/PC1-2all_family_myst_lm.txt")
print("PC1")
summary(reg_PC1all_family_myst)
anova(reg_PC1all_family_myst)
print("PC2")
summary(reg_PC2all_family_myst)
anova(reg_PC2all_family_myst)
sink() 

#Calculate regression for each component taking genus into account - Mysticeti only
reg_PC1all_genus_myst <- lm(Comp1[rows_mysticeti] ~ genus[rows_mysticeti], data = pcscores_all)
reg_PC2all_genus_myst <- lm(Comp2[rows_mysticeti] ~ genus[rows_mysticeti], data = pcscores_all)

#View results and p-value
summary(reg_PC1all_genus_myst)
summary(reg_PC2all_genus_myst)
anova(reg_PC1all_genus_myst)
anova(reg_PC2all_genus_myst)

#Save results of significant regression to file
sink("Output/PC1-2all_genus_myst_lm.txt")
print("PC1")
summary(reg_PC1all_genus_myst)
anova(reg_PC1all_genus_myst)
print("PC2")
summary(reg_PC2all_genus_myst)
anova(reg_PC2all_genus_myst)
sink() 

#Calculate regression for each component taking feeding mode and categories into account
reg_PC1all_feeding_cat <- lm(Comp1 ~ feeding * category, data = pcscores_all)
reg_PC2all_feeding_cat <- lm(Comp2 ~ feeding * category, data = pcscores_all)

#View results and p-value
summary(reg_PC1all_feeding_cat)
summary(reg_PC2all_feeding_cat)
anova(reg_PC1all_feeding_cat)
anova(reg_PC2all_feeding_cat)

#Save results of significant regression to file
sink("Output/PC1-2all_feeding_cat_lm.txt")
print("PC1")
summary(reg_PC1all_feeding_cat)
anova(reg_PC1all_feeding_cat)
print("PC2")
summary(reg_PC2all_feeding_cat)
anova(reg_PC2all_feeding_cat)
sink() 

#Save results of all regressions to 1 file
sink("Output/PC1-2_all_lm.txt")
print("PC1")
anova(reg_PC1all_size)
anova(reg_PC1all_family)
anova(reg_PC1all_genus)
anova(reg_PC1all_category)
anova(reg_PC1all_group)
anova(reg_PC1all_genus_myst)
anova(reg_PC1all_family_myst)
anova(reg_PC1all_feeding)
anova(reg_PC1all_feeding_myst)
anova(reg_PC1all_feeding_cat)
print("PC2")
anova(reg_PC2all_size)
anova(reg_PC2all_family)
anova(reg_PC2all_genus)
anova(reg_PC2all_category)
anova(reg_PC2all_group)
anova(reg_PC2all_genus_myst)
anova(reg_PC2all_family_myst)
anova(reg_PC2all_feeding)
anova(reg_PC2all_feeding_myst)
anova(reg_PC2all_feeding_cat)
sink()

###### 
#Next - ch. 4 - ANOVA feeding whole skull