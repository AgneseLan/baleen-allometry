
#===========================================================#
#                                                           #
#       CURVES AND POINTS ANALYSES - MYSTICETI ONLY         #
#                                                           #
#===========================================================#


#CH. 4 - Allometry correction for common allometry, PCA on residuals

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
library(rray)
library(abind)
library(reshape2)
library(scales)
library(mcp)

#require(devtools)
#install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")
#devtools::install_github("kassambara/easyGgplot2")
#remotes::install_github("r-lib/rray")

#ALLOMETRY CORRECTION ----
##Evaluate allometry and get the allometry-free shapes using LogCS, use this for analyses

#Regression shape on logCS size
allometry <- procD.lm(gdf$coords ~ gdf$size, turbo = F, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of allometry with logCS
summary(allometry) 

#Regression score of shape vs logCS - regression method with "RegScore" plotting
allometry_plot_regscore <- plot(allometry, type = "regression",predictor = gdf$size, reg.type = "RegScore",
                                main = "Shape vs logCS",xlab = "logCS", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics
text(x = gdf$size, y = allometry_plot_regscore$RegScore, labels = gdf$Id,
     pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels

##Add regression line with confidence intervals to plot
#Create object to use for linear model
allometry_regscores <- allometry_plot_regscore[["RegScore"]] 

#Linear model for line
allometry_regline <- lm(allometry_regscores ~ gdf$size)

##Make better allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_plot <- data.frame(logCS = allometry_plot_regscore[["plot.args"]][["x"]], RegScores = allometry_plot_regscore[["plot.args"]][["y"]])

#Convert data frame to tibble
allometry_plot <- as_tibble(allometry_plot)
#Add labels and other attributes to tibble as columns
allometry_plot <- allometry_plot %>% 
  mutate(specimens = gdf$Id, family = gdf$family, category = gdf$category, group = gdf$group, stage = gdf$stage,
         genus = gdf$genus, feeding = gdf$feeding, size = gdf$size, TL = gdf$TL_100, BZW = gdf$BZW_100)
glimpse(allometry_plot)

#Nice plot with specimens colored by age AND regression line with confidence intervals
allometry_ggplot <- ggplot(allometry_plot, aes(x = logCS, y = RegScores, label = specimens))+
  geom_point(aes(colour = category, fill = category, shape = group), size = 2, alpha = 0.8)+   
  geom_smooth(aes(x = size, y = RegScores), method = 'lm', inherit.aes = F,         #confidence intervals and reg line, before points
              colour = "darkblue", fill = 'gainsboro', linetype = "dashed", size = 0.8)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_colour_manual(name = "Growth stage", labels = levels(categories), 
                      values = mypalette_category, aesthetics = c("colour", "fill"))+      
  scale_shape_manual(values = shapes, name = "Group", labels = levels(groups))+
  theme_classic(base_size = 12)+
  ylab("Regression Score - p=0.001**")+
  theme(legend.key = element_blank(), legend.title = element_text(size = 11, face = "bold"), 
        legend.position = c(0.8,0.2),  legend.direction = "vertical", legend.justification = c(0,0))+
  geom_text_repel(colour = "black", size = 3,          #label last so that they are on top of fill
                  force_pull = 3, point.padding = 1, max.overlaps = 40)     #position of tables relative to point (proximity and distance) 
allometry_ggplot

##Create residuals array from null to then save as coordinates for analyses
allometry_array <- arrayspecs(allometry$residuals,p = dim(gdf$coords)[1], k = dim(gdf$coords)[2]) 

#New shapes adjusted for allometry with CS to use in analyses
allometry_residuals <- allometry_array + array(mean_shape, dim(allometry_array)) 

#Save mean shape of allometry-adjusted shapes to use later
mean_shape_residuals <- mshape(allometry_residuals)

####Gdf for allometry residuals ----
gdf_res <- geomorph.data.frame(coords = allometry_residuals,  Id = classifiers$code, genus = classifiers$genus2, family = classifiers$family, 
                               group = classifiers$group, category = classifiers$category, stage = gdf$stage,
                               TL_100 = classifiers$TL_100, BZW_100 = classifiers$BZW_100,
                               feeding = classifiers$Feeding_BL20)

#PCA ALLOMETRY RESIDUALS ----

#New PCA plot with data corrected for allometry
PCA_residuals <- gm.prcomp(allometry_residuals) 

#List of PC components and proportion of variations
PCA_residuals

#Save PCA results to file
sink("Output/PCA_residuals_components.txt")
print("PCA allometry residuals")
print(PCA_residuals)
sink() 

#Change row names to codes to make plot readable
row.names(PCA_residuals$x) <- gdf_res$Id

##View plot
plot(PCA_residuals, main = "PCA residuals - PC1-PC2",  pch = 21, #title and type of point to be used
     col = "deeppink",    bg = "deeppink",  cex = 1, font.main = 2)      #improve graphics
#Add quick labels to plot
text(x = PCA_residuals$x[,1], y = PCA_residuals$x[,2], labels = rownames(PCA_residuals$x), 
     pos = 1,   offset = 0.5,  cex = 0.75)    #improve graphics 

##View plot
plot(PCA_residuals, axis1 = 1, axis2 = 3, main = "PCA residuals - PC1-PC3",  pch = 21, #title and type of point to be used
     col = "deeppink",    bg = "deeppink",  cex = 1, font.main = 2)      #improve graphics
#Add quick labels to plot
text(x = PCA_residuals$x[,1], y = PCA_residuals$x[,3], labels = rownames(PCA_residuals$x), 
     pos = 1,   offset = 0.5,  cex = 0.75)    #improve graphics 

#Save PC scores as object to use later
pcscores_res <- PCA_residuals$x

#Min max shapes code
#Save shapes of extremes for axes used in plot
PC1min_res <- PCA_residuals[["shapes"]][["shapes.comp1"]][["min"]]
PC1max_res <- PCA_residuals[["shapes"]][["shapes.comp1"]][["max"]] 
PC2min_res <- PCA_residuals[["shapes"]][["shapes.comp2"]][["min"]] 
PC2max_res <- PCA_residuals[["shapes"]][["shapes.comp2"]][["max"]] 

#Show 3D deformation from mean with points overlay, do this for all 4 extremes - using spheres3D for points
#PC1min colors
#spheres3d(mean_shape, radius=.0005, color = "gray60", alpha = 0.5, fastTransparency = T) - plot mean specimens with transparency
PC1min_res_points <- spheres3d(PC1min_res, radius=.001, color = col_modules)

rgl.snapshot(filename = "Output/PC1min_res.png") 
rgl.snapshot(filename = "Output/PC1min_res1.png")
clear3d()

#PC1max colors
PC1max_res_points <- spheres3d(PC1max_res, radius=.001, color = col_modules)

rgl.snapshot(filename = "Output/PC1max_res.png") 
rgl.snapshot(filename = "Output/PC1max_res1.png")
clear3d()

#PC2min colors
PC2min_res_points <- spheres3d(PC2min_res, radius=.001, color = col_modules)

rgl.snapshot(filename = "Output/PC2min_res.png") 
rgl.snapshot(filename = "Output/PC2min_res1.png") 
clear3d()

#PC2max colors
PC2max_res_points <- spheres3d(PC2max_res, radius=.001, color = col_modules)

rgl.snapshot(filename = "Output/PC2max_res.png") 
rgl.snapshot(filename = "Output/PC2max_res1.png") 
clear3d()

##Make better PCA plot using ggplot
#Read PC scores as tibble
pcscores_res <- as_tibble(pcscores_res)
#Add labels and other attributes to tibble as columns
pcscores_res <- pcscores_res %>% 
  mutate(specimens = gdf_res$Id, group = gdf_res$group, category = gdf_res$category, family = gdf_res$family,
         genus = gdf_res$genus, feeding = gdf_res$feeding)
glimpse(pcscores_res)

#Nice PCA plot with stages and groups
PCA_res_ggplot <- ggplot(pcscores_res,  aes(x = Comp1, y = Comp2, label = specimens, colour = genus, fill = genus))+
  geom_point(size = 3, aes(shape = category))+
  geom_text_repel(colour = "black", size = 4, max.overlaps = 40)+
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = shapes_cat)+            #legend and color adjustments
  scale_colour_manual(name = "Genus", labels = levels(genera), #copy from as.factor(genera)
                      values = mypalette_taxa, aesthetics = c("colour","fill"), 
                      guide = guide_legend(label.theme = element_text(size =10, angle = 0, face = "italic")))+
  theme_bw()+
  xlab("PC 1 (50.8%)")+ #copy this from standard PCA plot (PCA_all_plot)
  ylab("PC 2 (11.9%)")+
  theme(legend.title = element_text(size = 12, face = "bold")) 

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_ggplot

#Make hulls for PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
#Divide by group first
hulls_res_category_myst <- pcscores_res %>%
  filter(group %in% c("mysticeti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

hulls_res_category_odont <- pcscores_res %>%
  filter(group %in% c("odontoceti")) %>%
  group_by(category) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around categories (4 stages: Early Fetus, Late Fetus/Neonate, Juvenile, Adult)
PCA_res_category_ggplot <- ggplot(pcscores_res, aes(x = Comp1, y = Comp2, colour = category, fill = category))+
  geom_point(size = 3, aes(shape = family))+
  scale_colour_manual(name = "Growth stage", labels = levels(categories), #to be ordered as they appear in tibble
                      values = mypalette_category)+            #legend and color adjustments
  geom_polygon(data = hulls_res_category_myst, aes(x = x, y = y, fill = category), 
               alpha = .1, size = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  geom_polygon(data = hulls_res_category_odont, aes(x = x, y = y, fill = category), 
               alpha = .1, size = 0.6, show.legend = FALSE)+ #colored hulls with transparency
  scale_fill_manual(name = "Growth stage", labels = levels(categories),
                    values =  mypalette_category)+ #must match scale_colour_manual
  scale_shape_manual(name = "Family", labels = levels(families),
                     values = shapes_fam)+
  theme_bw()+
  xlab("PC 1 (50.8%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (11.9%)")+
  theme(legend.title = element_text(size = 11, face = "bold"))+
  #Remove legend for a scale_ using guide
  guides(fill = guide_legend(label = F, title = NULL, override.aes = list(shape = NA)))

#scale_y_reverse() #reverse y scale to match traj analysis

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_category_ggplot

#Add phylopics for groups
PCA_res_category_ggplot <- 
  PCA_res_category_ggplot +
  add_phylopic(myst, alpha = 1, x = -0.05, y = 0.1, ysize = 0.015, color = "gray30")+
  add_phylopic(odont, alpha = 1, x = 0.2, y = -0.05, ysize = 0.025, color = "gray30")
PCA_res_category_ggplot

#Make hulls for PCA plot with hulls around genera
hulls_res_genus <- pcscores_res %>%
  group_by(genus) %>%
  slice(chull(Comp1, Comp2)) %>%
  rename(x = Comp1, y = Comp2)

#Nice PCA plot with hulls around genera 
PCA_res_genus_ggplot1 <- ggplot(pcscores_res, aes(x = Comp1, y = Comp2))+
  geom_polygon(data = hulls_res_genus, aes(x = x, y = y, group = genus, colour = genus, linetype = group), inherit.aes = F,
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
  xlab("PC 1 (50.8%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (11.9%)")+
  theme(legend.title = element_text(size = 11, face = "bold"))+
  guides(shape = guide_legend(override.aes = list(fill = "grey10", colour = "grey10")))

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_genus_ggplot1

##PCA mysticeti only----
#Nice PCA plot with hulls around genera - myst only
pcscores_res_myst <- pcscores_res %>% filter(group == "mysticeti")
hulls_res_genus_myst <- hulls_res_genus %>% filter(group == "mysticeti")

PCA_res_genus_myst_ggplot <- ggplot(pcscores_res_myst, aes(x = Comp1, y = Comp2))+
  geom_polygon(data = hulls_res_genus_myst, aes(x = x, y = y, group = genus, colour = genus, linetype = family), inherit.aes = F,
               linewidth = 1, alpha = 0.005, show.legend = FALSE)+ #colored hulls with transparency
  geom_point(size = 4, aes(shape = category, alpha = category, colour = genus, fill = genus))+
  scale_colour_manual(name = "Genus", labels = levels(as.factor(pcscores_res_myst$genus)), 
                      values = mypalette_taxa[-c(6,7,9)], aesthetics = c("colour","fill"),
                      guide = guide_legend(label.theme = element_text(size =10, angle = 0, face = "italic")))+ #to be ordered as they appear in tibble
  scale_alpha_manual(name = "Growth stage", labels =  levels(categories), #to be ordered as they appear in tibble
                     values = c(0.3, 0.5, 0.7, 1))+            #legend and color adjustments
  scale_shape_manual(name = "Growth stage", labels =  levels(categories), 
                     values = shapes_cat)+
  scale_linetype_manual(values = c(3,1,2))+
  theme_bw()+
  xlab("PC 1 (50.8%)")+ #copy this from standard PCA plot (PCA_res_plot)
  ylab("PC 2 (11.9%)")+
  theme(legend.title = element_text(size = 11, face = "bold"), 
        legend.position = c(0.85,0.2), legend.key = element_blank(),
        legend.background=element_blank())+
  guides(shape = guide_legend(override.aes = list(fill = "grey10", colour = "grey10")),
         color = "none", fill = "none")

#Visualize plot and save as PDF using menu in bar on the right
PCA_res_genus_myst_ggplot

#Add phylopics for groups
PCA_res_genus_myst_ggplot <- 
  PCA_res_genus_myst_ggplot +
  add_phylopic(Bala, alpha = 1, x = -0.05, y = 0.07, ysize = 0.009, color = mypalette_taxa[1])+
  add_phylopic(Ball, alpha = 1, x = 0.028, y = 0, ysize = 0.007, color = mypalette_taxa[2])+
  add_phylopic(Bals, alpha = 1, x = -0.07, y = -0.075, ysize = 0.007, color = mypalette_taxa[3])+
  add_phylopic(Cap, alpha = 1, x = -0.15, y = 0.025, ysize = 0.008, color = mypalette_taxa[4])+
  add_phylopic(Esch, alpha = 1, x = 0.025, y = 0.06, ysize = 0.006, color = mypalette_taxa[5])+
  add_phylopic(Meg, alpha = 1, x = -0.01, y = -0.05, ysize = 0.008, color = mypalette_taxa[8])
PCA_res_genus_myst_ggplot

#Get shapes at extremes Mysticeti
#Specify predicted shapes based on position in morphospace (here x/y coordinates)

PC_myst_res <- PCA_residuals$x[rows_mysticeti,1:2]
#User-picked spots can be anything, but it in this case evenly-spaced PCA coordinates
#User-picked spots can be anything, but it in this case evenly-spaced PCA coordinates
preds_res <- shape.predictor(gdf_res$coords[,,rows_mysticeti], x= PC_myst_res, Intercept = FALSE,
                         predPC1_PC2min = c(-0.16, -0.08),
                         predPC1_PC2max = c(0.08, 0.08))

#PC1_PC2max colors
PC1_PC2max_myst_res_points <- spheres3d(preds_res$predPC1_PC2max, radius=.0004, color = col_modules)

rgl.snapshot(filename = "Output/PC1_PC2max_myst_res.png") 
rgl.snapshot(filename = "Output/PC1_PC2max_myst_res1.png") 
clear3d()

#PC1_PC2min colors
PC1_PC2min_myst_res_points <- spheres3d(preds_res$predPC1_PC2min, radius=.003, color = col_modules)

rgl.snapshot(filename = "Output/PC1_PC2min_myst_res.png") 
rgl.snapshot(filename = "Output/PC1_PC2min_myst_res1.png") 
clear3d()


###Regression PC1 and PC2 ----
#Calculate regression for each component taking family into account
reg_PC1res_genus <- lm(Comp1 ~ genus, data = pcscores_res)
reg_PC2res_genus <- lm(Comp2 ~ genus, data = pcscores_res)

#View results and p-value
summary(reg_PC1res_genus)
summary(reg_PC2res_genus)
anova(reg_PC1res_genus)
anova(reg_PC2res_genus)

#Save results of significant regression to file
sink("Output/PC1-2res_genus_lm.txt")
print("PC1")
summary(reg_PC1res_genus)
anova(reg_PC1res_genus)
print("PC2")
summary(reg_PC2res_genus)
anova(reg_PC2res_genus)
sink() 


#Calculate regression for each component taking family into account
reg_PC1res_family <- lm(Comp1 ~ family, data = pcscores_res)
reg_PC2res_family <- lm(Comp2 ~ family, data = pcscores_res)

#View results and p-value
summary(reg_PC1res_family)
summary(reg_PC2res_family)
anova(reg_PC1res_family)
anova(reg_PC2res_family)

#Save results of significant regression to file
sink("Output/PC1-2res_family_lm.txt")
print("PC1")
summary(reg_PC1res_family)
anova(reg_PC1res_family)
print("PC2")
summary(reg_PC2res_family)
anova(reg_PC2res_family)
sink() 

#Calculate regression for each component taking category into account
reg_PC1res_category <- lm(Comp1 ~ category, data = pcscores_res)
reg_PC2res_category <- lm(Comp2 ~ category, data = pcscores_res)

#View results and p-value
summary(reg_PC1res_category)
summary(reg_PC2res_category)
anova(reg_PC1res_category)
anova(reg_PC2res_category)

#Save results of significant regression to file
sink("Output/PC1-2res_category_lm.txt")
print("PC1")
summary(reg_PC1res_category)
anova(reg_PC1res_category)
print("PC2")
summary(reg_PC2res_category)
anova(reg_PC2res_category)
sink() 

#Calculate regression for each component taking group into account
reg_PC1res_group <- lm(Comp1 ~ group, data = pcscores_res)
reg_PC2res_group <- lm(Comp2 ~ group, data = pcscores_res)

#View results and p-value
summary(reg_PC1res_group)
summary(reg_PC2res_group)
anova(reg_PC1res_group)
anova(reg_PC2res_group)

#Save results of significant regression to file
sink("Output/PC1-2res_group_lm.txt")
print("PC1")
summary(reg_PC1res_group)
anova(reg_PC1res_group)
print("PC2")
summary(reg_PC2res_group)
anova(reg_PC2res_group)
sink() 

#Calculate regression for each component taking feeding mode into account
reg_PC1res_feeding <- lm(Comp1 ~ feeding, data = pcscores_res)
reg_PC2res_feeding <- lm(Comp2 ~ feeding, data = pcscores_res)

#View results and p-value
summary(reg_PC1res_feeding)
summary(reg_PC2res_feeding)
anova(reg_PC1res_feeding)
anova(reg_PC2res_feeding)

#Save results of significant regression to file
sink("Output/PC1-2res_feeding_lm.txt")
print("PC1")
summary(reg_PC1res_feeding)
anova(reg_PC1res_feeding)
print("PC2")
summary(reg_PC2res_feeding)
anova(reg_PC2res_feeding)
sink() 

#Calculate regression for each component taking feeding mode into account - Mysticeti only
reg_PC1res_feeding_myst <- lm(Comp1[rows_mysticeti] ~ feeding[rows_mysticeti], data = pcscores_res)
reg_PC2res_feeding_myst <- lm(Comp2[rows_mysticeti] ~ feeding[rows_mysticeti], data = pcscores_res)

#View results and p-value
summary(reg_PC1res_feeding_myst)
summary(reg_PC2res_feeding_myst)
anova(reg_PC1res_feeding_myst)
anova(reg_PC2res_feeding_myst)

#Save results of significant regression to file
sink("Output/PC1-2res_feeding_myst_lm.txt")
print("PC1")
summary(reg_PC1res_feeding_myst)
anova(reg_PC1res_feeding_myst)
print("PC2")
summary(reg_PC2res_feeding_myst)
anova(reg_PC2res_feeding_myst)
sink() 

#Calculate regression for each component taking family into account - Mysticeti only
reg_PC1res_family_myst <- lm(Comp1[rows_mysticeti] ~ family[rows_mysticeti], data = pcscores_res)
reg_PC2res_family_myst <- lm(Comp2[rows_mysticeti] ~ family[rows_mysticeti], data = pcscores_res)

#View results and p-value
summary(reg_PC1res_family_myst)
summary(reg_PC2res_family_myst)
anova(reg_PC1res_family_myst)
anova(reg_PC2res_family_myst)

#Save results of significant regression to file
sink("Output/PC1-2res_family_myst_lm.txt")
print("PC1")
summary(reg_PC1res_family_myst)
anova(reg_PC1res_family_myst)
print("PC2")
summary(reg_PC2res_family_myst)
anova(reg_PC2res_family_myst)
sink() 

#Calculate regression for each component taking genus into account - Mysticeti only
reg_PC1res_genus_myst <- lm(Comp1[rows_mysticeti] ~ genus[rows_mysticeti], data = pcscores_res)
reg_PC2res_genus_myst <- lm(Comp2[rows_mysticeti] ~ genus[rows_mysticeti], data = pcscores_res)

#View results and p-value
summary(reg_PC1res_genus_myst)
summary(reg_PC2res_genus_myst)
anova(reg_PC1res_genus_myst)
anova(reg_PC2res_genus_myst)

#Save results of significant regression to file
sink("Output/PC1-2res_genus_myst_lm.txt")
print("PC1")
summary(reg_PC1res_genus_myst)
anova(reg_PC1res_genus_myst)
print("PC2")
summary(reg_PC2res_genus_myst)
anova(reg_PC2res_genus_myst)
sink() 

#Calculate regression for each component taking feeding mode and categories into account
reg_PC1res_feeding_cat <- lm(Comp1 ~ feeding * category, data = pcscores_res)
reg_PC2res_feeding_cat <- lm(Comp2 ~ feeding * category, data = pcscores_res)

#View results and p-value
summary(reg_PC1res_feeding_cat)
summary(reg_PC2res_feeding_cat)
anova(reg_PC1res_feeding_cat)
anova(reg_PC2res_feeding_cat)

#Save results of significant regression to file
sink("Output/PC1-2res_feeding_cat_lm.txt")
print("PC1")
summary(reg_PC1res_feeding_cat)
anova(reg_PC1res_feeding_cat)
print("PC2")
summary(reg_PC2res_feeding_cat)
anova(reg_PC2res_feeding_cat)
sink() 

#Save results of res regressions to 1 file
#Save results of res regressions to 1 file
sink("Output/PC1-2_res_lm.txt")
print("PC1")
anova(reg_PC1res_family)
anova(reg_PC1res_genus)
anova(reg_PC1res_category)
anova(reg_PC1res_group)
anova(reg_PC1res_genus_myst)
anova(reg_PC1res_family_myst)
anova(reg_PC1res_feeding)
anova(reg_PC1res_feeding_myst)
anova(reg_PC1res_feeding_cat)
print("PC2")
anova(reg_PC2res_family)
anova(reg_PC2res_genus)
anova(reg_PC2res_category)
anova(reg_PC2res_group)
anova(reg_PC2res_genus_myst)
anova(reg_PC2res_family_myst)
anova(reg_PC2res_feeding)
anova(reg_PC2res_feeding_myst)
anova(reg_PC2res_feeding_cat)
sink()

###### 
#Next - ch. 5 - Allometry analysis by genus
