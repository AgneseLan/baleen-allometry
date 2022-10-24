
#===========================================================#
#                                                           #
#       CURVES AND POINTS ANALYSES - MYSTICETI ONLY         #
#                                                           #
#===========================================================#

#CH. 7 - Mophological disparity and clustering analyses

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
library(reshape2)
library(scales)
require(grid)
library(Anthropometry)

#require(devtools)
#install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")
#devtools::install_github("kassambara/ggcorrplot")
#apropos("x") lists objects with matching part of name

#MORPHOLOGICAL DISPARITY ----
#Calculate Procrustes variances and distances between groups, with p-value for each pair of groups
#How different are each group shapes compared to other group shapes?

#Since most genera do not have all categories, use stages for feeding
#Do not use for genera, too few as well

#All data
#Disparity between categories, considering entire dataset and not different genera
disparity_category <- morphol.disparity(coords ~ 1, groups = ~ category, iter = 999, data = gdf)

#Results and significance
summary(disparity_category)

#Save results to file
sink("Output/disparity_category.txt")
print(summary(disparity_category))
sink() 

#Disparity between genera, considering entire dataset and not different categories
disparity_genus <- morphol.disparity(coords ~ 1, groups = ~ genus, iter = 999, data = gdf)

#Results and significance
summary(disparity_genus)

#Save results to file
sink("Output/disparity_genus.txt")
print(summary(disparity_genus))
sink() 

#Disparity between groups, considering entire dataset and not different categories
disparity_group <- morphol.disparity(coords ~ 1, groups = ~ group, iter = 999, data = gdf)

#Results and significance
summary(disparity_group)

#Save results to file
sink("Output/disparity_group.txt")
print(summary(disparity_group))
sink() 

#Disparity between categories, considering entire dataset AND groups
disparity_group_category <- morphol.disparity(coords ~ 1, groups = ~ group*category, iter = 999, data = gdf)

#Results and significance
summary(disparity_group_category)

#Save results to file
sink("Output/disparity_group_category.txt")
print(summary(disparity_group_category))
sink() 

#Allometry corrected
#Disparity between categories with size, considering entire dataset and not different groups
disparity_size_category <- morphol.disparity(coords ~ size, groups = ~ category, iter = 999, data = gdf)

#Results and significance
summary(disparity_size_category)

#Save results to file
sink("Output/disparity_size_category.txt")
print(summary(disparity_size_category))
sink() 

#Disparity between genera with size, considering entire dataset and not different categories
disparity_size_genus <- morphol.disparity(coords ~ size, groups = ~ genus, iter = 999, data = gdf)

#Results and significance
summary(disparity_size_genus)

#Save results to file
sink("Output/disparity_size_genus.txt")
print(summary(disparity_size_genus))
sink() 

#Disparity between groups with size, considering entire dataset and not different categories
disparity_size_group <- morphol.disparity(coords ~ size, groups = ~ group, iter = 999, data = gdf)

#Results and significance
summary(disparity_size_group)

#Save results to file
sink("Output/disparity_size_group.txt")
print(summary(disparity_size_group))
sink() 

#Disparity between categories with size, considering entire dataset AND groups
disparity_size_group_category <- morphol.disparity(coords ~ size, groups = ~ group*category, iter = 999, data = gdf)

#Results and significance
summary(disparity_size_group_category)

#Save results to file
sink("Output/disparity_size_group_category.txt")
print(summary(disparity_size_group_category))
sink() 

##Heatmaps plots for significant differences in disparity ----

#Create palette for heatmap trajectory plot
mypalette_disp <- brewer.pal(9,"PuRd")
image(1:9,1, as.matrix(1:9), col = mypalette_disp,xlab="Reds (sequential)",
      ylab = "", yaxt = "n")

###All data by group and category ----
#Save p-values as object
disp_corr <- disparity_group_category[["PV.dist"]]
disp_pvals <- disparity_group_category[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
vars <- rownames(disp_corr)

#Replace string names to make them shorter
vars <- str_replace_all(vars, "\\.", "_")

categories_list_short <- c("1", "2", "3", "4")

#Loop replacements categories
for (u in 1:length(categories_list)){
  vars <- str_replace_all(vars, categories_list[u], categories_list_short[u])
}

#Loop replacements genera
for (t in 1:length(groups_list)){
  vars <- str_replace_all(vars, groups_list[t], groups_list_short[t])
}

#Check it worked
vars

#Set correct row and col names for both
rownames(disp_corr) <- vars
rownames(disp_pvals) <- vars
colnames(disp_corr) <- vars
colnames(disp_pvals) <- vars

#Get upper triangles only - half matrix, eliminates redundant info
disp_corr_upper_tri <- get_upper_tri(disp_corr)
disp_pvals_upper_tri <- get_upper_tri(disp_pvals)

#Melt to make table in the format needed for heatmap
disp_corr_melt <- melt(disp_corr_upper_tri, value.name = "corr", na.rm = TRUE)
disp_pvals_melt <- melt(disp_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
disp_pvals_melt$corr <- disp_corr_melt$corr

#Create columns where only significant values are shown
disp_pvals_melt <- disp_pvals_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                              p_if_sig = ifelse(sig_p, p, NA),
                                              corr_if_sig = ifelse(sig_p, corr, NA))
disp_pvals_melt

#Nice heatmap plot
disparity_group_category_heatmap_ggplot <- ggplot(data = disp_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  scale_fill_gradient2(low = mypalette_disp[9], high = mypalette_disp[2], mid = mypalette_disp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_disp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Morphological disparity by group and growth stage")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
disparity_group_category_heatmap_ggplot

###Mysticeti by genus ----
#Save p-values as object
disp_corr_genus <- disparity_genus[["PV.dist"]]
disp_pvals_genus <- disparity_genus[["PV.dist.Pval"]]

#Get upper triangles only - half matrix, eliminates redundant info
disp_corr_genus_upper_tri <- get_upper_tri(disp_corr_genus)
disp_pvals_genus_upper_tri <- get_upper_tri(disp_pvals_genus)

#Melt to make table in the format needed for heatmap
disp_corr_genus_melt <- melt(disp_corr_genus_upper_tri, value.name = "corr", na.rm = TRUE)
disp_pvals_genus_melt <- melt(disp_pvals_genus_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
disp_pvals_genus_melt$corr <- disp_corr_genus_melt$corr

#Create columns where only significant values are shown
disp_pvals_genus_melt <- disp_pvals_genus_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                          p_if_sig = ifelse(sig_p, p, NA),
                                                          corr_if_sig = ifelse(sig_p, corr, NA))
disp_pvals_genus_melt

#Select mysticeti only
genera_list_mysticeti <- levels(as.factor(gdf$genus[rows_mysticeti]))

disp_pvals_genus_melt_myst <- disp_pvals_genus_melt %>% filter(Var1 %in% genera_list_mysticeti) %>% 
  filter(Var2 %in% genera_list_mysticeti)

#Nice heatmap plot
disparity_genus_heatmap_ggplot <- ggplot(data = disp_pvals_genus_melt_myst, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  scale_fill_gradient2(low = mypalette_disp[9], high = mypalette_disp[2], mid = mypalette_disp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_disp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Morphological disparity by genus in Mysticeti")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
disparity_genus_heatmap_ggplot


###Allometry corrected by group and category ----
#Save p-values as object
disp_size_corr <- disparity_size_group_category[["PV.dist"]]
disp_size_pvals <- disparity_size_group_category[["PV.dist.Pval"]]

#Save row and col names as variables to change string - colnames = rownames for both
vars_size <- rownames(disp_size_corr)

#Replace string names to make them shorter
vars_size <- str_replace_all(vars_size, "\\.", "_")

#Loop replacements categories
for (u in 1:length(categories_list)){
  vars_size <- str_replace_all(vars_size, categories_list[u], categories_list_short[u])
}

#Loop replacements genera
for (t in 1:length(groups_list)){
  vars_size <- str_replace_all(vars_size, groups_list[t], groups_list_short[t])
}

#Check it worked
vars_size

#Set correct row and col names for both
rownames(disp_size_corr) <- vars_size
rownames(disp_size_pvals) <- vars_size
colnames(disp_size_corr) <- vars_size
colnames(disp_size_pvals) <- vars_size

#Get upper triangles only - half matrix, eliminates redundant info
disp_size_corr_upper_tri <- get_upper_tri(disp_size_corr)
disp_size_pvals_upper_tri <- get_upper_tri(disp_size_pvals)

#Melt to make table in the format needed for heatmap
disp_size_corr_melt <- melt(disp_size_corr_upper_tri, value.name = "corr", na.rm = TRUE)
disp_size_pvals_melt <- melt(disp_size_pvals_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
disp_size_pvals_melt$corr <- disp_size_corr_melt$corr

#Create columns where only significant values are shown
disp_size_pvals_melt <- disp_size_pvals_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                        p_if_sig = ifelse(sig_p, p, NA),
                                                        corr_if_sig = ifelse(sig_p, corr, NA))
disp_size_pvals_melt

#Nice heatmap plot
disparity_size_group_category_heatmap_ggplot <- ggplot(data = disp_size_pvals_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  scale_fill_gradient2(low = mypalette_disp[9], high = mypalette_disp[2], mid = mypalette_disp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_disp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Morphological disparity size corr. by group and growth stage")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
disparity_size_group_category_heatmap_ggplot

###Allometry corrected Mysticeti by genus ----
#Save p-values as object
disp_size_corr_genus <- disparity_size_genus[["PV.dist"]]
disp_size_pvals_genus <- disparity_size_genus[["PV.dist.Pval"]]

#Get upper triangles only - half matrix, eliminates redundant info
disp_size_corr_genus_upper_tri <- get_upper_tri(disp_size_corr_genus)
disp_size_pvals_genus_upper_tri <- get_upper_tri(disp_size_pvals_genus)

#Melt to make table in the format needed for heatmap
disp_size_corr_genus_melt <- melt(disp_size_corr_genus_upper_tri, value.name = "corr", na.rm = TRUE)
disp_size_pvals_genus_melt <- melt(disp_size_pvals_genus_upper_tri, value.name = "p", na.rm = TRUE)

#Add column to main table
disp_size_pvals_genus_melt$corr <- disp_size_corr_genus_melt$corr

#Create columns where only significant values are shown
disp_size_pvals_genus_melt <- disp_size_pvals_genus_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                    p_if_sig = ifelse(sig_p, p, NA),
                                                                    corr_if_sig = ifelse(sig_p, corr, NA))
disp_size_pvals_genus_melt

#Select mysticeti only
disp_size_pvals_genus_melt_myst <- disp_size_pvals_genus_melt %>% filter(Var1 %in% genera_list_mysticeti) %>% 
  filter(Var2 %in% genera_list_mysticeti)

#Nice heatmap plot
disparity_size_genus_heatmap_ggplot <- ggplot(data = disp_size_pvals_genus_melt_myst, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  scale_fill_gradient2(low = mypalette_disp[9], high = mypalette_disp[2], mid = mypalette_disp[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_disp[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  ggtitle ("Morphological disparity size corr. by genus in Mysticeti")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(size = 10),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.3,
                               title.position = "top", title.hjust = 0.5))
disparity_size_genus_heatmap_ggplot


#CLUSTERING SHAPES ----
#Test if fetal specimens closer to adults of species or other fetal specimens

#Whole skull
gpa2 <- gpagen(shape_array, verbose = T)
#Align mysticeti separately to look at clustering among them
gpa2_myst <- gpagen(shape_array[,,rows_mysticeti], verbose = T)

#Save proc distances
proc_dist <- gpa2$procD
proc_dist_myst <- gpa2_myst$procD

#Clustering ward
cluster_pDist_ward <- hclust(proc_dist, method = "ward.D2")

cluster_pDist_ward_myst <- hclust(proc_dist_myst, method = "ward.D2")

#All specs
plot(cluster_pDist_ward, labels = Ids, main = "Cluster Proc. Distances",
     xlab = NULL, ylab = "ProcDist") #labels specimens
plot(cluster_pDist_ward, labels = categories, main = "Cluster Proc. Distances",
     xlab = NULL, ylab = "ProcDist") #labels categories
rect.hclust(cluster_pDist_ward, k = 10)

clusters_all_ward <- rect.hclust(cluster_pDist_ward, k = 10)

#Myst only
plot(cluster_pDist_ward_myst, labels = Ids[rows_mysticeti], main = "Cluster Proc. Distances",
     xlab = NULL, ylab = "ProcDist")
plot(cluster_pDist_ward_myst, labels = categories[rows_mysticeti], main = "Cluster Proc. Distances",
     xlab = NULL, ylab = "ProcDist")
rect.hclust(cluster_pDist_ward_myst, k = 5)

clusters_myst_ward_myst <- rect.hclust(cluster_pDist_ward_myst, k = 5)

##Calculate kmeans grouping using shape data ----
#Use k optimized using clustering plots

#All specs
shape_clusters <- LloydShapes(gpa2$coords, numClust = 10, simul = F, verbose = T)

asig <- shape_clusters$asig 
table(asig) 

#Plot with kmeans cluster labels
plot(cluster_pDist_ward, labels = asig, main = "Cluster Proc. Distances",
     xlab = NULL, ylab = "ProcDist")
rect.hclust(cluster_pDist_ward, k = 10)

#Myst only
shape_clusters_myst <- LloydShapes(gpa2_myst$coords, numClust = 5, simul = F, verbose = T)

asig_myst <- shape_clusters_myst$asig
table(asig_myst)

#Plot with kmeans cluster labels
plot(cluster_pDist_ward_myst, labels = asig_myst, main = "Cluster Proc. Distances",
     xlab = NULL, ylab = "ProcDist")
rect.hclust(cluster_pDist_ward_myst, k = 5)

#Save results of significant regression to file
sink("Output/cluster_analyses.txt")
print("All specimens - k=10")
print("ward")
clusters_all_ward
print("kmeans shape")
table(asig)
setNames(asig, Ids)

print("Mysticeti only - k=5")
print("ward")
clusters_myst_ward
print("kmeans shape")
table(asig_myst)
setNames(asig_myst, Ids[rows_mysticeti])
sink() 

