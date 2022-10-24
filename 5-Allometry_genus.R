
#===========================================================#
#                                                           #
#       CURVES AND POINTS ANALYSES - MYSTICETI ONLY         #
#                                                           #
#===========================================================#


#CH. 5 - Allometry by genus, slopes for ASR of allometry

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

##Test different allometry between genera ----
#Divide by group only for plots
allometry_genus_comb <-  procD.lm(gdf$coords ~ gdf$size + gdf$genus, turbo = F, iter=999, print.progress = TRUE) 
allometry_genus_int <-  procD.lm(gdf$coords ~ gdf$size * gdf$genus, turbo = F, iter=999, print.progress = TRUE) 

#Main results of ANOVA analysis of allometry with logCS
summary(allometry_genus_comb)
summary(allometry_genus_int) 

#Save results of significant regression to file
sink("Output/allometry_models_genus.txt")
print("Null")
summary(allometry)

print("Combination +")
summary(allometry_genus_comb) 

print("Interaction *")
summary(allometry_genus_int)
sink() 

###Pairwise test difference allometry by genus ----

#ANOVAs - is a model significantly better than the others?
anova_allometry_models_genus <- anova(allometry, allometry_genus_comb, allometry_genus_int)
anova_allometry_models_genus

#Pairwise comparison for the combination and interaction model
#Helps determine if there is a significant difference in slope (int model) in the allometry trajectory on top of difference in intercept (comb model)
pairwise_allometry_genus <- pairwise(allometry_genus_int, fit.null = allometry_genus_comb,
                                     groups = gdf$genus, 
                                     covariate = gdf$size, print.progress = FALSE) 
pairwise_allometry_genus

#Distances between slope vectors (end-points) - absolute difference between slopes of groups
#if significant means int model better than comb
pairwise_allometry_genus_dist <- summary(pairwise_allometry_genus, confidence = 0.95, test.type = "dist") 
pairwise_allometry_genus_dist

#Correlation between slope vectors (and angles) - similarity of vector orientation or angle,
#if significant means the vectors of the groups are oriented in different ways 
pairwise_allometry_genus_VC <- summary(pairwise_allometry_genus, confidence = 0.95, test.type = "VC",
                                       angle.type = "deg") 
pairwise_allometry_genus_VC

#Absolute difference between slope vector lengths - difference in rate of change per covariate unit (size),
#if significant means there is a significant rate of change difference in shape between groups during growth
pairwise_allometry_genus_DL <-summary(pairwise_allometry_genus, confidence = 0.95, test.type = "DL") 
pairwise_allometry_genus_DL 

#Compare the dispersion around genus slopes - fit of the data to the regression
#if significant difference might be problem as it means the groups are not evenly sampled or one of them contains relevant outliers
pairwise_allometry_genus_var <-summary(pairwise_allometry_genus, confidence = 0.95, test.type = "var")
pairwise_allometry_genus_var

#Save results to file
sink("Output/pairwise_allometry_genus.txt")
print("ANOVA models")
print(anova_allometry_models_genus)

print("1-Pairwise absolute distances slopes")
pairwise_allometry_genus_dist

print("2-Distance between angles (slope directions)")
pairwise_allometry_genus_VC

print("3-Difference in slope vector length (difference in rate of change of shape per unit of size)")
pairwise_allometry_genus_DL

print("4-Difference in dispersion around mean slope")
pairwise_allometry_genus_var
sink()

###Heatmaps plots for significant differences in pairwise ----

#Create palette for heatmap plot
mypalette_seq_reds <- brewer.pal(9,"Reds")
image(1:9,1, as.matrix(1:9), col = mypalette_seq_reds,xlab="Reds (sequential)",
      ylab = "", yaxt = "n")

#Save p-values as object
pairwise_allometry_dist <- pairwise_allometry_genus_dist[["pairwise.tables"]][["D"]]
pairwise_allometry_dist_p <- pairwise_allometry_genus_dist[["pairwise.tables"]][["P"]]
pairwise_allometry_angle <- pairwise_allometry_genus_VC[["pairwise.tables"]][["angle"]]
pairwise_allometry_angle_p <- pairwise_allometry_genus_VC[["pairwise.tables"]][["P"]]
pairwise_allometry_length <- pairwise_allometry_genus_DL[["pairwise.tables"]][["D"]]
pairwise_allometry_length_p <- pairwise_allometry_genus_DL[["pairwise.tables"]][["P"]]

#Make list to change tables faster
pairwise_allometry_list <- list(pairwise_allometry_dist, pairwise_allometry_dist_p, pairwise_allometry_angle, pairwise_allometry_angle_p, 
                                pairwise_allometry_length, pairwise_allometry_length_p)

#Make list of shorter genera
genera_list <- levels(genera)
genera_list_short <- str_sub(genera_list, 1, 4)
genera_list_short <- c("Bala", "Bals" ,"Ball" ,"Cape", "Esch" ,"Kogi" ,
                       "Lage", "Mega","Phoc") #change ones that are the same by printing list

#Set correct row and col names for both
#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  rownames(pairwise_allometry_list[[l]]) <- genera_list_short 
  colnames(pairwise_allometry_list[[l]]) <- genera_list_short 
}

#Save only lower triangle for each
pairwise_allometry_lower_tri_list <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  pairwise_allometry_lower_tri_list[[l]] <- get_upper_tri(pairwise_allometry_list[[l]])
}

#Melt to make table in the format needed for heatmap
pairwise_allometry_melt <- list()

#Loop
for (l in 1:6){   #number of variable is fixed, given by parwise results
  pairwise_allometry_melt[[l]] <- melt(pairwise_allometry_lower_tri_list[[l]], na.rm = TRUE)
}

#Create single data frames 
pairwise_allometry_dist_melt <- data.frame(pairwise_allometry_melt[[1]], p = pairwise_allometry_melt[[2]][[3]])
pairwise_allometry_angle_melt <- data.frame(pairwise_allometry_melt[[3]], p = pairwise_allometry_melt[[4]][[3]])
pairwise_allometry_length_melt <- data.frame(pairwise_allometry_melt[[5]], p = pairwise_allometry_melt[[6]][[3]])

#Create columns for group values for Var1 and Var2 - useful for plotting by group heatmaps
pairwise_allometry_dist_melt$group1 <- if_else(pairwise_allometry_dist_melt$Var1 %in% c("Bala", "Bals", "Ball" ,"Cape","Esch", "Mega"), 
                                               "myst", "odon")
pairwise_allometry_dist_melt$group2 <- if_else(pairwise_allometry_dist_melt$Var2 %in% c("Bala", "Bals", "Ball" ,"Cape","Esch", "Mega"), 
                                               "myst", "odon")

pairwise_allometry_angle_melt$group1 <- if_else(pairwise_allometry_angle_melt$Var1 %in% c("Bala", "Bals", "Ball" ,"Cape","Esch", "Mega"), 
                                                "myst", "odon")
pairwise_allometry_angle_melt$group2 <- if_else(pairwise_allometry_angle_melt$Var2 %in% c("Bala", "Bals", "Ball" ,"Cape","Esch", "Mega"), 
                                                "myst", "odon")

pairwise_allometry_length_melt$group1 <- if_else(pairwise_allometry_length_melt$Var1 %in% c("Bala", "Bals", "Ball" ,"Cape","Esch", "Mega"), 
                                                 "myst", "odon")
pairwise_allometry_length_melt$group2 <- if_else(pairwise_allometry_length_melt$Var2 %in% c("Bala", "Bals", "Ball" ,"Cape","Esch", "Mega"), 
                                                 "myst", "odon")

#Create columns where only significant values are shown
pairwise_allometry_dist_melt <- pairwise_allometry_dist_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                        p_if_sig = ifelse(sig_p, p, NA),
                                                                        value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 2)))

pairwise_allometry_angle_melt <- pairwise_allometry_angle_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                          p_if_sig = ifelse(sig_p, p, NA),
                                                                          value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 2)))

pairwise_allometry_length_melt <- pairwise_allometry_length_melt %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                            p_if_sig = ifelse(sig_p, p, NA),
                                                                            value_if_sig = ifelse(sig_p, value, NA)) %>%
  mutate_at(vars(starts_with("value")), list(~ round(., 2)))

#Heatmaps will give error if no variables significant!!
#Check NAs first - if TRUE do not plot
all(is.na(pairwise_allometry_dist_melt$p_if_sig))

all(is.na(pairwise_allometry_angle_melt$p_if_sig))

all(is.na(pairwise_allometry_length_melt$p_if_sig))

#Nice heatmap plot for each variable
pairwise_allometry_dist_heatmap_ggplot <- ggplot(data = pairwise_allometry_dist_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_seq_reds[9], high = mypalette_seq_reds[2], mid = mypalette_seq_reds[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_reds[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope distance")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(angle = 45, size = 12, vjust = 0.9, hjust = 0.9, margin = NULL),
        axis.text.y =  element_text(size = 12, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.5, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 13), legend.text = element_text(size = 11))+
  guides(fill = guide_colorbar(barwidth = 8, barheight = 1.2,
                               title.position = "top", title.hjust = 0.5))
pairwise_allometry_dist_heatmap_ggplot

#Nice heatmap plot for each variable
pairwise_allometry_angle_heatmap_ggplot <- ggplot(data = pairwise_allometry_angle_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_seq_reds[9], high = mypalette_seq_reds[2], mid = mypalette_seq_reds[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_reds[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope angle difference")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(angle = 45, size = 12, vjust = 0.9, hjust = 0.9, margin = NULL),
        axis.text.y =  element_text(size = 12, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 13), legend.text = element_text(size = 11))+
  guides(fill = "none")
pairwise_allometry_angle_heatmap_ggplot

#Nice heatmap plot for each variable
pairwise_allometry_length_heatmap_ggplot <- ggplot(data = pairwise_allometry_length_melt, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  geom_text(aes(Var2, Var1, label = value_if_sig), color = "white", size = 4) +
  scale_fill_gradient2(low = mypalette_seq_reds[9], high = mypalette_seq_reds[2], mid = mypalette_seq_reds[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0.001, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_reds[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  ggtitle ("Slope length difference")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(angle = 45, size = 12, vjust = 0.9, hjust = 0.9, margin = NULL),
        axis.text.y =  element_text(size = 12, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.justification = c(1, 0), legend.position = c(0.4, 0.7),  legend.direction = "horizontal",
        legend.title = element_text(size = 13), legend.text = element_text(size = 11))+
  guides(fill = "none")
pairwise_allometry_length_heatmap_ggplot

ggarrange(pairwise_allometry_dist_heatmap_ggplot, pairwise_allometry_angle_heatmap_ggplot,
          pairwise_allometry_length_heatmap_ggplot,
          ncol = 3, nrow = 1, common.legend = F)

####Plots by group####
#Make data frame for each stage
pairwise_allometry_dist_melt_mysticeti <- pairwise_allometry_dist_melt %>% filter(str_detect(group1, "myst")) %>% 
  filter(str_detect(group2, "myst"))

pairwise_allometry_angle_melt_mysticeti <- pairwise_allometry_angle_melt %>% filter(str_detect(group1, "myst")) %>% 
  filter(str_detect(group2, "myst"))

pairwise_allometry_length_melt_mysticeti <- pairwise_allometry_length_melt %>% filter(str_detect(group1, "myst")) %>% 
  filter(str_detect(group2, "myst"))

#Create object for labels
genera_list_mysticeti <- c("Bala", "Bals", "Ball" ,"Cape","Esch", "Mega")

#Nice heatmap plot
pairwise_allometry_dist_heatmap_ggplot_mysticeti  <- ggplot(data = pairwise_allometry_dist_melt_mysticeti, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  scale_fill_gradient2(low = mypalette_seq_reds[9], high = mypalette_seq_reds[2], mid = mypalette_seq_reds[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_reds[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  scale_x_discrete(labels = genera_list_mysticeti)+
  scale_y_discrete(labels = genera_list_mysticeti)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  ggtitle ("Distance - Mysticeti")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(angle = 45, size = 10,hjust = 0.9),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.position = c(0.18,0.85),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.1,
                               title.position = "top", title.hjust = 0.5))
pairwise_allometry_dist_heatmap_ggplot_mysticeti 

pairwise_allometry_angle_heatmap_ggplot_mysticeti  <- ggplot(data = pairwise_allometry_angle_melt_mysticeti, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  scale_fill_gradient2(low = mypalette_seq_reds[9], high = mypalette_seq_reds[2], mid = mypalette_seq_reds[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_reds[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+   
  scale_x_discrete(labels = genera_list_mysticeti)+   
  scale_y_discrete(labels = genera_list_mysticeti)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  ggtitle ("Angle - Mysticeti")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(angle = 45, size = 10,hjust = 0.9),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.position = c(0.18,0.85),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.1,
                               title.position = "top", title.hjust = 0.5))
pairwise_allometry_angle_heatmap_ggplot_mysticeti 

pairwise_allometry_length_heatmap_ggplot_mysticeti  <- ggplot(data = pairwise_allometry_length_melt_mysticeti, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  scale_fill_gradient2(low = mypalette_seq_reds[9], high = mypalette_seq_reds[2], mid = mypalette_seq_reds[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_reds[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+   
  scale_x_discrete(labels = genera_list_mysticeti)+   
  scale_y_discrete(labels = genera_list_mysticeti)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  ggtitle ("Length - Mysticeti")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(angle = 45, size = 10,hjust = 0.9),
        axis.text.y =  element_text(size = 10, vjust = -0.2, margin = NULL), panel.grid.major = element_blank(),
        legend.position = c(0.18,0.85),  legend.direction = "horizontal",
        legend.title = element_text(size = 11), legend.text = element_text(size = 10))+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1.1,
                               title.position = "top", title.hjust = 0.5))
pairwise_allometry_length_heatmap_ggplot_mysticeti 

ggarrange(pairwise_allometry_dist_heatmap_ggplot_mysticeti, pairwise_allometry_angle_heatmap_ggplot_mysticeti,
          pairwise_allometry_length_heatmap_ggplot_mysticeti, 
          ncol = 3, nrow = 1, common.legend = T, heights = c(1,1.2))

###Plot regression by genus ----
#Regression score of shape vs logCS and comb or int (best model)- regression method with "RegScore" plotting
allometry_genus_plot_regscore <- plot(allometry_genus_int, type = "regression",predictor = gdf$size, reg.type = "RegScore",
                                      main = "Shape vs logCS * genus",xlab = "logCS", pch = 21, col = "chartreuse4", bg = "chartreuse4", cex = 1.2, font.main = 2)   #improve graphics
text(x = gdf$size, y = allometry_genus_plot_regscore$RegScore, labels = Ids,
     pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels

##Make better allometry plot with ggplot
#Create data frame object that ggplot can read - use data from plot object you want to improve
allometry_genus_plot <- data.frame(logCS = allometry_genus_plot_regscore[["plot.args"]][["x"]], 
                                   RegScores = allometry_genus_plot_regscore[["plot.args"]][["y"]])
#Convert data frame to tibble
allometry_genus_plot <- as_tibble(allometry_genus_plot)
#Add labels and other attributes to tibble as columns
allometry_genus_plot <- allometry_genus_plot %>% 
  mutate(specimens = gdf$Id, group = gdf$group, category = gdf$category, family = gdf$family, stage = gdf$stage,
         genus = gdf$genus, size = gdf$size, TL = gdf$TL_100, BZW = gdf$BZW_100)
glimpse(allometry_genus_plot)

#Plot allometry regression by genus
allometry_genus_ggplot <- ggplot(allometry_genus_plot, aes(x = logCS, y = RegScores))+
  geom_point(size = 4, aes(fill = genus, colour = genus, shape = group), alpha = 0.5)+  
  geom_smooth(method = "lm", aes(x = logCS, y = RegScores, linetype = group, colour = genus), inherit.aes = F,        
              size = 1, alpha = 0.4, se = F, show.legend = T)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_colour_manual(name = "Genus", labels = levels(genera), 
                      values = mypalette_taxa, aesthetics = c("colour","fill"))+         
  scale_shape_manual(name = "Group", labels =  levels(groups), 
                     values = shapes)+
  scale_linetype_manual(name = "Group", labels =  levels(groups), 
                        values = c(1,2))+
  theme_classic(base_size = 12)+
  ylab("Regression Score - p = 0.001**")+
  theme(legend.key = element_blank(), legend.background = element_blank(),
        legend.title = element_text(size = 11, face = "bold"), 
        legend.position = c(0.7,0.1),  legend.direction = "vertical", legend.justification = c(0,0))+
  guides(colour = guide_legend(label.theme = element_text(size =10, angle = 0, face = "italic"),
                               override.aes = list(size = 3.5, linetype = 0)), 
         linetype = guide_legend(override.aes = list(colour = "black")), shape = guide_legend(override.aes = list(fill = "gray80")))
allometry_genus_ggplot

#Plot groups separately
allometry_genus_plot_mysticeti <- allometry_genus_plot %>% filter(group == "mysticeti")
allometry_genus_plot_odontoceti <- allometry_genus_plot %>% filter(group == "odontoceti")

#Plot allometry regression by genus mysticeti
allometry_genus_ggplot_mysticeti <- ggplot(allometry_genus_plot_mysticeti, aes(x = logCS, y = RegScores))+
  geom_point(size = 3, aes(fill = genus, colour = genus, shape = family), alpha = 0.3)+  
  geom_smooth(method = "lm", aes(x = logCS, y = RegScores, colour = genus), linetype = 1, inherit.aes = F,        
              size = 0.8, alpha = 0.4, se = F, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_colour_manual(name = "Genus", labels = levels(as.factor(allometry_genus_plot_mysticeti$genus)), 
                      values = mypalette_taxa[-c(6,7,9)], aesthetics = c("colour","fill"),
                      guide = guide_legend(label.theme = element_text(size =11, angle = 0, face = "italic")))+ 
  scale_shape_manual(name = "Family", labels = levels(families)[c(1:2,4)],
                     values = shapes_fam[c(1:2,4)], 
                     guide = guide_legend(override.aes = list(size = 4)))+
  theme_classic(base_size = 12)+
  ylab("Regression Score - p = 0.001**")+
  theme(legend.key = element_blank(), legend.background = element_blank(),
        legend.title = element_text(size = 11, face = "bold"), 
        legend.position = c(0.6,0.1),  legend.direction = "vertical", legend.box = "horizontal", 
        legend.justification = c(0,0))
allometry_genus_ggplot_mysticeti

#Add phylopic
allometry_genus_ggplot_mysticeti <- 
  allometry_genus_ggplot_mysticeti +
  add_phylopic(Bals, alpha = 1, x = 3.5, y = -0.02, ysize = 0.04, color = mypalette_taxa[3])+
  add_phylopic(Ball, alpha = 1, x = 4.4, y = 0.15, ysize = 0.04, color = mypalette_taxa[2])+
  add_phylopic(Bala, alpha = 1, x = 3.8, y = 0.2, ysize = 0.05, color = mypalette_taxa[1])+
  add_phylopic(Cap, alpha = 1, x = 3.8, y = 0.14, ysize = 0.05, color = mypalette_taxa[4])+
  add_phylopic(Esch, alpha = 1, x = 4.2, y = 0.2, ysize = 0.04, color = mypalette_taxa[5])+
  add_phylopic(Meg, alpha = 1, x = 4.1, y = 0.09, ysize = 0.05, color = mypalette_taxa[8])
allometry_genus_ggplot_mysticeti

#Plot allometry regression by genus odontoceti
allometry_genus_ggplot_odontoceti <- ggplot(allometry_genus_plot_odontoceti, aes(x = logCS, y = RegScores))+
  geom_point(size = 3, aes(fill = genus, colour = genus, shape = family), alpha = 0.8)+  
  geom_smooth(method = "lm", aes(x = logCS, y = RegScores, colour = genus), linetype = 2, inherit.aes = F,        
              size = 0.8, alpha = 0.4, se = F, show.legend = F)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_colour_manual(name = "Genus", labels = levels(as.factor(allometry_genus_plot_odontoceti$genus)), 
                      values = mypalette_taxa[c(6,7,9)], aesthetics = c("colour","fill"),
                      guide = guide_legend(label.theme = element_text(size = 11, angle = 0, face = "italic")))+ 
  scale_shape_manual(name = "Family", labels = levels(families)[-c(1:2,4)],
                     values = shapes_fam[-c(1:2,4)])+
  theme_classic(base_size = 12)+
  ylab("Regression Score - p = 0.001**")+
  theme(legend.key = element_blank(), legend.background = element_blank(),
        legend.title = element_text(size = 11, face = "bold"), legend.box = "horizontal",
        legend.position = c(0.6,0.1),  legend.direction = "vertical", legend.justification = c(0,0))
allometry_genus_ggplot_odontoceti

#Add phylopic
allometry_genus_ggplot_odontoceti <- 
  allometry_genus_ggplot_odontoceti +
  add_phylopic(Kog, alpha = 1, x = 3.2, y = -0.12, ysize = 0.03, color = mypalette_taxa[6])+
  add_phylopic(Lage, alpha = 1, x = 3.2, y = 0.05, ysize = 0.08, color = mypalette_taxa[7])+
  add_phylopic(Phoc, alpha = 1, x = 2.7, y = -0.07, ysize = 0.04, color = mypalette_taxa[9])
allometry_genus_ggplot_odontoceti

ggarrange(allometry_genus_ggplot_mysticeti,allometry_genus_ggplot_odontoceti, 
          ncol = 2, nrow = 1, widths = c(1.2,1))

#SLOPES AND INTERCEPTS FOR ANCESTRAL STATE RECONTRUCTION OF ALLOMETRY ----

#Get coefficents for each genus per stage to use later in anc state
#Regression score of shape vs logCS - regression method with "RegScore" plotting
allometry_genus_plot_regscore <- plot(allometry_genus_int, type = "regression",predictor = gdf$size, reg.type = "RegScore", 
                                      main = "Shape vs logCS  by genus",xlab = "logCS", pch = 21, cex = 1.2, font.main = 2,
                                      col = c(mypalette_taxa, mypalette_taxa), bg = c(mypalette_taxa, mypalette_taxa))
text(x = gdf$size, y = allometry_genus_plot_regscore$RegScore, labels = gdf$Id,
     pos = 3, offset = 0.5, cex = 0.75)    #improve appearance of labels  
#Create object to use for linear model
allometry_genus_regscores <- allometry_genus_plot_regscore[["RegScore"]] 

#Linear model for line by genus
allometry_genus_regline_df <- data.frame(RegScores = allometry_genus_regscores, logCS = gdf$size, genus = gdf$genus)

allometry_genus_regline <- lm(RegScores ~ logCS * genus, data = allometry_genus_regline_df)

#Get coeffs - the first 2 are reference intercept and slopes, other values are differences!
allometry_genus_regline_coeffs <- as.matrix(allometry_genus_regline$coefficients)

#Save intercepts and slopes separately
allometry_genus_regline_intercepts <- as.matrix(allometry_genus_regline_coeffs[c(1, 3:(length(genera_list)+1)),])
allometry_genus_regline_slopes <- as.matrix(allometry_genus_regline_coeffs[c(2, length(genera_list)+2:(length(genera_list))),])

#Calculate real intercepts and slopes
allometry_genus_regline_intercepts_ok <- as.matrix(c(allometry_genus_regline_intercepts[1,], allometry_genus_regline_intercepts[1,]+
                                                       allometry_genus_regline_intercepts[2:length(allometry_genus_regline_intercepts),]))

allometry_genus_regline_slopes_ok <- as.matrix(c(allometry_genus_regline_slopes[1,], allometry_genus_regline_slopes[1,]+
                                                   allometry_genus_regline_slopes[2:length(allometry_genus_regline_slopes),]))

#Save as data frame with grouping variables
allometry_genus_coeffs <- data.frame(Slope = allometry_genus_regline_slopes_ok, Intercept = allometry_genus_regline_intercepts_ok, 
                                     row.names = levels(as.factor(gdf$genus)))
#Check for NA and other issues
allometry_genus_coeffs

allometry_genus_coeffs <- allometry_genus_coeffs %>% mutate(genus = rownames(allometry_genus_coeffs))

allometry_genus_coeffs <- merge(allometry_genus_coeffs, families_genera_groups_list, by = "genus")
allometry_genus_coeffs 

###### 
#Next - ch. 6 - Ancestral state reconstruction allometry
