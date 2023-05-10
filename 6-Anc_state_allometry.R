
#===========================================================#
#                                                           #
#       CURVES AND POINTS ANALYSES - MYSTICETI ONLY         #
#                                                           #
#===========================================================#

#CH. 6 - Ancestral state reconstruction allometry slope

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
library(rray)
library(abind)
library(mcp)
library(emmeans)
library(grateful)
library(ape)
library(viridis)
set.seed(17)

#require(devtools)
#install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")
#remotes::install_github("r-lib/rray")
#remotes::install_github("Pakillo/grateful")

#cite_packages(style = "pnas", out.format = "docx")

#apropos("x") lists objects with matching part of name

#ALLOMETRY WITH ANC STATE REC (modified from Morris et al., 2019) ----

#Import trees in Nexus format - branch lengths needed!!
tree_genera <- "Data/tree_myst.txt"   #tree with all genera
#Taxa with both stages very few, not worth it
#tree_genera_stages <- "Data/tree_stages.txt" #tree only with genera that have both prenatal and postnatal stages

##Read the trees for analysis
tree_all1 <- read.nexus(tree_genera) #tree with all genera
plot(tree_all1)

#Make sure tip labels match taxa names in data frame
wrong_tips_all <- sort(tree_all1$tip.label)
tree_all <- tree_all1

tree_all$tip.label <- genera_list[match(tree_all1$tip.label, wrong_tips_all)]
plot(tree_all)

#Check names of each tree in object
summary(tree_all)

##ASR tree plot ----
#Check names correspond
setdiff(tree_all$tip.label, allometry_genus_coeffs$genus)

#Reconstruct ancestral slopes and intercepts
slope_genus <- setNames(allometry_genus_coeffs$Slope, allometry_genus_coeffs$genus)
int_genus <- setNames(allometry_genus_coeffs$Intercept, allometry_genus_coeffs$genus)

slope_fit_genus <- fastAnc(tree_all, slope_genus, vars=TRUE,CI=TRUE)
print(slope_fit_genus, printlen=10)

int_fit_genus <- fastAnc(tree_all, int_genus, vars=TRUE,CI=TRUE)
print(int_fit_genus, printlen=10)

#Plot on tree
slope_mapping_genus <- contMap(tree_all, slope_genus, plot=FALSE)

plot(slope_mapping_genus,legend=0.7*max(nodeHeights(tree_all)),
     sig=2,fsize=c(0.7,0.9))

int_mapping_genus <- contMap(tree_all, int_genus, plot=FALSE)

plot(int_mapping_genus,legend=0.7*max(nodeHeights(tree_all)),
     sig=2,fsize=c(0.7,0.9))

#Costum colors tree mapping
my_colors_genus <- viridis(10, option = "turbo", direction = 1)
slope_mapping_genus <- setMap(slope_mapping_genus, my_colors_genus)

plot(slope_mapping_genus,legend=FALSE,
     ylim=c(1-0.09*(Ntip(slope_mapping_genus$tree)-1),Ntip(slope_mapping_genus$tree)),
     mar=c(5.1,0.4,0.4,0.4))
add.color.bar(40,slope_mapping_genus$cols,title="Slope",
              lims=slope_mapping_genus$lims,digits=3,prompt=FALSE,x=0,
              y=1-0.08*(Ntip(slope_mapping_genus$tree)-1),lwd=4,fsize=1,subtitle="")

int_mapping_genus <- setMap(int_mapping_genus, my_colors_genus)

plot(int_mapping_genus,legend=FALSE,
     ylim=c(1-0.09*(Ntip(int_mapping_genus$tree)-1),Ntip(int_mapping_genus$tree)),
     mar=c(5.1,0.4,0.4,0.4))
add.color.bar(40,int_mapping_genus$cols,title="Intercept",
              lims=int_mapping_genus$lims,digits=3,prompt=FALSE,x=0,
              y=1-0.08*(Ntip(int_mapping_genus$tree)-1),lwd=4,fsize=1,subtitle="")

#Save results to file
sink("Output/anc slope interecept whole skull.txt")
print("Taxa")
print(slope_genus)
print(int_genus)
print("Anc nodes")
print(slope_fit_genus$ace)
print(int_fit_genus$ace)
sink() 

##ASR regression plot ----

#Calculate ancestral states with phytools - slope and intercept for each component
slope_anc.ML <- fastAnc(tree_all, allometry_genus_coeffs[,"Slope"])
int_anc.ML <- fastAnc(tree_all, allometry_genus_coeffs[,"Intercept"])

#Prepare data frame for plot
#Plot tree with node labels to check what node is where
plot(tree_all, show.node.label = T)
nodelabels() #save tree figure
#Save node numbers
nodes_all <- as.numeric(names(slope_anc.ML))
#Create col names for ancestral data frame
anc_col_names <- c("Slope","Intercept")

#Create data frames with ancestral state values
anc_values_traj <- matrix(nrow = length(nodes_all), ncol = 2, dimnames = list(nodes_all,anc_col_names))
anc_values_traj[,"Slope"] <- slope_anc.ML
anc_values_traj[,"Intercept"] <- int_anc.ML
anc_values_traj
#Add nodes_all
anc_values_traj <- as.data.frame(anc_values_traj)
anc_values_traj$genus <- as.character(nodes_all) 
#Order based on node number
anc_values_traj <- anc_values_traj[order(anc_values_traj$genus),]
anc_values_traj

##Plot allometric trajectories with abline for each group and ancestral node
#Check groups and variables
glimpse(allometry_genus_plot)

#Make tibble and add columns for plotting
anc_values_traj <- anc_values_traj %>% as_tibble() %>% mutate(group = "anc_node")

#Create palette including anc nodes_all
mypalette_greys <- as.matrix(ggthemes_data[["tableau"]][["color-palettes"]][["ordered-sequential"]][["Gray Warm"]][["value"]])
mypalette_greys[7] <- "#b3aba6" #change cause there is a wrong value in palette
image(1:20,1, as.matrix(1:20), col = mypalette_greys, xlab="Greys (sequential)",
      ylab = "", yaxt = "n")
mypalette_greys2 <- mypalette_greys[c(4:20)]

mypalette_browns <- as.matrix(ggthemes_data[["tableau"]][["color-palettes"]][["ordered-sequential"]][["Brown"]][["value"]])
image(1:20,1, as.matrix(1:20), col = mypalette_browns, xlab="Browns (sequential)",
      ylab = "", yaxt = "n")
mypalette_browns2 <- mypalette_browns[c(15:20)]

#Palette for nodes_all
mypalette_nodes_all <- colorRampPalette(c(rev(mypalette_greys2[c(TRUE, FALSE)]),rev(mypalette_browns2)))
mypalette_nodes_all(length(nodes_all))
plot(rep(1,length(nodes_all)),col=mypalette_nodes_all(length(nodes_all)),pch=19,cex=3)
mypalette_nodes_all1 <- mypalette_nodes_all(length(nodes_all))
#reorder colors to make consecutive nodes easier to tell apart
mypalette_nodes_all2 <- c("#59504E", "#A63D32","#BD6036", "#A9A09D", "#736967", "#8D8481",  "#B34D34", "#C5BDB9")

mypalette_taxa_nodes_all <- c(mypalette_taxa, mypalette_nodes_all2)
plot(rep(1,length(mypalette_taxa_nodes_all)),col=mypalette_taxa_nodes_all ,pch=19,cex=3, 
     main = "Taxa and nodes colors", ylab = "", xlab = "" ,cex.main = 1.5)
title(xlab = "1-Balaena 2-Balaenoptera_l 3-Balaenoptera_s 4-Caperea 5-Eschrichtius 
   6-Kogia 7-Lagenorhynchus 8-Megaptera 9-Phocoena, 10-17 anc nodes", cex.lab = 1, font.lab = 3, line = -2)
text(x = seq_along(1:length(mypalette_taxa_nodes_all)), y = 1.05, 
     labels = seq_along(1:length(mypalette_taxa_nodes_all)))

mypalette_nodes_all_only <- mypalette_taxa_nodes_all[length(mypalette_taxa)+1:length(mypalette_taxa_nodes_all)]
mypalette_nodes_all_only <- as.vector(na.omit(mypalette_nodes_all_only))
plot(rep(1,length(mypalette_nodes_all_only)),col=mypalette_nodes_all_only ,pch=19,cex=3, 
     main = "Nodes colors", ylab = "", xlab = "" ,cex.main = 1.5)
title(xlab = "1 - anc node all Cetacea 10
2-6 - anc nodes Mysticeti 11-15
7-8 - anc nodes Odontoceti 16-17", cex.lab = 1, font.lab = 3, line = -2)
text(x = seq_along(1:length(mypalette_nodes_all_only)), y = 1.05, 
     labels = seq_along(1:length(mypalette_nodes_all_only)))

#Plot both stages together - abline to check it works and line are where expected
allometry_anc_nodes_all_ggplot <- allometry_genus_ggplot +
  #line on plot
  geom_abline(data = anc_values_traj, 
              aes(intercept = Intercept, slope = Slope, group = genus), linetype = 2, size  = 1)
allometry_anc_nodes_all_ggplot

##TEST DIFFERENCES BETWEEN ANC STATES AND GROUPS IN SLOPE/INTERCEPT - PAIRWISE ----

###Create estimated values for anc and genera ----

#Make x values (size) based on modern
anc_X <- expand.grid(logCS = seq(min(gdf$size), 
                                 max(gdf$size), #use min and max of x values as limits  
                                 length.out = length(genera_list)*length(nodes_all)))  #length of sequence must be multiple of taxa for later

#Calculate Y values (reg scores) for ancestral nodes_all using anc.ML intercept and slope values
#List slopes and intercepts to make loop work
#Empty objects
slopes_anc <- as.list(slope_anc.ML)
intercepts_anc <- as.list(int_anc.ML)

#Empty object
anc_Y_list <- list()

#Loop
for (a in 1:length(nodes_all)){
  anc_Y_list[[a]] <- intercepts_anc[[a]] + 
    (anc_X * slopes_anc[[a]])
}

#Make vector for Y values
anc_Y <- as.vector(unlist(anc_Y_list))

#Create groups and order variables for ancestral states data
anc_nodes_all <- cbind(rep(nodes_all, each =  dim(anc_X)[[1]]))
anc_groups <- rep("anc_node", times = dim(anc_X)[[1]])

#Create data frame with ReScores, logCS estimated and groups and orders for anc to match allometry tibble
allometry_nodes_all <- data.frame(logCS = anc_X, RegScores = anc_Y, 
                                  genus = as.character(anc_nodes_all), group = anc_groups)
glimpse(allometry_nodes_all)

#Combine with allometry tibble
colnames(allometry_genus_plot) #check columns
colnames(allometry_nodes_all)
#Get column names to delete
col_diffs <- setdiff(colnames(allometry_genus_plot),colnames(allometry_nodes_all))

#Create new data frame with taxa and anc data
allometry_genus_all <- select(allometry_genus_plot, -all_of(col_diffs))

allometry_anc_all <- bind_rows(allometry_genus_all, allometry_nodes_all)

##Order values by genus, useful for plot legend
#Make factor for variable
allometry_anc_all$genus <- factor(allometry_anc_all$genus, 
                                  levels = c(genera_list, as.character(nodes_all))) #copy from string printed with the code above
#Order
allometry_anc_all <- allometry_anc_all[order(allometry_anc_all$genus),]
#Check
glimpse(allometry_anc_all)

###Pairwise comparison of regression model between genera and ancestral nodes ----
#Create models, with different slopes and int or just int
allometry_anc_all_null <- lm.rrpp(RegScores ~ logCS,   #null model
                                  data = allometry_anc_all, print.progress = FALSE, iter = 999) 
allometry_anc_all_comb <- lm.rrpp(RegScores ~ logCS + genus,
                                  data = allometry_anc_all, print.progress = FALSE, iter = 999) 
allometry_anc_all_int <- lm.rrpp(RegScores ~ logCS * genus,
                                 data = allometry_anc_all, print.progress = FALSE, iter = 999) 

#Check results
summary(allometry_anc_all_null)
summary(allometry_anc_all_comb)
summary(allometry_anc_all_int)

#Anova for difference between models
anova_allometry_anc_models_all <- anova(allometry_anc_all_null, allometry_anc_all_comb, allometry_anc_all_int)
anova_allometry_anc_models_all

#Check slopes different for best model
anova(allometry_anc_all_int)

#If slopes different, calculate p values for pairs using emmeans package
#First recreate model with lm() 
allometry_anc_all_int1 <- lm(RegScores ~ logCS * genus,
                             data = allometry_anc_all) 
#Check anova still ok
anova(allometry_anc_all_int1)

#Get pairwise comparisons of slopes
allometry_anc_all_emms <- emmeans(allometry_anc_all_int1, "genus")

#to make graph, confusing for lots of groups - pwpp(allometry_anc_all_emms)

#Visualize table with p values for differences between groups - for heatmaps
allometry_anc_all_ems_table <- pwpm(allometry_anc_all_emms, diffs = F, means = F) #eliminate differences and means in the diagonal, only keep p values
#Transform in table
allometry_anc_all_ems_table <- char(allometry_anc_all_ems_table)

#Save results to file
sink("Output/pairwise_allometry_anc_nodes_all_genera.txt")
print("ANOVA models")
print(anova_allometry_anc_models_all)

print("summary best model - lmrpp")
anova(allometry_anc_all_int)

print("summary model used for comparisons - lm")
summary(allometry_anc_all_int1)
anova(allometry_anc_all_int1)

print("Pairwise comparison using emmeans")
summary(allometry_anc_all_emms)

print("Full results table emmeans pairwise comparions")
pwpm(allometry_anc_all_emms)
sink()

###Clean up environment before proceeding
save(allometry_anc_all_int1, file = "allometry_anc_all_int_model.RData")

rm("allometry_anc_all_int1")

##Heatmaps plots for significant differences in slopes ----

#Create palette for heatmap plot
mypalette_seq_blues <- brewer.pal(9,"Blues")
image(1:9,1, as.matrix(1:9), col = mypalette_seq_blues,xlab="Blues (sequential)",
      ylab = "", yaxt = "n")

#Make already melted table with p-values
allometry_anc_all_ems_pvals <- as.table(allometry_anc_all_ems_table)
allometry_anc_all_ems_pvals <- as.data.frame(allometry_anc_all_ems_pvals)

#Delete extra columns
allometry_anc_all_ems_pvals <- allometry_anc_all_ems_pvals[,-c(3:4)]
#Delete empty rows for p vals
allometry_anc_all_ems_pvals <- allometry_anc_all_ems_pvals[-which(allometry_anc_all_ems_pvals$Freq.Freq == ""), ]

#Change last column name
names(allometry_anc_all_ems_pvals)[3] <- "p"

#Delete < so that it can be converted to numeric
allometry_anc_all_ems_pvals$p <- str_replace(allometry_anc_all_ems_pvals$p, "<", "0")
allometry_anc_all_ems_pvals$p <- as.numeric(allometry_anc_all_ems_pvals$p)

#Make vars names shorter
#Loop replacements genera
for (t in 1:length(genera_list)){
  allometry_anc_all_ems_pvals$Var1 <- str_replace_all(allometry_anc_all_ems_pvals$Var1, genera_list[t], genera_list_short[t])
  allometry_anc_all_ems_pvals$Var2 <- str_replace_all(allometry_anc_all_ems_pvals$Var2, genera_list[t], genera_list_short[t])
}

#Loop replacements nodes_all to match text-simpler numbering
nodes_all_tree  <- as.character(nodes_all)

#Make factor for variable
allometry_anc_all_ems_pvals$Var1 <- factor(allometry_anc_all_ems_pvals$Var1, 
                                           levels = c(genera_list_short, as.character(nodes_all))) #copy from string printed with the code above
#Order
allometry_anc_all_ems_pvals <- allometry_anc_all_ems_pvals[order(allometry_anc_all_ems_pvals$Var1),]
#Check
glimpse(allometry_anc_all_ems_pvals)

#Subset data frame to keep only anc nodes_all for plots
allometry_anc_all_ems_pvals <- subset(allometry_anc_all_ems_pvals, Var2 %in% nodes_all)

#Create columns where only significant values are shown
allometry_anc_all_ems_pvals <-  allometry_anc_all_ems_pvals %>% mutate(sig_p = ifelse(p < .05, T, F),
                                                                       p_if_sig = ifelse(sig_p, p, NA))
allometry_anc_all_ems_pvals

#Nice heatmap plot
pairwise_allometry_anc_all_heatmap_ggplot  <- ggplot(data = allometry_anc_all_ems_pvals, aes(Var2, Var1, fill = p_if_sig))+
  geom_tile(colour = "gray80")+
  scale_fill_gradient2(low = mypalette_seq_blues[9], high = mypalette_seq_blues[2], mid = mypalette_seq_blues[5], #negative correlations are in blue color and positive correlations in red. 
                       midpoint = 0.03, limit = c(0, 0.049), space = "Lab", #scale is from min to max p-values
                       na.value =  mypalette_seq_blues[1], name = "P-values < 0.05") + 
  theme_minimal()+ 
  coord_fixed()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  ggtitle ("Slope differences anc states and taxa")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x =  element_text(angle = 45, size = 11,hjust = 0.9),
        axis.text.y =  element_text(size = 12, margin = NULL), panel.grid.major = element_blank(),
        legend.position = c(0.25,0.85),  legend.direction = "horizontal",
        legend.title = element_text(size = 10), legend.text = element_text(size = 8))+
  guides(fill = guide_colorbar(barwidth = 6, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
pairwise_allometry_anc_all_heatmap_ggplot

###Plot allometry tips and anc nodes_all ----
#Improve plot using new dataset
allometry_anc_all_ggplot  <- ggplot(allometry_anc_all, aes(x = logCS, y = RegScores, colour = genus))+
  geom_smooth(data = allometry_anc_all, aes(x = logCS, y = RegScores, colour = genus), method = 'lm',          #confidence intervals and reg line, before points
              se = F, show.legend = T, size = 1.2)+      #put col and other graphics OUTSIDE of aes()!!!
  #points after, so they are on top
  scale_color_manual(name = "Nodes", labels = c(genera_list,
                                                nodes_all_tree), values = mypalette_taxa_nodes_all)+         
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle ("Allometry by genus with ancestral nodes")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.text = element_text(size = 10), 
        legend.title = element_text(size = 11), legend.position = c(0.7,0), legend.direction = "vertical", legend.justification = c(0.4,0),
        legend.key = element_blank(), legend.background = element_blank(), legend.box.margin = margin(0,0,0,0 ,unit = "pt"))+
  guides(colour =  guide_legend(ncol=4), linetype = guide_legend(override.aes = list(color = "black", size = 0.8)))
allometry_anc_all_ggplot 

#Make new tibbles with only coeffs for genera and anc nodes_all
allometry_anc_all_coeffs <- bind_rows(allometry_genus_coeffs,anc_values_traj)

#Make factor for variable
allometry_anc_all_coeffs$genus <- factor(allometry_anc_all_coeffs$genus, 
                                         levels = c(genera_list, as.character(nodes_all))) #copy from string printed with the code above
#Order
allometry_anc_all_coeffs <- allometry_anc_all_coeffs[order(allometry_anc_all_coeffs$genus),]
#Check 
allometry_anc_all_coeffs$genus

#Plot with abline to get full lenght - easier comparison
allometry_anc_all_ggplot  <- ggplot(allometry_anc_all, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = allometry_anc_all_coeffs, 
              aes(intercept = Intercept, slope = Slope,  colour = genus, linetype = group, alpha = group), size  = 1.2)+
  #points after, so they are on top
  scale_linetype_manual(name = "Group", labels = c("ancestral node","Mysticeti", "Odontoceti"), values = c(4, 1:2))+
  scale_color_manual(name = "Nodes", labels = c(genera_list,
                                                nodes_all_tree), values = mypalette_taxa_nodes_all)+ 
  scale_alpha_manual(values = c(1, 0.5, 0.5))+
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  ggtitle ("Allometry by genus with ancestral nodes")+ 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.text = element_text(size = 10), 
        legend.title = element_text(size = 11), legend.position = c(0.7,0), legend.direction = "vertical", 
        legend.justification = c(0.4,0), legend.box = "horizontal",
        legend.key = element_blank(), legend.background = element_blank(), legend.box.margin = margin(0,0,0,0 ,unit = "pt"))+
  guides(colour =  guide_legend(ncol=4), alpha = "none",
         linetype = guide_legend(override.aes = list(color = "black", size = 0.8)))
allometry_anc_all_ggplot

###Plot with selected nodes_all/taxa - clearer plots ----

#The lines are too different in length, use abline by stage
#Plot myst odont + anc node 
allometry_anc_all_coeffs_groups <- allometry_genus_coeffs %>% group_by(group) %>% 
  summarise(Slope = mean(Slope), Intercept = mean(Intercept))
allometry_anc_all_coeffs_groups <- bind_rows(allometry_anc_all_coeffs_groups, anc_values_traj[1,])
allometry_anc_all_coeffs_groups

#Plot prenatal with abline to get full length - easier comparison
allometry_anc_all_groups_ggplot  <- ggplot(allometry_anc_all, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = allometry_anc_all_coeffs_groups, 
              aes(intercept = Intercept, slope = Slope,  colour = group, alpha = group, linetype = group), size  = 1.2)+
  #points after, so they are on top
  scale_color_manual(values = c(mypalette_taxa_nodes_all[10], mypalette_groups))+ 
  scale_alpha_manual(values = c(1, 0.5, 0.5))+
  scale_linetype_manual(values = c(2, 1, 1))+
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14), legend.text = element_text(size = 10), 
        legend.title = element_text(size = 11, face = "bold"), legend.position = c(0.7,0), legend.direction = "vertical", 
        legend.justification = c(0.4,0), legend.box = "horizontal",
        legend.key = element_blank(), legend.background = element_blank(), legend.box.margin = margin(0,0,0,0 ,unit = "pt"))+
  guides(colour = "none", alpha = "none", linetype = "none")
allometry_anc_all_groups_ggplot

#Add phylopic
allometry_anc_all_groups_ggplot <- 
  allometry_anc_all_groups_ggplot +
  add_phylopic(myst, alpha = 0.8, x = 4, y = 0.2, ysize = 0.07, color = mypalette_groups[1])+
  add_phylopic(odont, alpha = 0.8, x = 2.5, y = -0.12, ysize = 0.13, color = mypalette_groups[2])+
  annotate("text", x = 3.5, y = -0.05, label = "common ancestor \n Neoceti", color = mypalette_taxa_nodes_all[10],
           family = "", fontface = 3, size = 6)
allometry_anc_all_groups_ggplot

#Plot myst + genus odont and anc node all
#Check which nodes_all/taxa apply
allometry_anc_all_coeffs_odont <- allometry_genus_coeffs %>% filter(group == "odontoceti")
allometry_anc_all_coeffs_odont <- bind_rows(allometry_anc_all_coeffs_odont, anc_values_traj[1,]) %>%
  mutate(genus = coalesce(genus,group), family = coalesce(family,group))
allometry_anc_all_coeffs_odont
allometry_anc_all_coeffs_odont_myst <- bind_rows(allometry_anc_all_coeffs_odont,allometry_anc_all_coeffs_groups[1,]) %>%
  mutate(genus = coalesce(genus,group), family = coalesce(family,group))
allometry_anc_all_coeffs_odont_myst 
allometry_anc_all_coeffs_odont_myst$genus[5] <- "Mysticeti"

#Add additional grouping for linetype plot
allometry_anc_all_coeffs_odont_myst$line <- c(rep(2, times = 3), 4,1)
allometry_anc_all_coeffs_odont_myst

allometry_anc_all_coeffs_odont_myst$genus <- as.factor(allometry_anc_all_coeffs_odont_myst$genus)
allometry_anc_all_coeffs_odont_myst <- allometry_anc_all_coeffs_odont_myst[order(allometry_anc_all_coeffs_odont_myst$genus),]
allometry_anc_all_coeffs_odont_myst 

#Plot with abline to get full length - easier comparison
allometry_anc_all_odont_myst_ggplot  <- ggplot(allometry_anc_all, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = allometry_anc_all_coeffs_odont_myst, 
              aes(intercept = Intercept, slope = Slope,  colour = genus, alpha = group, linetype = genus), size  = 1.2)+
  #points after, so they are on top
  scale_color_manual(name = "Nodes", 
                     values = c(mypalette_taxa_nodes_all[c(10,6,7, 1,9)]))+ 
  scale_alpha_manual(values = c(1, 1, 1))+
  theme_classic(base_size = 12)+
  scale_linetype_manual(name = "Nodes", values = as.vector(allometry_anc_all_coeffs_odont_myst$line))+
  ylab("Regression Score")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14), legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12, face = "bold"), legend.position = c(0.7,0.05), legend.direction = "horizontal", 
        legend.justification = c(0.4,0),
        legend.key = element_blank(), legend.background = element_blank())+
  guides(colour = guide_legend(nrow = 2, byrow = T, title.position = "top", keywidth = unit(2, "char"),
                               label.theme = element_text(angle = 0, face = "italic")),
         alpha = "none")
allometry_anc_all_odont_myst_ggplot

#Add phylopic
allometry_anc_all_odont_myst_ggplot <- 
  allometry_anc_all_odont_myst_ggplot +
  add_phylopic(myst, alpha = 1, x = 4, y = 0.18, ysize = 0.07, color = mypalette_groups[1])+
  add_phylopic(Kog, alpha = 1, x = 3.2, y = -0.12, ysize = 0.05, color = mypalette_odont[1])+
  add_phylopic(Lage, alpha = 1, x = 2.5, y = -0.08, ysize = 0.12, color = mypalette_odont[2])+
  add_phylopic(Phoc, alpha = 1, x = 4.4, y = 0.17, ysize = 0.06, color = mypalette_odont[3])
allometry_anc_all_odont_myst_ggplot

#Plot myst by genus
#Plot myst + anc nodes_all
#Check which nodes_all/taxa apply
allometry_anc_all_coeffs_myst <- allometry_genus_coeffs %>% filter(group == "mysticeti")
allometry_anc_all_coeffs_myst <- bind_rows(allometry_anc_all_coeffs_myst, anc_values_traj[2:6,]) %>%
  mutate(genus = coalesce(genus,group), family = coalesce(family,group))
allometry_anc_all_coeffs_myst

#Add additional grouping for linetype plot
line_nodes_myst <- rep(c(4,3,2), times = 1)
allometry_anc_all_coeffs_myst$line <- c(rep(1, times = 6), line_nodes_myst, 3, 2)
allometry_anc_all_coeffs_myst

#Plot 1 - common anc and non-Balaenopteridae
#Plot prenatal with abline to get full length - easier comparison
allometry_anc_all_myst_ggplot1  <- ggplot(allometry_anc_all, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = allometry_anc_all_coeffs_myst[c(1,4,7:8),], 
              aes(intercept = Intercept, slope = Slope,  colour = genus, alpha = group, linetype = genus), size  = 1.2)+
  #points after, so they are on top
  scale_color_manual(name = "Nodes", 
                     values = c(mypalette_taxa_nodes_all[11:12], mypalette_myst2[c(1,4)]))+ 
  scale_alpha_manual(values = c(1, 0.5))+
  theme_classic(base_size = 12)+
  scale_linetype_manual(name = "Nodes", values = rev(as.vector(allometry_anc_all_coeffs_myst$line[c(1,4,7:8)])))+
  ylab("Regression Score")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14), legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12, face = "bold"), legend.position = c(0.7,0.05), legend.direction = "horizontal", 
        legend.justification = c(0.4,0),
        legend.key = element_blank(), legend.background = element_blank())+
  guides(colour = guide_legend(nrow = 2, byrow = T, title.position = "top", keywidth = unit(2, "char"),
                               label.theme = element_text(angle = 0, face = "italic")),
         alpha = "none")
allometry_anc_all_myst_ggplot1

#Add phylopic
allometry_anc_all_myst_ggplot1 <- 
  allometry_anc_all_myst_ggplot1 +
  add_phylopic(Bala, alpha = 0.8, x = 3.6, y = 0.18, ysize = 0.08, color = mypalette_myst2[1])+
  add_phylopic(Cap, alpha = 0.8, x = 3.1, y = -0.12, ysize = 0.08, color = mypalette_myst2[4])
allometry_anc_all_myst_ggplot1

#Plot 2 - Balaenopteridae
#Plot postnatal with abline to get full length - easier comparison
allometry_anc_all_myst_ggplot2  <- ggplot(allometry_anc_all, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = allometry_anc_all_coeffs_myst[-c(1,4,7:8),], 
              aes(intercept = Intercept, slope = Slope,  colour = genus, alpha = group, linetype = genus), size  = 1.2)+
  #points after, so they are on top
  scale_color_manual(name = "Nodes", 
                     values = c(mypalette_taxa_nodes_all[13:15], mypalette_myst2[-c(1,4)]))+ 
  scale_alpha_manual(values = c(1, 0.5))+
  theme_classic(base_size = 12)+
  scale_linetype_manual(name = "Nodes", values = rev(as.vector(allometry_anc_all_coeffs_myst$line[-c(1,4,7:8)])))+
  ylab("Regression Score")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14), legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12, face = "bold"), legend.position = c(0.7,0.05), legend.direction = "horizontal", 
        legend.justification = c(0.4,0),
        legend.key = element_blank(), legend.background = element_blank())+
  guides(colour = guide_legend(nrow = 2, byrow = T, title.position = "top", keywidth = unit(2, "char"),
                               label.theme = element_text(angle = 0, face = "italic")),
         alpha = "none")
allometry_anc_all_myst_ggplot2

#Add phylopic
allometry_anc_all_myst_ggplot2 <- 
  allometry_anc_all_myst_ggplot2 +
  add_phylopic(Ball, alpha = 0.8, x = 2.6, y = -0.05, ysize = 0.06, color = mypalette_myst2[2])+
  add_phylopic(Bals, alpha = 0.8, x = 2.5, y = -0.22, ysize = 0.06, color = mypalette_myst2[3])+
  add_phylopic(Esch, alpha = 0.8, x = 4.2, y = 0.28, ysize = 0.06, color = mypalette_myst2[5])+
  add_phylopic(Meg, alpha = 0.8, x = 4.3, y = 0.1, ysize = 0.07, color = mypalette_myst2[6])
allometry_anc_all_myst_ggplot2

#Basic plot by family with all nodes_all to grab legend
allometry_anc_all_myst_legend  <- ggplot(allometry_anc_all, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = allometry_anc_all_coeffs_myst, 
              aes(intercept = Intercept, slope = Slope,  colour = genus, alpha = group, linetype = genus), size  = 1.2)+
  #points after, so they are on top
  scale_color_manual(name = "Nodes", 
                     values = c(mypalette_taxa_nodes_all[c(29:33)], mypalette_myst2))+ 
  scale_alpha_manual(values = c(1, 0.5))+
  scale_linetype_manual(name = "Nodes", values = c(3, 4 ,2, 3, 2, rep(1, times = 6)))+
  theme_classic(base_size = 12)+
  ylab("Regression Score")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14), legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12, face = "bold"), legend.position = "bottom", legend.direction = "horizontal", 
        legend.key = element_blank(), legend.background = element_blank())+
  guides(colour = guide_legend(nrow = 1, byrow = T, title.position = "top", 
                               label.theme = element_text(angle = 0, face = "italic")),
         alpha = "none")
allometry_anc_all_myst_legend

leg_anc_all_myst <- get_legend(allometry_anc_all_myst_legend)

ggarrange(allometry_anc_all_myst_ggplot1,allometry_anc_all_myst_ggplot2, 
          ncol = 2, nrow = 1, common.legend = T, legend = "bottom", legend.grob = leg_anc_all_myst)

#Plot myst by genus
#Plot myst + anc nodes interesting
allometry_anc_all_coeffs_myst_nodes <- allometry_anc_all_coeffs_myst[1:9,]

#Plot with abline to get full length - easier comparison
allometry_anc_all_myst_nodes_ggplot  <- ggplot(allometry_anc_all, aes(x = logCS, y = RegScores))+
  geom_point(size = 0, colour = "white")+
  #line on plot
  geom_abline(data = allometry_anc_all_coeffs_myst_nodes, 
              aes(intercept = Intercept, slope = Slope,  colour = genus, alpha = group, linetype = genus), size  = 1)+
  #points after, so they are on top
  scale_color_manual(name = "Nodes", 
                     values = c(mypalette_taxa_nodes_all[11:13], mypalette_taxa_nodes_all[c(1:5,8)]))+ 
  scale_alpha_manual(values = c(1, 1))+
  theme_classic(base_size = 12)+
  scale_linetype_manual(name = "Nodes", values = c(4, 3 ,2, 2, 1 ,1, 2, 1 ,1))+
  ylab("Regression Score")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14), legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12, face = "bold"), legend.position = c(0.7,0.05), legend.direction = "horizontal", 
        legend.justification = c(0.4,0),
        legend.key = element_blank(), legend.background = element_blank())+
  guides(colour = guide_legend(nrow = 2, byrow = T, title.position = "top", keywidth = unit(2, "char"),
                               label.theme = element_text(angle = 0, face = "italic")),
         alpha = "none")
allometry_anc_all_myst_nodes_ggplot 

#Add phylopic
allometry_anc_all_myst_nodes_ggplot  <- 
  allometry_anc_all_myst_nodes_ggplot  +
  add_phylopic(Bala, alpha = 0.8, x = 4, y = 0.25, ysize = 0.06, color = mypalette_taxa_nodes_all[1])+
  add_phylopic(Cap, alpha = 0.8, x = 2.6, y = -0.21, ysize = 0.05, color = mypalette_taxa_nodes_all[4])+
  add_phylopic(Ball, alpha = 0.8, x = 2.5, y = -0.1, ysize = 0.045, color = mypalette_taxa_nodes_all[2])+
  add_phylopic(Bals, alpha = 1, x = 2.5, y = -0.15, ysize = 0.04, color = mypalette_taxa_nodes_all[3])+
  add_phylopic(Esch, alpha = 0.8, x = 2.95, y = -0.2, ysize = 0.04, color = mypalette_taxa_nodes_all[5])+
  add_phylopic(Meg, alpha = 0.8, x = 4.3, y = 0.12, ysize = 0.06, color = mypalette_taxa_nodes_all[8])
allometry_anc_all_myst_nodes_ggplot 

###### 
#Next - ch. 7 - Disparity and cluster analyses