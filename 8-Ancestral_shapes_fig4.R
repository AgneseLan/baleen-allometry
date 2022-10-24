
#===========================================================#
#                                                           #
#       CURVES AND POINTS ANALYSES - MYSTICETI ONLY         #
#                                                           #
#===========================================================#


#CH.8 - Ancestral shapes prenatal and postnatal for Figure 4

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

#require(devtools)
#install_github("JeroenSmaers/evomap")
#devtools::install_github("wabarr/ggphylomorpho")
#devtools::install_github("aphanotus/borealis")

#apropos("x") lists objects with matching part of name

#MEAN SHAPES PER TAXON PER STAGE -----
#To be used in phyloPCA to estimate mean shapes

#1-Search for rows for each stage and genus
#Empty object
stage_rows <- list()

#Loop
for (s in 1:length(stages_list)){
  stage_rows[[s]] <- which(gdf_mysticeti$stage == stages_list[s])
}

#2 - Select coords for each stage
#Empty object
stage_coords <- list()

#Loop
for (s in 1:s){
  stage_coords[[s]] <- gdf_mysticeti$coords[,,stage_rows[[s]]]
}

#3 - Select genera for each stage
#Empty object
stage_genera <- list()

#Loop
for (s in 1:s){
  stage_genera[[s]] <- gdf_mysticeti$genus[stage_rows[[s]]]
}

#4 - Combine coords and genera in new gdf per stage
#Empty object
gdf_stage <- list()

#Loop
for (s in 1:s){
  gdf_stage[[s]] <- geomorph.data.frame(coords = stage_coords[[s]], genus = stage_genera[[s]])
}

#5 - Select rows for each genus in each stage
#Empty object
stage_genus_rows_1 <- list()
stage_genus_rows_2 <- list()

#Loop
for (k in 1:length(genera_list_mysticeti)){
  stage_genus_rows_1[[k]] <- which(gdf_stage[[1]]$genus == genera_list_mysticeti[k])
  stage_genus_rows_2[[k]] <- which(gdf_stage[[2]]$genus == genera_list_mysticeti[k])
}

#5 - Make arrays with mshape for each genus in each stage
#Empty object
#First number = number of landmark points, second number = 3 dimensions, third number = number of stage
coords_prenatal <- array(0, dim = c(dim(gdf_mysticeti$coords)[[1]],3,length(genera_list_mysticeti))) 
coords_postnatal <- array(0, dim = c(dim(gdf_mysticeti$coords)[[1]],3,length(genera_list_mysticeti)))

#Loop
for (k in 1:k){
  coords_prenatal[,,k] <- mshape(gdf_stage[[1]]$coords[,,stage_genus_rows_1[[k]]])
  coords_postnatal[,,k] <- mshape(gdf_stage[[2]]$coords[,,stage_genus_rows_2[[k]]])
}

#Assign dimnames as genera names
dimnames(coords_prenatal) [[3]] <- as.list(genera_list_mysticeti)
dimnames(coords_postnatal) [[3]] <- as.list(genera_list_mysticeti)

#PHYLOPCA PER STAGE - obtain ancestral shapes, only for illustration purposes ----

#Prune tree to keep only Mysticeti taxa
tree_mysticeti <- drop.tip(tree_all, setdiff(tree_all$tip.label, genera_list_mysticeti))
plot(tree_mysticeti, show.node.label = T)
nodelabels() #save tree figure

##Prenatal ----

#PCA analysis with phylogeny
#Phylogeny is used to calculate principal components of variation - GLS method
PCA_prenatal_phylo <- gm.prcomp(coords_prenatal, phy = tree_mysticeti, GLS = TRUE)

#View results
summary(PCA_prenatal_phylo)

#Plot PCA
plot(PCA_prenatal_phylo, phylo = TRUE, #add phylogeny
     main = "phyloPCA prenatal", pch = 21, col = mypalette_category[2], bg = mypalette_category[2], cex = 2, font.main = 2,
     phylo.par = list(node.txt.col  ="gray10", tip.txt.adj = c(0.5, -1)))

#Save anc shapes
anc_prenatal <- arrayspecs(PCA_prenatal_phylo[["ancestors"]], dim(gdf_mysticeti$coords)[[1]], 3)


##Postnatal ----

#PCA analysis with phylogeny
#Phylogeny is used to calculate principal components of variation - GLS method
PCA_postnatal_phylo <- gm.prcomp(coords_postnatal, phy = tree_mysticeti, GLS = TRUE)

#View results
summary(PCA_postnatal_phylo)

#Plot PCA
plot(PCA_postnatal_phylo, phylo = TRUE, #add phylogeny
     main = "phyloPCA postnatal", pch = 21, col = mypalette_category[4], bg = mypalette_category[4], cex = 2, font.main = 2,
     phylo.par = list(node.txt.col  ="gray10", tip.txt.adj = c(0.5, -1)))

#Save anc shapes
anc_postnatal <- arrayspecs(PCA_postnatal_phylo[["ancestors"]], dim(gdf_mysticeti$coords)[[1]], 3)


#PLOT SHAPES PER NODE - prenatal and postnatal ----

#Create simplified color palette for prenatal and postnatal shapes
#Get first palette
mypalette_blue_red <- brewer.pal(11,"RdBu")
image(1:11, 1, as.matrix(1:11), col = mypalette_blue_red, xlab = "Blue Red",
      ylab = "", yaxt = "n")

#Color scale - simplified rostrum, braincase and nasals
#Different colors for each stage
col_prenatal <- c(LMmodules, curvemodules)

col_prenatal <- as.factor(col_prenatal)

levels(col_prenatal) <- c(mypalette_blue_red[10], mypalette_blue_red[10], mypalette_blue_red[10], #basioccipital
                         mypalette_blue_red[10], mypalette_blue_red[10], #condyles
                         mypalette_blue_red[10], #exoccipital
                         mypalette_blue_red[10], mypalette_blue_red[10], mypalette_blue_red[10], #frontal, interparietal, jugal
                         mypalette_blue_red[8], mypalette_blue_red[9],  mypalette_blue_red[9], #maxilla, nasals
                         mypalette_blue_red[8], mypalette_blue_red[8], mypalette_blue_red[8], #palatines, premax
                         mypalette_blue_red[10], mypalette_blue_red[10], mypalette_blue_red[8]) #squamosal, socc., vomer


col_postnatal <- c(LMmodules, curvemodules)

col_postnatal <- as.factor(col_postnatal)

levels(col_postnatal) <- c(mypalette_blue_red[2], mypalette_blue_red[2], mypalette_blue_red[2], #basioccipital
                          mypalette_blue_red[2], mypalette_blue_red[2], #condyles
                          mypalette_blue_red[2], #exoccipital
                          mypalette_blue_red[2], mypalette_blue_red[2], mypalette_blue_red[2], #frontal, interparietal, jugal
                          mypalette_blue_red[4], mypalette_blue_red[3],  mypalette_blue_red[3], #maxilla, nasals
                          mypalette_blue_red[4], mypalette_blue_red[4], mypalette_blue_red[4], #palatines, postmax
                          mypalette_blue_red[2], mypalette_blue_red[2], mypalette_blue_red[4]) #squamosal, socc., vomer


#Plot shapes using vector relative to mean shape prenatal or postnatal for modern taxa
#Node numbers as tree

node7_prenatal <- c(plotRefToTarget(mshape(coords_prenatal),anc_prenatal[,,1],
                                  method="vector", radius=.00001, mag=1),
                    spheres3d(anc_prenatal[,,1], radius=.001, color = col_prenatal),
                    spheres3d(mshape(coords_prenatal), radius=.0005, color = "gray50",  alpha = 0.8, fastTransparency = T))

rgl.snapshot(filename = "Output/node7_prenatal.png") 
rgl.snapshot(filename = "Output/node7_prenatal1.png") 
clear3d()

node7_postnatal <- c(plotRefToTarget(mshape(coords_postnatal),anc_postnatal[,,1],
                                    method="vector", radius=.00001, mag=1),
                    spheres3d(anc_postnatal[,,1], radius=.001, color = col_postnatal),
                    spheres3d(mshape(coords_postnatal), radius=.0005, color = "gray50",  alpha = 0.8, fastTransparency = T))

rgl.snapshot(filename = "Output/node7_postnatal.png") 
rgl.snapshot(filename = "Output/node7_postnatal1.png") 
clear3d()

node8_prenatal <- c(plotRefToTarget(mshape(coords_prenatal),anc_prenatal[,,2],
                                    method="vector", radius=.00001, mag=1),
                    spheres3d(anc_prenatal[,,2], radius=.001, color = col_prenatal),
                    spheres3d(mshape(coords_prenatal), radius=.0005, color = "gray50",  alpha = 0.8, fastTransparency = T))

rgl.snapshot(filename = "Output/node8_prenatal.png") 
rgl.snapshot(filename = "Output/node8_prenatal1.png") 
clear3d()

node8_postnatal <- c(plotRefToTarget(mshape(coords_postnatal),anc_postnatal[,,2],
                                     method="vector", radius=.00001, mag=1),
                     spheres3d(anc_postnatal[,,2], radius=.001, color = col_postnatal),
                     spheres3d(mshape(coords_postnatal), radius=.0005, color = "gray50",  alpha = 0.8, fastTransparency = T))

rgl.snapshot(filename = "Output/node8_postnatal.png") 
rgl.snapshot(filename = "Output/node8_postnatal1.png") 
clear3d()

node9_prenatal <- c(plotRefToTarget(mshape(coords_prenatal),anc_prenatal[,,3],
                                    method="vector", radius=.00001, mag=1),
                    spheres3d(anc_prenatal[,,3], radius=.001, color = col_prenatal),
                    spheres3d(mshape(coords_prenatal), radius=.0005, color = "gray50",  alpha = 0.8, fastTransparency = T))

rgl.snapshot(filename = "Output/node9_prenatal.png") 
rgl.snapshot(filename = "Output/node9_prenatal1.png") 
clear3d()

node9_postnatal <- c(plotRefToTarget(mshape(coords_postnatal),anc_postnatal[,,3],
                                     method="vector", radius=.00001, mag=1),
                     spheres3d(anc_postnatal[,,3], radius=.001, color = col_postnatal),
                     spheres3d(mshape(coords_postnatal), radius=.0005, color = "gray50",  alpha = 0.8, fastTransparency = T))

rgl.snapshot(filename = "Output/node9_postnatal.png") 
rgl.snapshot(filename = "Output/node9_postnatal1.png") 
clear3d()

#Plot mean shape prenatal and postnatal each taxon
#Check taxon order
dimnames(coords_prenatal)[[3]]

balaena_prenatal <- spheres3d(coords_prenatal[,,1], radius=.001, color = col_prenatal)

rgl.snapshot(filename = "Output/balaena_prenatal.png") 
rgl.snapshot(filename = "Output/balaena_prenatal1.png") 
clear3d()

balaenopteraL_prenatal <- spheres3d(coords_prenatal[,,2], radius=.001, color = col_prenatal)

rgl.snapshot(filename = "Output/balaenopteraL_prenatal.png") 
rgl.snapshot(filename = "Output/balaenopteraL_prenatal1.png") 
clear3d()

balaenopteraS_prenatal <- spheres3d(coords_prenatal[,,3], radius=.001, color = col_prenatal)

rgl.snapshot(filename = "Output/balaenopteraS_prenatal.png") 
rgl.snapshot(filename = "Output/balaenopteraS_prenatal1.png") 
clear3d()

caperea_prenatal <- spheres3d(coords_prenatal[,,4], radius=.001, color = col_prenatal)

rgl.snapshot(filename = "Output/caperea_prenatal.png") 
rgl.snapshot(filename = "Output/caperea_prenatal1.png") 
clear3d()

eschrichtius_prenatal <- spheres3d(coords_prenatal[,,5], radius=.001, color = col_prenatal)

rgl.snapshot(filename = "Output/eschrichtius_prenatal.png") 
rgl.snapshot(filename = "Output/eschrichtius_prenatal1.png") 
clear3d()

megaptera_prenatal <- spheres3d(coords_prenatal[,,6], radius=.001, color = col_prenatal)

rgl.snapshot(filename = "Output/megaptera_prenatal.png") 
rgl.snapshot(filename = "Output/megaptera_prenatal1.png") 
clear3d()


balaena_postnatal <- spheres3d(coords_postnatal[,,1], radius=.001, color = col_postnatal)

rgl.snapshot(filename = "Output/balaena_postnatal.png") 
rgl.snapshot(filename = "Output/balaena_postnatal1.png") 
clear3d()

balaenopteraL_postnatal <- spheres3d(coords_postnatal[,,2], radius=.001, color = col_postnatal)

rgl.snapshot(filename = "Output/balaenopteraL_postnatal.png") 
rgl.snapshot(filename = "Output/balaenopteraL_postnatal1.png") 
clear3d()

balaenopteraS_postnatal <- spheres3d(coords_postnatal[,,3], radius=.001, color = col_postnatal)

rgl.snapshot(filename = "Output/balaenopteraS_postnatal.png") 
rgl.snapshot(filename = "Output/balaenopteraS_postnatal1.png") 
clear3d()

caperea_postnatal <- spheres3d(coords_postnatal[,,4], radius=.001, color = col_postnatal)

rgl.snapshot(filename = "Output/caperea_postnatal.png") 
rgl.snapshot(filename = "Output/caperea_postnatal1.png") 
clear3d()

eschrichtius_postnatal <- spheres3d(coords_postnatal[,,5], radius=.001, color = col_postnatal)

rgl.snapshot(filename = "Output/eschrichtius_postnatal.png") 
rgl.snapshot(filename = "Output/eschrichtius_postnatal1.png") 
clear3d()

megaptera_postnatal <- spheres3d(coords_postnatal[,,6], radius=.001, color = col_postnatal)

rgl.snapshot(filename = "Output/megaptera_postnatal.png") 
rgl.snapshot(filename = "Output/megaptera_postnatal1.png") 
clear3d()
