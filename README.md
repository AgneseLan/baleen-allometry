# Testing heterochrony in baleen whales 🐳 🦐 🔬📈
### Connecting skull shape ontogeny and evolution of different feeding adaptations in Mysticeti using 3D geometric morphometrics 

Authors: [Agnese Lanzetti](mailto:agnese.lanzetti@gmail.com?subject=[GitHub]%20Ontogeny%20Baleen%20Whales%20Paper%20Code), Roberto Portela Miguez,
Vincent Fernandez, Anjali Goswami

To cite the paper: 

Available at: https://github.com/AgneseLan/baleen-allometry

If using any of this code or data please cite the paper above and this repo

To cite this repo:



## Data :floppy_disk: 

The data are provided in the Data folder. Mesh files (PLY) needed to test postioning when importing landmarks are available at Phenome10k (https://www.phenome10k.org/). 

- __Landmark data__: *pts folder* <br />
Text files with landmark coordinates for each specimen in PTS format. 

- __Surface data__: *ply folder* <br />
Empty folder where mesh files from Phenome10k need to be saved to reproduce the code.

- __Specimens' classifiers, landmark/curves lists__: *absent_curves_myst.csv, absent_LMs_myst.csv, curves_all.csv, LMs_all.csv, specimens_myst.csv* <br />
Spreadsheets with additional inforation for analyses: list of absent landmarks and curves, list of curves, list of landmarks, classifiers for specimens (specimen names, ID, group, family, genera, common name, length, bizygomatic width, age, growth stage, feeding mode).

- __Phylogentic tree__: *tree_myst.txt* <br />
Tree with branch lenghts for the taxa used in the analyses in Nexus format. Branch lenghts and topology from McGowen et al., 2020, Syst. Biol. 69, 479-501.

- __Reference mesh for plotting__: *refmesh_myst.zip* <br />
Reduced mesh in PLY format used for plotting landmarks.

- __Silhouettes of taxa for plots__: *b.bona.png, b.physalus.png, balaena.png, caperea.png, eschrichtius.png, kogia.png, lagenorhynchus.png, megaptera.png, phocoena.png*

## Analysis :computer:
In this repository you will find raw data (.csv and data files) and code for analyses (code supplied as .R files)

📁 Data

As described above. Meshes used to collect and test landmarks available on Phenome10k (https://www.phenome10k.org/). 

⌨ Code for analyses - .R files

*1-Import-resample-slide.R, 1-Slider3d_2.R, 2-Absent_bones.R, 3-GPA_PCA.R, 4-Allometry.R, 5-Allometry_genus.R, 6-Anc_state_allometry.R, 7-Disparity_Cluster.R, 8-Ancestral_shapes_fig4.R*

Code files are numbered providing the order the analyses need to be performed in.
Before running analyses, save Data folder in the same directory as the R project. This will allow to import the data as detailed in the code provided.

## License 📃
This project is licensed under the MIT License - see the LICENSE.md file for details

## Session Info 📋
For reproducibility purposes, here is the output of devtools::session_info() used to perform the analyses in the publication.

```
R version 4.1.3 (2022-03-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252  
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C  LC_TIME=English_United States.1252  

attached base packages:
[1] grid  stats  graphics  grDevices utils  datasets  methods  base  

other attached packages:
[1] Anthropometry_1.17  viridis_0.6.2  viridisLite_0.4.1  grateful_0.0.3  emmeans_1.8.1-1  
[6] mcp_0.3.2 reshape2_1.4.4  rray_0.1.0.9000  scales_1.2.1  car_3.1-1  
[11] carData_3.0-5  evomap_0.0.0.9000  phytools_1.2-0  maps_3.4.0  gridExtra_2.3  
[16] png_0.1-7  rphylopic_0.3.0  ggplotify_0.1.0  ggpubr_0.4.0  ggthemes_4.2.4  
[21] borealis_2021.02.04 ggfortify_0.4.14  ggphylomorpho_0.1.0 gginnards_0.1.1  ggrepel_0.9.1  
[26] RColorBrewer_1.1-3  magick_2.7.3  SURGE_0.1.0  devtools_2.4.5  usethis_2.1.6  
[31] abind_1.4-5  geiger_2.0.10  ape_5.6-2  qgraph_1.9.2  EMMLi_0.0.3  
[36] paleomorph_0.1.4  Rvcg_0.21  geomorph_4.0.4  Matrix_1.5-1  rgl_0.110.2  
[41] RRPP_1.3.1  Morpho_2.10  forcats_0.5.2  stringr_1.4.1  dplyr_1.0.10  
[46] purrr_0.3.5  readr_2.1.3  tidyr_1.2.1  tibble_3.1.8  ggplot2_3.3.6  
[51] tidyverse_1.3.2  

loaded via a namespace (and not attached):
  [1] utf8_1.2.2  tidyselect_1.2.0  htmlwidgets_1.5.4  combinat_0.0-8  additivityTests_1.1-4.1
  [6] biclust_2.0.3  munsell_0.5.0  codetools_0.2-18  interp_1.1-3  miniUI_0.1.1.1  
 [11] withr_2.5.0  colorspace_2.0-3  knitr_1.40  rstudioapi_0.14  stats4_4.1.3  
 [16] ggsignif_0.6.4  mnormt_2.1.1  optimParallel_1.0-2  depth_2.1-1.1  coda_0.19-4  
 [21] vctrs_0.5.0  generics_0.1.3  circular_0.4-95  clusterGeneration_1.3.7 xfun_0.34  
 [26] fastcluster_1.2.3  archetypes_2.2-0.1  R6_2.5.1  doParallel_1.0.17  cachem_1.0.6  
 [31] gridGraphics_0.5-1  assertthat_0.2.1  promises_1.2.0.1  nnet_7.3-18  googlesheets4_1.0.1  
 [36] gtable_0.3.1  processx_3.7.0  phangorn_2.10.0  rlang_1.0.6  scatterplot3d_0.3-42  
 [41] splines_4.1.3  rstatix_0.7.0  gargle_1.2.1  broom_1.0.1  checkmate_2.1.0  
 [46] yaml_2.3.6  modelr_0.1.9  backports_1.4.1  httpuv_1.6.6  Hmisc_4.7-1  
 [51] tools_4.1.3  psych_2.2.9  lavaan_0.6-12  shapes_1.2.6  gridBase_0.4-7  
 [56] ellipsis_0.3.2  sessioninfo_1.2.2  Rcpp_1.0.9  plyr_1.8.7  base64enc_0.1-3  
 [61] ps_1.7.1  prettyunits_1.1.1  rpart_4.1.19  deldir_1.0-6  pbapply_1.5-0  
 [66] urlchecker_1.0.1  deSolve_1.34  haven_2.5.1  cluster_2.1.4  colorRamps_2.3.1  
 [71] fs_1.5.2  crul_1.3  magrittr_2.0.3  data.table_1.14.4  reprex_2.0.2  
 [76] googledrive_2.0.0  mvtnorm_1.1-3  matrixStats_0.62.0  flexclust_1.4-1  pkgload_1.3.0  
 [81] evaluate_0.17  patchwork_1.1.2  hms_1.1.2  mime_0.12  xtable_1.8-4  
 [86] jpeg_0.1-9  readxl_1.4.1  compiler_4.1.3  crayon_1.5.2  htmltools_0.5.3  
 [91] corpcor_1.6.10  later_1.3.0  tzdb_0.3.0  Formula_1.2-4  expm_0.999-6  
 [96] lubridate_1.8.0  DBI_1.1.3  dbplyr_2.2.1  subplex_1.8  MASS_7.3-58.1  
[101] boot_1.3-28  cli_3.4.1  quadprog_1.5-8  parallel_4.1.3  igraph_1.3.5  
[106] pkgconfig_2.0.3  numDeriv_2016.8-1.1  foreign_0.8-83  xml2_1.3.3  foreach_1.5.2  
[111] pbivnorm_0.6.0  minpack.lm_1.2-2  estimability_1.4.1  rvest_1.0.3  yulab.utils_0.0.5  
[116] bezier_1.1.2  callr_3.7.2  digest_0.6.30  ICGE_0.4.2  httpcode_0.3.0  
[121] rmarkdown_2.17  cellranger_1.1.0  fastmatch_1.1-3  htmlTable_2.4.1  curl_4.3.3  
[126] modeltools_0.2-23  shiny_1.7.2  gtools_3.9.3  lifecycle_1.0.3  nlme_3.1-160  
[131] glasso_1.11  jsonlite_1.8.3  fansi_1.0.3  pillar_1.8.1  loo_2.5.1  
[136] lattice_0.20-45  fastmap_1.1.0  httr_1.4.4  plotrix_3.8-2  pkgbuild_1.3.1  
[141] survival_3.4-0  glue_1.6.2  remotes_2.4.2  FNN_1.1.3.1  fdrtool_1.2.17  
[146] iterators_1.0.14  class_7.3-20  nnls_1.4  stringi_1.7.8  profvis_0.3.7  
[151] latticeExtra_0.6-30  memoise_2.0.1  renv_0.16.0  
```
