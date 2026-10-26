#############################################
# The block is loading the packages used to read and prepare the tables
#############################################
library(dplyr)      
library(tidyr)       
library(readr)       
##############################################################################
# The block is loading the packages used to set up the structure of the figures
##############################################################################
library(ggplot2)   
library(ggpubr)  
library(magick)      
library(parallel)    

##############################################################################
# The block is loading the packages used for the geomorphometric analyses
##############################################################################
library(Momocs)
library(sfsmisc)
library(geomorph)
library(morphospace)
library(mvMORPH)
#######################################################################################
# The block is loading the packages used to read and manage phylogenetic tree data
#######################################################################################
library(ape)
library(phytools)
library(nlme)
########################################################################
# The block is loading the packages used for the functional analysis
########################################################################
library(alphahull)
library(pracma)
library(MASS)
########################################################################
# The block is loading the packages used for the statistical analyses
########################################################################
library(dunn.test)
library(pairwiseAdonis)
library(knitr)
library(tibble)
library(randomForest)
library(pls)
library(forcats)
#######################################################
# The block is loading the packages used for plotting
######################################################
library(gridExtra)
library(grid)
library(RColorBrewer)
library(scales)
library(viridis)
library(ggnewscale)
library(ggrepel)
library(cowplot)
library(akima)
library(mgcv)
library(magrittr)
library(patchwork)
library(viridisLite)

###############################################################
# The block is saveguarding function to properly read the tables
###############################################################
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate
summarise <- dplyr::summarise
summarize <- dplyr::summarize

cat(" All libraries loaded\n\n")
##################################################################
# The lines are loading the custom functions used to extract wing outlines
###########################################################
source("R/MorphometricExtraction_Functions.r")
source("R/MorphoFiles_Function.r")

#################################
# The block is helperring to set exact decimals
#################################
format_fixed <- function(x, digits = 3) {
ifelse(is.na(x), NA, formatC(as.numeric(x), format = "f", digits = digits, flag = "0"))
}