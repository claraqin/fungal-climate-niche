# R code for Qin et al. on fungal niches and climate sensitivity
# 0-setup.R: Load libraries, constants, and custom functions


# Load libraries
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(glmnet)
library(raster)
library(sp)
library(sf)
library(rnaturalearth)
library(speedyseq)
library(doParallel)
library(pracma)
library(viridis)
library(grid)
library(gridExtra)
library(caret)
library(egg)
library(cowplot)
library(neonUtilities)
library(pracma)
library(rgdal)
library(sf)
library(plotbiomes)
library(grid)
library(rnaturalearth)
library(plotbiomes)
library(maxnet)
library(enmSdm)
theme_set(theme_bw())

# Set constants
N_CORES <- 16 # socs-stats has a maximum of 88 cores

# Load custom functions
source("R/utils.R")

nested_combine <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

# Set custom themes
theme_custom <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

theme_custom_map <- theme_custom +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

