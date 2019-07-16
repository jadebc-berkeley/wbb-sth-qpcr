#######################################
# WASH Benefits Bangladesh STH KK qPCR validation

# configure data directories
# source base functions
# load libraries
#######################################
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(washb)
library(VennDiagram)
library(venneuler)
library(irr)
library(clusrank)
library(assertthat)

# define directories
data_dir = "~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/"
# data_dir = "~/Dropbox/WASH-B-STH-Add-on/TFGH/Data/RData/public/"
fig_dir = paste0(here::here(),"/3-figures/")
tab_dir = paste0(here::here(),"/4-tables/")

# source base functions  
source(paste0(here::here(), "/2-fig-tab-scripts/0-base-plot-functions.R"))
source(paste0(here::here(), "/2-fig-tab-scripts/0-base-table-functions.R"))

