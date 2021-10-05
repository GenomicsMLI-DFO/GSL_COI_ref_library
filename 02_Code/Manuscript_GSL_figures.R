# Info --------------------------------------------------------------------

# Figures, tables and stats related to the manuscript
# 
# Audrey Bourret
# 2021-10-04
#

# Library -----------------------------------------------------------------

library(here)
library(tidyverse)




library(readxl)

# Internal functions
for(i in 1:length( list.files("./01_Code/00_Functions/eDNA_fct") )){
  source(file.path("./01_Code/00_Functions/eDNA_fct",  list.files("./01_Code/00_Functions/eDNA_fct")[i]))  
}

library(ggpubr)
library(ggtext)


# Data --------------------------------------------------------------------

#load(file.path("01_Raw_data/MS_All.assignments.Rdata"))
#
#write_csv(RES.total, file.path(here::here(), "01_Raw_data", "MS_AllAssignments.csv"))
#

RES.total <- read_csv(file.path(here::here(), "01_Raw_data", "MS_AllAssignments.csv"))
RES.total

