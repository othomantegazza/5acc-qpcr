library(tidyverse)
library(readxl)
source("helper_functions_fluidigm.R")

# read data ---------------------------------------------------------------

chip3 <- read_fluidigm("data/Fluidgmdata_publi.xlsx", sheet = "CHIP3")
chip4 <- read_fluidigm("data/Fluidgmdata_publi.xlsx", sheet = "CHIP4")

# Select targets and normalizers ------------------------------------------

norms <- c("HK04", "HK09", "HKludo01", "HKludo02")
targets <- c("FZP|BFL1|SGDP7", "LAX1", "G1L5 TAW", "OsMADS01 LHS1",
             "EREBP129 - OsRAV2", "AP37|OsERF3", "ERF33", 
             "ERF93 / OsERF1")

keep <- c(norms, targets)

# filter fluidigm files ---------------------------------------------------

chip3 <- chip3 %>% filter(target_name %in% keep)
chip4 <- chip4 %>% filter(target_name %in% keep)

# check -------------------------------------------------------------------

chip3$target_name %>% unique()
chip4$target_name %>% unique()

# save as csv -------------------------------------------------------------

chip3 %>% write_csv("data/table-SX-fluidigm-chip3.csv")
chip4 %>% write_csv("data/table-SX-fluidigm-chip4.csv")
