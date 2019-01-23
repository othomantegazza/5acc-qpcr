library(tidyverse)

source("helper_functions_fluidigm.R")

chip3 <- read_csv("data/table-SX-fluidigm-chip3.csv")
chip4 <- read_csv("data/table-SX-fluidigm-chip4.csv")

norms <- c("HK04", "HK09", "HKludo01", "HKludo02")

normalize_chip <- function(chip) 
{
  norms_chip <- chip %>%
    filter(target_name %in% norms) %>%
    group_by(sample_name) %>%
    summarize(norm_geom_mean = gm_mean(ct_value))
  
  chip_exp <- chip %>%
    filter(! target_name %in% norms) %>%
    left_join(norms_chip) %>%
    # subtract normalizer and take exponential to estimate expression
    # So the formula is 2^-(Ct_gene - Ct_norm)
    mutate(expression = 2^(-(ct_value - norm_geom_mean))) %>%
    # Low expressed genes (Ct 999) to 0
    mutate(expression = round(expression, digits = 5))
}

chip3_exp <- chip3 %>% normalize_chip()
chip4_exp <- chip4 %>% normalize_chip()

save(chip3_exp, chip4_exp, file = "data/fluidigm-confim-sampling.Rdata")
