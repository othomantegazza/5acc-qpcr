library(tidyverse)
library(readxl)

source("helper_functions_fluidigm.R")

fluidigms <- c(ap2_exp = "data/AP2ok_CHIP3.xlsx",
               ref_exp = "data/HK_CHIP3.xlsx") %>%
  map(., read_fluidigm)


norms <- fluidigms$ref_exp %>%
  group_by(sample_name) %>%
  summarize(norm_geom_mean = gm_mean(ct_value))

ap2_exp <- fluidigms$ap2_exp %>%
  left_join(norms)  %>%
  # subtract normalizer and take exponential to estimate expression
  # note!!! Skip calibration, is it legit????
  # So the formula is 2^-(Ct_gene - Ct_norm)
  mutate(expression = 2^(-(ct_value - norm_geom_mean))) %>%
  # Low expressed genes (Ct 999) to 0
  mutate(expression = round(expression, digits = 5))


# Save --------------------------------------------------------------------

save(ap2_exp, file = "data/ap2_fluidigm.Rdata")
