library(tidyverse)
library(magrittr)
library(cowplot)


source("helper_functions_fluidigm.R")  

chip3 <- readRDS("data/chip3-normalized.Rds")
chip4 <- readRDS("data/chip4-normalized.Rds")


# Select marker genes -----------------------------------------------------

genes <- c("FZP|BFL1|SGDP7", "G1L5 TAW", "LAX1", "OsMADS01 LHS1")

# Load file ---------------------------------------------------------------

qpcr_exp <- 
  bind_rows(chip3, chip4) %>%
  filter(target_name %in% genes) %>% 
  mutate(target_name = case_when(target_name == "FZP|BFL1|SGDP7" ~ "FZP",
                                 TRUE ~ target_name)) %>% 
  scale_tidy_fluidigm() %>%
  # filter(species != "Osj") %>% 
  mutate(species = case_when(species == "Or" ~ "O. rufipogon",
                             species == "Osi" ~ "O. sativa indica",
                             species == "Osj" ~ "O. sativa japonica",
                             species == "Ob" ~ "O. barthii",
                             species == "Og" ~ "O. glaberrima"),
         species = factor(species,
                          levels = c("O. rufipogon",
                                     "O. sativa indica",
                                     "O. sativa japonica",
                                     "O. barthii",
                                     "O. glaberrima")),
         stage = case_when(stage == "stage_1" ~ "Stage 1: RM",
                           stage == "stage_2" ~ "Stage 2: IM",
                           stage == "stage_3" ~ "Stage 3: DM",
                           stage == "stage_4" ~ "Stage 4: Floret"))
         
# Impute 0/ND -------------------------------------------------------------

# What is a the minimal expression besides 0?
# btw 0 must be intended as Non Detected.

# set it to one order of magnitude smaller than
# the minimum that was analyzed for a given amplicon

impute_zero <- function(expr) {
  min_expr <- 
    expr[expr > 0] %>%
    min()
  return(min_expr/10)
} 


qpcr_imputed <- 
  qpcr_exp %>% 
  group_by(target_name) %>%
  mutate(expression = case_when(expression == 0 ~ impute_zero(expression),
                                TRUE ~ expression)) 
  

# Plot --------------------------------------------------------------------

p <- 
  qpcr_imputed %>%
  filter(target_name != "OsMADS14") %>% 
  {lineplot_fluidigm(nm = "Sampling Check", dat = .)} +
  scale_y_log10() +
  theme(strip.text = element_text(face = "italic")) +
  labs(title = "qPCR expression of selected marker genes [semi-log plot]")


pdf("fig/suppl-fig-qpcr-confirms-sampling.pdf",
    height = 6.5, width = 6.2,
    paper = "a4")

p %>% print()

dev.off()

# Add tiff figure to plot -------------------------------------------------

# # Note
# #
# # compressed with:
# # tiffcp -c zip \
# #     data-raw/FigMeristemCollect.tif \
# #     data-raw/FigMeristemCollect-compr.tif 
# 
# p_img <- ggdraw() +  
#   draw_image("../data-raw/FigXX MeristemCollect.tif")
# 
# p_comb <- plot_grid(p_img, p,
#                     labels = c("1.", "2."),
#                     nrow = 2,
#                     rel_heights = c(1, 2))
# 
# pdf("../fig/suppl-fig-qpcr-confirms-sampling.pdf",
#     height = 10, width = 6.2,
#     paper = "a4")
# 
# p_comb %>%
#   ggdraw() %>%
#   print()
# 
# dev.off()
