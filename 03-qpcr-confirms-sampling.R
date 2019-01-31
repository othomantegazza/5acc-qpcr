#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")


library(tidyverse)
library(magrittr)
library(cowplot)

source("helper_functions_fluidigm.R")  

chip3 <- readRDS(snakemake@input[["chip3"]])
chip4 <- readRDS(snakemake@input[["chip4"]])

# dev
# chip3 <- readRDS("data/chip3-normalized.Rds")
# chip4 <- readRDS("data/chip4-normalized.Rds")

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
           stage = case_when(stage == "stage_1" ~ "RM",
                             stage == "stage_2" ~ "IM",
                             stage == "stage_3" ~ "DM",
                             stage == "stage_4" ~ "FM"))

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

# fix gene names
rename <- c("G1L5 TAW" = "LOC_Os10g33780\n(G1L5 / TAW1)",
            "FZP" = "LOC_Os07g47330\n(FZP)",
            "LAX1" = "LOC_Os01g61480\n(LAX1)",
            "OsMADS01 LHS1" = "LOC_Os03g11614\n(LHS1)")
qpcr_imputed2 <- qpcr_imputed %>%
    ungroup() %>% 
    mutate(
        target_name = recode(target_name, !!!rename),
        stage = factor(stage, levels = c("RM", "IM", "DM", "FM")))


# Plot --------------------------------------------------------------------

p <- 
    qpcr_imputed2 %>%
    filter(target_name != "OsMADS14") %>% 
    {lineplot_fluidigm(nm = "Sampling Check", dat = .)} +
    scale_y_log10() +
    theme_minimal(base_size = 8, base_family = "Helvetica") +
    theme(axis.text.x = element_text(angle = 0,
                                     hjust = 0.5),
          strip.text = element_text(face = "italic"),
          panel.background = element_rect(colour = "black")) +
    xlab(NULL) + ggtitle(NULL)

# Add tiff figure to plot -------------------------------------------------

# # Note
# #
# # compressed with:
# # tiffcp -c zip \
# #     data-raw/FigMeristemCollect.tif \
# #     data-raw/FigMeristemCollect-compr.tif 
# 
p_img <- ggdraw() +
    draw_image("data/FigXX MeristemCollect.tif")

p_comb <- plot_grid(p_img, p,
                    labels = c("(a)", "(b)"),
                    nrow = 2,
                    rel_heights = c(1, 3),
                    label_size = 10,
                    label_fontfamily = "Helvetica")

ggsave(snakemake@output[["f1"]],
       device = cairo_pdf,
       p_comb,
       width = 178,
       height = 178,
       units = "mm")

# Log
sessionInfo()


