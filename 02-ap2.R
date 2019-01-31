#!/usr/bin/env Rscript

# set log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

library(tidyverse)
library(grid)
library(gridExtra)

source("helper_functions_fluidigm.R")  

chip3 <- readRDS(snakemake@input[["chip3"]])
chip4 <- readRDS(snakemake@input[["chip4"]])

# dev
# chip3 <- readRDS("data/chip3-normalized.Rds")
# chip4 <- readRDS("data/chip4-normalized.Rds")

# select ap2 genes -------------------------------------------------------

ap2s <- c("AP37|OsERF3", "EREBP129 - OsRAV2", "ERF33", "ERF93 / OsERF1")

ap2_exp <- 
    bind_rows(chip3, chip4) %>% 
    filter(target_name %in% ap2s)

# Dot plot ---------------------------------------------------------------

dat <-  
    ap2_exp %>%
    scale_tidy_fluidigm() %>%
    filter(species != "Osj") %>% 
    mutate(species = case_when(species == "Or" ~ "O. rufipogon",
                               species == "Osi" ~ "O. sativa indica",
                               # species == "Osj" ~ "O. sativa japonica",
                               species == "Ob" ~ "O. barthii",
                               species == "Og" ~ "O. glaberrima"),
           species = factor(species,
                            levels = c("O. rufipogon",
                                       "O. sativa indica",
                                       # "O. sativa japonica",
                                       "O. barthii",
                                       "O. glaberrima")),
           stage = case_when(stage == "stage_1" ~ "RM",
                             stage == "stage_2" ~ "IM",
                             stage == "stage_3" ~ "DM",
                             stage == "stage_4" ~ "FM"))

# fix gene names
rename <- c("ERF33" = "LOC_Os04g46400\n(ERF33)",
            "EREBP129 - OsRAV2" = "LOC_Os04g46400\n(AP2/EREBP129)",
            "AP37|OsERF3" = "LOC_Os01g58420\n(ERF3)",
            "ERF93 / OsERF1" = "LOC_Os04g46220\n(ERF93)")
dat2 <- dat %>% mutate(
    target_name = recode(target_name, !!!rename),
    stage = factor(stage, levels = c("RM", "IM", "DM", "FM")))

# Plot -----------------------------------------------------


plts <- 
    dat2 %>% 
    {lineplot_fluidigm(nm = "",
                       dat = .,
                       alpha = .8)} +
    theme(strip.text = element_text(size = 8,
                                    face = "italic"))

gp <- plts +
    xlab(NULL) +
    theme_minimal(base_size = 8, base_family = "Helvetica") +
    theme(axis.text.x = element_text(angle = 0,
                                     hjust = 0.5),
          strip.text = element_text(face = "italic"),
          panel.background = element_rect(colour = "black"))


# +
#   theme(text = element_text(size = 8),
#         strip.text = element_text(size = 4))

# p <- cowplot::plot_grid(plts[[2]], plts[[3]],
#                         nrow = 2,
#                         rel_heights = c(3, 4) + 1.5,
#                         labels = c("A", "C")) %>%
#   cowplot::plot_grid(., plts[[1]],
#                      labels = c("", "B")) %>%
#   cowplot::add_sub(., str_wrap("qPCR confirms the behaviour of selected
#                                genes of cluster 4 and cluster 5.
#                                We have also measured gene expressio in two
#                                additional developmental stages. Stage 1 is a Rachis
#                                Meristem, Stage2 is a Branch Meristem, Stage 3
#                                is a Spikelet Meristem, Stage 4 is a Developing
#                                Spikelet."),
#                    size = 11) %>%
#   cowplot::ggdraw()

ggsave(snakemake@output[["f1"]],
       device = cairo_pdf,
       gp,
       width = 178,
       height = 178,
       units = "mm")

# Log
sessionInfo()
