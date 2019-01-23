library(tidyverse)
library(grid)
library(gridExtra)

source("helper_functions_fluidigm.R")  

chip3 <- readRDS("data/chip3-normalized.Rds")
chip4 <- readRDS("data/chip4-normalized.Rds")


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
         stage = case_when(stage == "stage_1" ~ "Stage 1: RM",
                           stage == "stage_2" ~ "Stage 2: IM",
                           stage == "stage_3" ~ "Stage 3: DM",
                           stage == "stage_4" ~ "Stage 4: Floret"))



# Plot -----------------------------------------------------


plts <- 
  dat %>% 
  {lineplot_fluidigm(nm = "",
                    dat = .,
                    alpha = .8)} +
  theme(strip.text = element_text(size = 8,
                                  face = "italic"))

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


pdf("fig/Figure_S8-suppl-fig-fluidigm-ap2.pdf",
    height = 9, width = 9,
    paper = "a4")
plts %>% print()
dev.off()


  