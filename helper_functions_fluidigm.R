# library(tidyverse)
# library(gridExtra)
# library(grid)

#################################################################################

#' Estimate Geometric Means
#' 
#' for normalizers
#' 
#' 
#' function from here:
#' https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in

gm_mean <-  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


# Load fluidigm data ------------------------------------------------------


read_fluidigm <- function(path, sheet = 1) {
  out <- read_excel(path, sheet = sheet) %>%
    select(sampleName, GeneName,
           LOC_Name, Value, Type__1) %>%
    dplyr::rename(sample_name = sampleName,
           target_name = GeneName,
           locus_id = LOC_Name,
           ct_value = Value,
           type = Type__1) %>%
    mutate(ct_value = as.numeric(ct_value)) %>%
    # remove standard curves
    filter(sample_name != "H20",
           !grepl("Mix", sample_name))
  out 
}


# scale and tidy -----------------------------------------------------

scale_tidy_fluidigm <- function(dat) {
  
  # Store descriptive names for stages
  
  stage_names <- c(N1 = "stage_1",
                   N2 = "stage_2",
                   N3 = "stage_3",
                   N4 = "stage_4")
  
  # Store descriptive names for species
  
  species_names <- c(B = "Ob",
                     G = "Og",
                     I = "Osi",
                     J = "Osj",
                     R = "Or")
  
  # Scale gene expression and scale gene expression separating 
  
  dat <- dat %>%
    # make a variable that encodes the species
    mutate(species = str_sub(sample_name, 1, 1)) %>%
    # make a variable that encodes for the stage
    mutate(stage = str_sub(sample_name, 3, 4)) %>%
    # make a variable that encode replicate
    mutate(replicate = str_sub(sample_name, 6, 6)) %>%
    # scale expression by gene
    group_by(target_name) %>%
    mutate(expr_scale_gene = scale(expression)) %>%
    ungroup() %>%
    # scale expression by gene end species
    group_by(target_name, species) %>%
    mutate(expr_scale_gene_spec = scale(expression),
           expression_NA = case_when(expression == 0 ~ NA_real_,
                                     TRUE ~ expression),
           expr_scale_gene_spec_0_10 = scales::rescale(expression_NA,
                                                       to = c(0.1, 10),
                                                       # from = range(c(0, expression_NA),
                                                       #              na.rm = TRUE,
                                                       #              finite = TRUE)
                                                       )
           ) %>%
    ungroup() %>%
    mutate(expr_scale_gene_spec = ifelse(is.na(expr_scale_gene_spec),
                                         0,
                                         expr_scale_gene_spec)) %>%
    # last, save descriptive names for stages and species
    mutate(stage = stage_names[stage]) %>%
    # mutate(stage = factor(stage, 
    #                       levels = c("Rachis_Meristem",
    #                                  "Primary_Branch_Meristem",
    #                                  "Spikelet_meristem",
    #                                  "Young_Flower_Development"))) %>%
    mutate(species = species_names[species])
  
  return(dat)
}


# Prepare for heatmap -----------------------------------------------------


prepare_for_heat <- function(dat, scaled_val = expr_scale_gene_spec){
  
  scaled_val <- enquo(scaled_val)
  
  dat <- dat %>%
    # use median of spec-scaled
    mutate(locus_id = paste(locus_id, target_name)) %>%
    group_by(locus_id, species, stage) %>%
    summarise(median_expr = median(!!scaled_val)) %>%
    ungroup() %>% 
    # spread for heatmap
    mutate(species_stage = paste(species, stage, sep = "_")) %>%
    select(locus_id, species_stage, median_expr) %>% 
    spread(key = species_stage, value = median_expr) %>%
    # reorder columns
    # select(locus_id,
    #        Osj_stage_1, Osj_stage_2, Osj_stage_3, Osj_stage_4,
    #        Ob_stage_1,  Ob_stage_2,  Ob_stage_3,  Ob_stage_4,
    #        Og_stage_1,  Og_stage_2,  Og_stage_3,  Og_stage_4,
    #        Or_stage_1,  Or_stage_2,  Or_stage_3,  Or_stage_4,
    #        Osi_stage_1, Osi_stage_2, Osi_stage_3, Osi_stage_4) %>%
    select(locus_id,
           Osj_stage_1,
           Ob_stage_1,  
           Og_stage_1,  
           Or_stage_1,  
           Osi_stage_1,
           Osj_stage_2, 
           Ob_stage_2,  
           Og_stage_2,  
           Or_stage_2,  
           Osi_stage_2, 
           Osj_stage_3,
           Ob_stage_3,  
           Og_stage_3,  
           Or_stage_3,  
           Osi_stage_3,
           Osj_stage_4,
           Ob_stage_4,
           Og_stage_4,
           Or_stage_4,
           Osi_stage_4) %>%
    # select(target_name,
    #        Japonica_Rachis_Meristem,
    #        Japonica_Primary_Branch_Meristem,
    #        Japonica_Spikelet_meristem,
    #        Japonica_Young_Flower_Development,
    #        Barthii_Rachis_Meristem,
    #        Barthii_Primary_Branch_Meristem,
    #        Barthii_Spikelet_meristem,
    #        Barthii_Young_Flower_Development,
    #        Glaberrima_Rachis_Meristem,
    #        Glaberrima_Primary_Branch_Meristem,
    #        Glaberrima_Spikelet_meristem,
    #        Glaberrima_Young_Flower_Development,
    #        Rufipogon_Rachis_Meristem,
    #        Rufipogon_Primary_Branch_Meristem,
    #        Rufipogon_Spikelet_meristem,
    #        Rufipogon_Young_Flower_Development,
    #        Indica_Rachis_Meristem,
    #        Indica_Primary_Branch_Meristem,
    #        Indica_Spikelet_meristem,
    #        Indica_Young_Flower_Development) %>%
    as.data.frame(.)
  
  return(dat)
}

prepare_for_heat_spec <- function(dat){
  
  dat <- dat %>%
    # use median of spec-scaled
    mutate(locus_id = paste(locus_id, target_name)) %>%
    group_by(locus_id, species, stage) %>%
    summarise(median_expr = median(expr_scale_gene_spec)) %>%
    ungroup() %>% 
    # spread for heatmap
    mutate(species_stage = paste(species, stage, sep = "_")) %>%
    select(locus_id, species_stage, median_expr) %>% 
    spread(key = species_stage, value = median_expr) %>%
    # reorder columns
    select(locus_id,
           Osj_stage_1, Osj_stage_2, Osj_stage_3, Osj_stage_4,
           Ob_stage_1,  Ob_stage_2,  Ob_stage_3,  Ob_stage_4,
           Og_stage_1,  Og_stage_2,  Og_stage_3,  Og_stage_4,
           Or_stage_1,  Or_stage_2,  Or_stage_3,  Or_stage_4,
           Osi_stage_1, Osi_stage_2, Osi_stage_3, Osi_stage_4) %>%
    # select(locus_id,
    #        Osj_stage_1,
    #        Ob_stage_1,  
    #        Og_stage_1,  
    #        Or_stage_1,  
    #        Osi_stage_1,
    #        Osj_stage_2, 
    #        Ob_stage_2,  
    #        Og_stage_2,  
    #        Or_stage_2,  
    #        Osi_stage_2, 
    #        Osj_stage_3,
    #        Ob_stage_3,  
    #        Og_stage_3,  
    #        Or_stage_3,  
    #        Osi_stage_3,
    #        Osj_stage_4,
    #        Ob_stage_4,
    #        Og_stage_4,
    #        Or_stage_4,
    #        Osi_stage_4) %>%
    # select(target_name,
    #        Japonica_Rachis_Meristem,
    #        Japonica_Primary_Branch_Meristem,
    #        Japonica_Spikelet_meristem,
    #        Japonica_Young_Flower_Development,
    #        Barthii_Rachis_Meristem,
    #        Barthii_Primary_Branch_Meristem,
    #        Barthii_Spikelet_meristem,
    #        Barthii_Young_Flower_Development,
    #        Glaberrima_Rachis_Meristem,
    #        Glaberrima_Primary_Branch_Meristem,
  #        Glaberrima_Spikelet_meristem,
  #        Glaberrima_Young_Flower_Development,
  #        Rufipogon_Rachis_Meristem,
  #        Rufipogon_Primary_Branch_Meristem,
  #        Rufipogon_Spikelet_meristem,
  #        Rufipogon_Young_Flower_Development,
  #        Indica_Rachis_Meristem,
  #        Indica_Primary_Branch_Meristem,
  #        Indica_Spikelet_meristem,
  #        Indica_Young_Flower_Development) %>%
  as.data.frame(.)
  
  return(dat)
}


# Plot RNAseq and fluidigms -----------------------------------------------

#' Plot both
#' 
#' requires tidyverse grid and gridextra

plot_both <- function(id_fluidigm,
                      rnaseq_dat,
                      fluidigm_dat) {
  print(id_fluidigm)
  
  fluidigm_dat <- fluidigm_dat %>% 
    filter(locus_id == id_fluidigm)

  locus_id <- unique(fluidigm_dat$locus_id)
  target_name <- unique(fluidigm_dat$target_name)
    
  p_fluid <- ggplot(fluidigm_dat,
                    aes(x = stage,
                        y = expression)) +
    geom_point() +
    geom_smooth(aes(group = 1), se = FALSE) +
    # stat_summary(fun.y=mean, colour="red", geom="line") +
    ggtitle(paste(locus_id, target_name)) +
    facet_wrap(facets = "species", ncol = 5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  p_rnaseq <- ggplot(rnaseq_dat %>% filter(locus_id == id_fluidigm), 
                     aes(x = stage,
                         y = `normalized expression`))  +
    geom_point() +
    geom_smooth(aes(group = 1), se = FALSE) +
    ggtitle(id_fluidigm) +
    facet_wrap(facets = "species", ncol = 5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # tt <- textGrob(paste("id", unique(rnaseq_dat$DESCRIPTION[rnaseq_dat$locus_id == id])))
  lay <- rbind(c(1, 1, 1, 1, 1),
               c(1, 1, 1, 1, 1),
               c(1, 1, 1, 1, 1),
               c(2, 2, 2, 2, 2),
               c(2, 2, 2, 2, 2))
  grid.arrange(grobs = list(p_fluid, p_rnaseq),
               layout_matrix = lay)
}


# Line plot fluidigm ------------------------------------------------------

lineplot_fluidigm <- function(nm, dat,
                              alpha = .5,
                              expr_val = expression)
  {
  expr_val <- enquo(expr_val)
  
  p <- ggplot(dat, aes(x = stage,
                       y = !!expr_val)) +
    geom_point(size = 2,
               col = "darkgrey",
               alpha = alpha) +
    # geom_smooth(aes(x = stage %>%
    #                   as_factor(.) %>%
    #                   as.numeric(.)),
    #             se = FALSE,
    #             colour = "black") +
    stat_summary(aes(x = stage %>%
                       as_factor(.) %>%
                       as.numeric(.)),
                 fun.y = mean,
                 geom="line",
                 size = 1.2,
                 linejoin = "round") +
    facet_grid(target_name ~ species,
               scales = "free_y") +
    # facet_grid(species ~ target_name,
    #            scales = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_text(hjust = 0,
                                     vjust = .5,
                                     angle = 270),
          strip.text = element_text(size = 6)) +
    labs(title = nm,
         y = "Relative Expression")
  return(p)
}

# old ---------------------------------------------------------------------



#' load and tidy fluidigm data
#' 
#' requires tidyverse
#' 
#' # test
#' library(tidyverse)
#' path <- "data/data9696 puce1essaiV1.csv"
#' tst <- load_fluidigm(path)

load_fluidigm_csv <- function(path) {
  dat <- read.csv2(path,
                      skip = 11,
                      stringsAsFactors = FALSE) %>%
    # the last two rows 
    .[, 1:12] %>%
    # Give descriptive names
    rename_all(.funs = ~make.names(names = c("Chamber", "Sample_Name", "Sample_Type",
                                             "Sample_rConc", "target_name", "target_type",
                                             "Ct_Value", "Ct_Quality", "Ct_Call", 
                                             "Ct_Threshold", "tm_in_range", "tm_out_range"))) %>%
    # routine stuff
    as_tibble() 
  
  if(sum(!complete.cases(dat))> 0) print("NA detected")

  # remove NA
  dat <- dat %>% na.omit()
  
  # Tidy data ---------------------------------------------------------------
  
  
  # How many measures per target?
  table(dat$target_name)
  
  # remove negative control (water?)
  # No, some have Ct, what to do with them?
  summary(dat %>%
            filter(Sample_Name == "H20") %>%
            select(Ct_Value))
  # for now remove.
  dat <- dat %>% filter(Sample_Name != "H20")
  
  # Remove standard curves?
  # Do all sample in standard curves contain the word "Mix" in their names? Yes
  dat <- dat %>% filter(!grepl("Mix", Sample_Name))
  
  # How many measurement now?
  table(dat$target_name)
  table(dat$Sample_Name)
  table(dat$Sample_Type)
  
  # These are descriptive names for stages
  stage_names <- c(N1 = "Rachis Meristem",
                   N2 = "Primary Branch Meristem",
                   N3 = "Spikelet meristem",
                   N4 = "Young Flower Development")
  
  # These are descriptive names for species
  species_names <- c(B = "Barthii",
                     G = "Glaberrima",
                     I = "Sativa Indica",
                     J = "Sativa Japonica",
                     R = "Rufipogon")
  
  # We must store stages and species names in separate variable, 
  # They must be ordered for plotting
  dat <- dat %>%
    # make a variable that encodes for the stage
    mutate(stage = str_sub(Sample_Name, 3, 4)) %>%
    mutate(stage = unname(stage_names[stage])) %>%
    # Order is chronological
    mutate(stage = factor(stage, 
                          levels = c("Rachis Meristem",
                                     "Primary Branch Meristem",
                                     "Spikelet meristem",
                                     "Young Flower Development"))) %>%
    # make a variable that encodes the species
    mutate(species = str_sub(Sample_Name, 1, 1)) %>%
    mutate(species = unname(species_names[species])) %>%
    # Order is arbitrary for plotting
    mutate(species = factor(species,
                            levels = c("Sativa Japonica",
                                       "Barthii",
                                       "Glaberrima",
                                       "Rufipogon",
                                       "Sativa Indica")))
  
  return(dat)
}


#' Normalize Fluidigm
#' 
#' Takes Geometric mean of normalizer and estimates expression
#' as 2^deltaCt
#' 
#' test
#' path <- "data/data9696 puce1essaiV1.csv" 
#' tst <- load_fluidigm(path)
#' tst <- normalize_fluidigm(tst, c("ACT2", "HK09", "HK04"))

normalize_fluidigm <- function(dat, normalizers)
{
  if(!all(normalizers %in% dat$target_name)) {
    stop("At least one normalizer is not contained in the target_name column of the fluidigm dataset")
  }
  norms <- dat %>% filter(target_name %in% normalizers)
  
  table(norms$target_name)
  table(norms$Sample_Name)
  
  # Take the geomoetric mean for every sample
  # Check that the function for geometric mean works properly!
  norms <- norms %>%
    group_by(Sample_Name) %>%
    summarize(norm_geom_mean = gm_mean(Ct_Value))
  
  # Delta Ct: Correct for normalizers and take exponential
  
  dat <- dat %>%
    # merge back normalizers to dataset
    left_join(norms, by = "Sample_Name") %>%
    # subtract normalizer and take exponential to estimate expression
    # note!!! Skip calibration, is it legit????
    # So the formula is 2^-(Ct_gene - Ct_norm)
    mutate(expression = 2^(-(Ct_Value - norm_geom_mean))) %>%
    # Low expressed genes (Ct 999) to 0
    mutate(expression = round(expression, digits = 5))
  
  return(dat)
} 
