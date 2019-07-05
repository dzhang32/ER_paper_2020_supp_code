library(derfinder)
library(rtracklayer)
library(GenomicRanges)
library(recount)
library(tidyverse)
library(stringr)
library(ggpubr)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

dis_between_ERs_all_tissues <- 
  read_delim("results/optimise_derfinder_cut_off/dis_between_ERs_amyg_hipp_crbl_blood.csv", delim = ",")

# Functions -------------------------------------------------------------------------------------------

tissues <- dis_between_ERs_all_tissues %>% filter(!is.na(tissue)) %>% .[["tissue"]] %>% unique()

tissues_corrected_all <- c()

for(i in 1:(length(tissues)-1)){
  
  tissue_1 <- tissues[i]
  tissue_2 <- tissues[i+1]
  
  print(str_c(tissue_1, " - ", tissue_2))
  
  index_start_1 <- which(dis_between_ERs_all_tissues$tissue == tissue_1) %>% min()
  index_start_2 <- which(dis_between_ERs_all_tissues$tissue == tissue_2) %>% min()
  
  tissues_corrected <- rep(tissue_1, times =(index_start_2-index_start_1))
  
  tissues_corrected_all <- c(tissues_corrected_all, tissues_corrected)
  
  
}


# Main ------------------------------------------------------------------------------------------------

distribution_dis_between_ERs_overlap_exons <- 
  ggplot(dis_between_ERs_all_tissues %>% filter(type == "multiple_ERs_overlapping_1_exon", cut_off %in% c(0.2, seq(1, 10))), 
         aes(x = dis_between_ERs)) +
  geom_step(aes(y=..y.., 
                colour = as.character(cut_off) %>% factor() %>% fct_relevel(as.character(c(0.2, seq(1, 10))))), stat="ecdf") +
  scale_y_continuous(name = "Cumulative proportion") +
  scale_x_log10(name = "Distance between ERs", breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  scale_colour_manual(name = "Cut off", values = get_palette("Blues", 11)) + 
  facet_wrap(~tissue) +
  theme_bw()

ggsave(plot = distribution_dis_between_ERs_overlap_exons, filename = "distribution_dis_between_ERs_overlap_exons_amyg_hipp_crbl_blood.png", 
       path = "results/optimise_derfinder_cut_off/", width = 14, height = 10, dpi = 600)
