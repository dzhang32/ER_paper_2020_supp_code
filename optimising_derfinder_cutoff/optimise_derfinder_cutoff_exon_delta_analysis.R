library(tidyverse)
library(stringr)
library(ggpubr)
library(forcats)
library(recount)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------
 
exon_details_all_tissues <- 
  read_delim("results/optimise_derfinder_cut_off/exon_delta_details_all_tissues_cut_off_1_10_0.2_maxgap_0_100_10.csv", delim = ",")

gtex_tissue_name_formatting <- read_delim("raw_data/gtex_tissue_name_formatting/OMIM_gtex_tissue_name_formatting.csv", delim = ",")

# Functions -------------------------------------------------------------------------------------------

##### First level #####

plot_exon_delta_diff_tissues <- function(exon_details_all_tissues){
  
  exon_details_all_tissues_w_propor <- 
  exon_details_all_tissues %>% 
    mutate(propor_exon_delta_eq_0 = num_exon_delta_eq_0/num_ERs_tested_for_exon_delta, 
           propor_ERs_for_ab_1_ER_overlapping_1_exon = num_ERs_for_ab_1_ER_overlapping_1_exon/num_ERs_tested_for_exon_delta) 
  
  propor_ERs_for_ab_1_ER_overlapping_1_exon_plot <- 
    ggplot(exon_details_all_tissues_w_propor %>% 
             filter(tissue %in% c("brain_cerebellum", "hippocampus", "amygdala", "hypothalamus", "putamenbasalganglia", "muscle_skeletal", 
                                  "heart_left_ventricle", "cells_transformed_fibroblasts", "whole_blood")) %>% 
             mutate(tissue = tissue %>% factor() %>% fct_relevel(c("brain_cerebellum", "hippocampus", "amygdala", "hypothalamus", "putamenbasalganglia",
                                                                   "muscle_skeletal", "heart_left_ventricle", "cells_transformed_fibroblasts", "whole_blood"))), 
           aes(x = cut_off, y = propor_ERs_for_ab_1_ER_overlapping_1_exon)) +
    geom_line(aes(colour = as.character(maxgap) %>% factor() %>% fct_relevel(as.character(seq(0, 100, 10))))) +
    scale_x_continuous(name = "Cut off") +
    scale_y_continuous(name = "Proportion multiple ERs overlapping 1 exon") +
    scale_colour_manual("Maxgap", values = c(get_palette("Blues", 14)[5:14], "red")) +
    facet_wrap(~ tissue) +
    theme_bw()
  
  ggplot(exon_details_all_tissues_w_propor %>% 
           filter(tissue %in% c("brain_cerebellum", "hippocampus", "amygdala", "hypothalamus", "putamenbasalganglia", "muscle_skeletal", 
                                "heart_left_ventricle", "cells_transformed_fibroblasts", "whole_blood", "testis")) %>% 
           mutate(tissue = tissue %>% factor() %>% fct_relevel(c("brain_cerebellum", "hippocampus", "amygdala", "hypothalamus", "putamenbasalganglia",
                                                                 "muscle_skeletal", "heart_left_ventricle", "cells_transformed_fibroblasts", "whole_blood"))), 
         aes(x = cut_off, y = total_ERs)) +
    geom_line(aes(colour = as.character(maxgap) %>% factor() %>% fct_relevel(as.character(seq(0, 100, 10))))) +
    scale_x_continuous(name = "Cut off") +
    scale_y_continuous(name = "Proportion multiple ERs overlapping 1 exon") +
    scale_colour_manual("Maxgap", values = c(get_palette("Blues", 14)[5:14], "red")) +
    facet_wrap(~ tissue) +
    theme_bw()
  
  exon_delta_median_plot <- 
    ggplot(exon_details_all_tissues_w_propor %>% 
           filter(tissue %in% c("brain_cerebellum", "hippocampus", "amygdala", "hypothalamus", "putamenbasalganglia", "muscle_skeletal", 
                                "heart_left_ventricle", "cells_transformed_fibroblasts", "whole_blood")) %>% 
           mutate(tissue = tissue %>% factor() %>% fct_relevel(c("brain_cerebellum", "hippocampus", "amygdala", "hypothalamus", "putamenbasalganglia", "muscle_skeletal", 
                                                                 "heart_left_ventricle", "cells_transformed_fibroblasts", "whole_blood"))), 
         aes(x = cut_off, y = exon_delta_median)) +
    geom_line(aes(colour = as.character(maxgap) %>% factor() %>% fct_relevel(as.character(seq(0, 100, 10))))) +
    scale_x_continuous(name = "Cut off") +
    scale_y_continuous(name = "Exon delta median") +
    scale_colour_manual("Maxgap", values = c(get_palette("Blues", 14)[5:14], "red")) +
    facet_wrap(~ tissue) +
    theme_bw()
  
  propor_exon_delta_eq_0_plot <- 
    ggplot(exon_details_all_tissues_w_propor %>% 
             filter(tissue %in% c("brain_cerebellum", "hippocampus", "amygdala", "hypothalamus", "putamenbasalganglia", "muscle_skeletal", 
                                  "heart_left_ventricle", "cells_transformed_fibroblasts", "whole_blood")) %>% 
             mutate(tissue = tissue %>% factor() %>% fct_relevel(c("brain_cerebellum", "hippocampus", "amygdala", "hypothalamus", "putamenbasalganglia",
                                                                   "muscle_skeletal", "heart_left_ventricle", "cells_transformed_fibroblasts", "whole_blood"))), 
           aes(x = cut_off, y = propor_exon_delta_eq_0)) +
    geom_line(aes(colour = as.character(maxgap) %>% factor() %>% fct_relevel(as.character(seq(0, 100, 10))))) +
    scale_x_continuous(name = "Cut off") +
    scale_y_continuous(name = "Propor exon delta = 0") +
    scale_colour_manual("Maxgap", values = c(get_palette("Blues", 14)[5:14], "red")) +
    facet_wrap(~ tissue) +
    theme_bw()
  
  num_exon_delta_eq_0_plot <- 
    ggplot(exon_details_all_tissues_w_propor %>% 
             filter(tissue %in% c("brain_cerebellum", "hippocampus", "amygdala", "hypothalamus", "putamenbasalganglia", "muscle_skeletal", 
                                  "heart_left_ventricle", "cells_transformed_fibroblasts", "whole_blood")) %>% 
             mutate(tissue = tissue %>% factor() %>% fct_relevel(c("brain_cerebellum", "hippocampus", "amygdala", "hypothalamus", "putamenbasalganglia",
                                                                   "muscle_skeletal", "heart_left_ventricle", "cells_transformed_fibroblasts", "whole_blood"))), 
           aes(x = cut_off, y = num_exon_delta_eq_0)) +
    geom_line(aes(colour = as.character(maxgap) %>% factor() %>% fct_relevel(as.character(seq(0, 100, 10))))) +
    scale_x_continuous(name = "Cut off") +
    scale_y_continuous(name = "Num exon delta = 0") +
    scale_colour_manual("Maxgap", values = c(get_palette("Blues", 14)[5:14], "red")) +
    facet_wrap(~ tissue) +
    theme_bw()

  return(list(propor_ERs_for_ab_1_ER_overlapping_1_exon_plot, exon_delta_median_plot, propor_exon_delta_eq_0_plot, num_exon_delta_eq_0_plot))
  
}

get_min_exon_delta_max_exon_delta_eq_0 <- function(exon_details_all_tissues){
  
  # get maxgap for the highest num exon delta = 0
  all_tissues_optimised_cutoff_maxgap <- 
  exon_details_all_tissues %>% 
    group_by(tissue) %>% 
    filter(exon_delta_median == min(exon_delta_median)) %>% 
    filter(num_exon_delta_eq_0 == max(num_exon_delta_eq_0)) %>% 
    filter(maxgap == min(maxgap)) %>% 
    filter(cut_off == min(cut_off))
  
  return(all_tissues_optimised_cutoff_maxgap)

}

get_list_ERs_each_tissue_optimal_cut_off_maxgap <- function(ER_min_exon_delta_max_exon_delta_eq_0){

  ERs_per_tissue_maxgap_paths_df <- 
    data_frame(ERs_per_tissue_maxgap_paths = list.files(path = "/data/recount/GTEx_SRP012682/gtex_ERs_varying_maxgaps/", full.names = T, recursive = T), 
             tissue = 
               ERs_per_tissue_maxgap_paths %>% 
               str_replace(".*/", "") %>% 
               str_replace("_ERs_all_chrs_cutoff_1_to_10_by_0.2_maxgap_0_to_100.rda", ""))
  
  list_ERs_each_tissue_optimal_cut_off_maxgap <- list()
  names_ERs_all_tissues_optimised_cut_off_maxgap <- character()
  
  for(i in 1:nrow(ER_min_exon_delta_max_exon_delta_eq_0)){
    
    optimal_cut_off <- ER_min_exon_delta_max_exon_delta_eq_0$cut_off[i]
    optimal_maxgap <- ER_min_exon_delta_max_exon_delta_eq_0$maxgap[i]
    tissue_to_filter <- ER_min_exon_delta_max_exon_delta_eq_0$tissue[i]
    tissue_path <- 
      ERs_per_tissue_maxgap_paths_df %>% 
      filter(tissue == tissue_to_filter) %>% 
      .[["ERs_per_tissue_maxgap_paths"]]
    
    print(str_c(Sys.time(), " - ", i, " - ", tissue_to_filter))

    load(tissue_path)
    
    name_ERs_one_tissue_optimised_cut_off_maxgap <- str_c(tissue_to_filter, " - cutoff:", optimal_cut_off, " - maxgap:", optimal_maxgap)
    
    names_ERs_all_tissues_optimised_cut_off_maxgap[[i]] <- name_ERs_one_tissue_optimised_cut_off_maxgap

    list_ERs_each_tissue_optimal_cut_off_maxgap[[i]] <- 
      ERs_tissue_all_cut_offs_all_maxgaps[[as.character(optimal_cut_off)]][[name_ERs_one_tissue_optimised_cut_off_maxgap]]

    
  }
  
  names(list_ERs_each_tissue_optimal_cut_off_maxgap) <- names_ERs_all_tissues_optimised_cut_off_maxgap
  
  return(list_ERs_each_tissue_optimal_cut_off_maxgap)
  
}

##### Second level #####

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_txDb_from_gtf.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/query_biomart.R")

get_non_overlapping_exons <- function(ensembl_grch38_v87_TxDb){
  
  ensembl_grch38_v87_exons_gr <- 
    ensembl_grch38_v87_TxDb %>% exons(columns = c("EXONNAME", "GENEID"))
  
  ensembl_grch38_v87_exons_gr_marked_overlapping <- 
    mark_overlapping_genes_gr(gr_1 = ensembl_grch38_v87_exons_gr, gr_2 = ensembl_grch38_v87_exons_gr, identical_gr = T, 
                              maxgap = -1L, minoverlap = 1L)
  
  ensembl_grch38_v87_exons_gr_non_overlapping <- 
    ensembl_grch38_v87_exons_gr_marked_overlapping[ensembl_grch38_v87_exons_gr_marked_overlapping$overlap_gr2 == F]
  
  return(ensembl_grch38_v87_exons_gr_non_overlapping)
  
}

##### Third level #####

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/mark_overlapping_genes_gr.R")

# Main ------------------------------------------------------------------------------------------------

##### plot exon delta metrics across different tissues to check optimisation #####

exon_delta_plots <- plot_exon_delta_diff_tissues(exon_details_all_tissues)

##### get list of ERs for each tissue using the optimal cut off #####

ER_min_exon_delta_max_exon_delta_eq_0 <- get_min_exon_delta_max_exon_delta_eq_0(exon_details_all_tissues)

list_ERs_each_tissue_optimal_cut_off_maxgap <- get_list_ERs_each_tissue_optimal_cut_off_maxgap(ER_min_exon_delta_max_exon_delta_eq_0)

# Save data -------------------------------------------------------------------------------------------

save(list_ERs_each_tissue_optimal_cut_off_maxgap,  
     file = "results/optimise_derfinder_cut_off/ERs_optimal_cut_off_maxgap_all_tissues.rda")

write_delim(ER_min_exon_delta_max_exon_delta_eq_0, 
            "results/optimise_derfinder_cut_off/exon_delta_details_optimised_maxgap_cutoff.csv", 
            delim = ",")
