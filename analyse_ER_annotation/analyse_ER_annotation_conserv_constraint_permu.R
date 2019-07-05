library(tidyverse)
library(stringr)
library(forcats)
library(ggpubr)
library(regioneR)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

load(file = "results/annotate_ERs/ERs_optimised_cut_off_max_gap_all_tissues_w_annot_df.rda")

OMIM_gtex_tissue_name_formatting <- read_delim("raw_data/gtex_tissue_name_formatting/OMIM_gtex_tissue_name_formatting.csv", delim = ",")

funicane_tstats <- read_delim("/home/dzhang/data/finucane_2018/tstats/GTEx.tstat.tsv", delim = "\t")

# Functions -------------------------------------------------------------------------------------------

source("/home/dzhang/projects/constraint_conservation_wd/constraint_conservation/conservation/conservation_general/conservation_general_functions_from_bw.R")

get_phast_cons_20_format_constraint_conserv_ERs <- function(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific, finucane_brain_specific_genes, tissue_to_filter){
  
  tissue_intron_inter_mean_CDTS <- 
  ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific %>% 
    filter(ensembl_grch38_v92_region_annot %in% c("intron", "intergenic"), percent_CDTS_bin_coverage >= 80, tissue == tissue_to_filter) %>% 
    group_by(tissue, ensembl_grch38_v92_region_annot) %>% 
    summarise(mean_CDTS = mean(mean_CDTS))
    
  tissue_intron_inter_mean_phastcons20 <- 
    ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific %>% 
    filter(ensembl_grch38_v92_region_annot %in% c("intron", "intergenic"), tissue == tissue_to_filter) %>% 
    dplyr::select(seqnames, start, end, tissue, ensembl_grch38_v92_region_annot) %>% 
    as.data.frame() %>% 
    toGRanges() %>% 
    get_conservation_score_for_regions_bw(bw_path = "/data/conservation/phastCons/hg38.phastCons20way.bw", gr = ., summaryFun = "mean") %>% 
    toDataframe() %>% 
    group_by(tissue, ensembl_grch38_v92_region_annot) %>% 
    summarise(mean_phast_cons_20_ER = mean(mean_phastCons20way, na.rm = T))
  
  tissue_intron_inter_mean_CDTS_phastcons20 <- 
    tissue_intron_inter_mean_CDTS %>% 
    left_join(tissue_intron_inter_mean_phastcons20)
  
  return(tissue_intron_inter_mean_CDTS_phastcons20)
  
}



plot_ERs_vs_rand_intron_inter <- function(randomised_intron_inter_conserv_constraint_per_tissue, tissue_intron_inter_mean_CDTS_phastcons20){
  
  randomised_intron_inter_conserv_constraint_per_tissue_to_plot <- 
  randomised_intron_inter_conserv_constraint_per_tissue %>% 
    left_join(tissue_intron_inter_mean_CDTS_phastcons20, by = c("theRegion" = "ensembl_grch38_v92_region_annot", "tissue" = "tissue"))
  
  ggplot(randomised_intron_inter_conserv_constraint_per_tissue_to_plot, aes(x = mean_phast_cons)) +
    geom_density() + 
    geom_vline(aes(xintercept = mean_phast_cons_20_ER)) +
    facet_wrap(~ theRegion)
  
  randomised_intron_inter_CDTS_plot <- 
    ggplot(randomised_intron_inter_conserv_constraint_per_tissue_to_plot, 
           aes(x = mean_CDTS)) +
    geom_area(aes(y = (..count..)/sum(..count..)), fill = "grey", colour = "black", stat = "bin", bins = 50) +
    geom_vline(aes(xintercept = mean_CDTS_ER), linetype = 2, colour = "#CD534CFF") +
    scale_x_continuous(name = "Mean CDTS") +
    scale_y_continuous(name = "Density") +
    facet_grid(theRegion ~ ., scales = "free_y") +
    theme_pubr(border = T, legend = "right") 
  
  randomised_intron_inter_phastcons_plot <- 
    ggplot(randomised_intron_inter_conserv_constraint_per_tissue_to_plot, 
           aes(x = mean_phast_cons)) +
    geom_area(aes(y = (..count..)/sum(..count..)), fill = "grey", colour = "black", stat = "bin", bins = 50) +
    geom_vline(aes(xintercept = mean_phast_cons_20_ER), linetype = 2, colour = "#CD534CFF") +
    scale_x_continuous(name = "Mean CDTS") +
    scale_y_continuous(name = "Density") +
    facet_grid(theRegion ~ ., scales = "free_y") +
    theme_pubr(border = T, legend = "right") 
 
  return(randomised_intron_inter_CDTS_plot) 
  
}

# Main ------------------------------------------------------------------------------------------------

feature_type_order_colors <- 
  data_frame(feature_type_order = c("exon", "exon, intron", "exon, intergenic", 
                                    "exon, intergenic, intron", 
                                    "intergenic", "intron"), 
             feature_type_colour = get_palette("jco", 9)[c(6, 1, 5, 3, 9, 4)])

ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific <-
  ERs_w_annotation_all_tissues %>%
  filter(width > 3, !str_detect(tissue, "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast"),
         ensembl_grch38_v92_region_annot != "exon, intergenic, intron")

finucane_brain_specific_genes <- get_list_brain_specific_genes(funicane_tstats, top_propor = 0.10)

randomised_intron_inter_conserv_constraint_per_tissue_paths <- 
  list.files(path = "results/generate_randomised_intron_inter_regions/intron_inter_randomised_regions_all_tissues_constraint_conserv/", 
             full.names = T)

for(i in seq_along(randomised_intron_inter_conserv_constraint_per_tissue_paths)){
  
  randomised_intron_inter_conserv_constraint_per_tissue_path <- randomised_intron_inter_conserv_constraint_per_tissue_paths[i]
  
  tissue_to_filter <- 
    randomised_intron_inter_conserv_constraint_per_tissue_path %>% 
    str_replace(".*/", "") %>% 
    str_replace("_rand_conserv_constraint_summarised.csv", "")
  
  randomised_intron_inter_conserv_constraint_per_tissue <- read_delim(randomised_intron_inter_conserv_constraint_per_tissue_path, delim = ",")
  
  tissue_intron_inter_mean_CDTS_phastcons20 <- 
    get_phast_cons_20_format_constraint_conserv_ERs(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific %>% 
                                                        filter(!is.na(p_annot_junc_count_samples_split_read_annot)), finucane_brain_specific_genes, 
                                                      tissue_to_filter)
  
  plot_ERs_vs_rand_intron_inter(randomised_intron_inter_conserv_constraint_per_tissue, tissue_intron_inter_mean_CDTS_phastcons20)
  

}



