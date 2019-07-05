library(tidyverse)
library(stringr)
library(forcats)
library(ggpubr)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(regioneR)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

load(file = "results/annotate_ERs/ERs_optimised_cut_off_max_gap_all_tissues_w_annot_df.rda")

OMIM_gtex_tissue_name_formatting <- read_delim("raw_data/gtex_tissue_name_formatting/OMIM_gtex_tissue_name_formatting.csv", delim = ",")

# Functions -------------------------------------------------------------------------------------------

##### First level #####

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_genomic_state.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/get_random_length_matched_regions.R")

get_intron_intergenic_space <- function(hg38_genome_mask_list, ensembl_grch38_v92_genomic_state, intron_inter){
  
  ensembl_grch38_v92_genomic_state_intron_inter <- 
    ensembl_grch38_v92_genomic_state[[1]][(ensembl_grch38_v92_genomic_state[[1]]$theRegion == intron_inter) & 
                                            (width(ensembl_grch38_v92_genomic_state[[1]]) > 3)]
  
  ensembl_grch38_v92_genomic_state_intron_inter_black_list_hits <- 
    findOverlaps(query = hg38_genome_mask_list$mask, subject = ensembl_grch38_v92_genomic_state_intron_inter, 
                 minoverlap = 1, maxgap = -1)
  
  ensembl_grch38_v92_genomic_state_intron_inter_no_black_list <- 
    ensembl_grch38_v92_genomic_state_intron_inter[-subjectHits(ensembl_grch38_v92_genomic_state_intron_inter_black_list_hits)]
  
  # remove unnecessary column save space
  ensembl_grch38_v92_genomic_state_intron_inter_no_black_list$tx_id <- NULL
  ensembl_grch38_v92_genomic_state_intron_inter_no_black_list$tx_name <- NULL
  ensembl_grch38_v92_genomic_state_intron_inter_no_black_list$gene <- NULL
  names(ensembl_grch38_v92_genomic_state_intron_inter_no_black_list) <- NULL
  
  return(ensembl_grch38_v92_genomic_state_intron_inter_no_black_list)
  
}

# generates output folder to store coloc results, does not overwrite the folder if it already exists
#' @param results_path path to folder 
#' @param folder_name name of folder to make
#' @return the full path to the new folder
make_results_dir <- function(results_path, folder_name){
  
  results_dir_path <- str_c(results_path, "/", folder_name)
  
  if(!dir.exists(results_dir_path)){
    
    dir.create(results_dir_path)
    
  }else {
    
    print(str_c(results_dir_path, " directory already exists.."))
    
  }
  
  return(results_dir_path)
  
}

# Main ------------------------------------------------------------------------------------------------

results_dir <- 
  make_results_dir(results_path = "/home/dzhang/projects/OMIM_wd/results/generate_randomised_intron_inter_regions/", 
                   folder_name = "intron_inter_randomised_regions_all_tissues_v92")

hg38_genome_mask_list <- getGenomeAndMask("hg38", mask = NULL)

ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific <-
  ERs_w_annotation_all_tissues %>%
  filter(width > 3, !str_detect(tissue, "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast"),
         ensembl_grch38_v92_region_annot != "exon, intergenic, intron")

ERs_w_annotation_intron_inter_gr <- 
  ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific %>% 
  filter(ensembl_grch38_v92_region_annot %in% c("intron", "intergenic")) %>% 
  dplyr::select(seqnames, start, end, ensembl_grch38_v92_region_annot, tissue) %>% 
  as.data.frame() %>% 
  toGRanges()

rm(list = c("ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific", "ERs_w_annotation_all_tissues"))

# get full intron or intergenic space from genomic state 
# removing those regions that overlap with known blacklist regions 
ensembl_grch38_v92_genomic_state <- 
  generate_genomic_state(output_path = "/data/references/ensembl/genomic_state/v92/ensembl_grch38_v92_genomic_state.rda")

ensembl_grch38_v92_genomic_state_intron_space <- 
  get_intron_intergenic_space(hg38_genome_mask_list, ensembl_grch38_v92_genomic_state, intron_inter = "intron")

ensembl_grch38_v92_genomic_state_intergenic_space <- 
  get_intron_intergenic_space(hg38_genome_mask_list, ensembl_grch38_v92_genomic_state, intron_inter = "intergenic")

tissues_to_filter <- ERs_w_annotation_intron_inter_gr$tissue %>% unique()

for(i in seq_along(tissues_to_filter)){
  
  tissue_to_filter <- tissues_to_filter[[i]]
  
  print(str_c(Sys.time(), " - ", tissue_to_filter))
  
  ERs_w_annotation_intron_inter_gr_tissue_filtered <- 
    ERs_w_annotation_intron_inter_gr %>% 
    subset(tissue == tissue_to_filter)
  
  results_dir_tissue <- 
    make_results_dir(results_path = results_dir, folder_name = tissue_to_filter)
  
  for(j in seq_along(c("intron", "intergenic"))){
    
    intron_inter_to_filter <- c("intron", "intergenic")[j]
    
    print(str_c(Sys.time(), " - ", intron_inter_to_filter))
    
    ERs_w_annotation_gr_tissue_intron_inter_filtered <- 
      ERs_w_annotation_intron_inter_gr_tissue_filtered %>% 
      subset(ensembl_grch38_v92_region_annot == intron_inter_to_filter)
    
    results_dir_tissue_intron_inter <- 
      make_results_dir(results_path = results_dir_tissue, folder_name = intron_inter_to_filter)
   
    # iterate across the number of permuations you want, generating 100 aggregates each time
    for(k in 1:(10000/100)){
      
      permus <- str_c((k * 100 - 99), "-", (k * 100))
      
      print(str_c("generating random regions for permutations: ", permus))
      
      if(intron_inter_to_filter == "intron"){
        
        ERs_w_annotation_gr_tissue_filtered_rand_regions <- 
          get_random_length_matched_regions(query_gr = ERs_w_annotation_gr_tissue_intron_inter_filtered, 
                                            subject_gr = ensembl_grch38_v92_genomic_state_intron_space,
                                            n_permuations = 100, n_aggregates = 100,
                                            min_rand_subject_size = 10000)
        
      }else if(intron_inter_to_filter == "intergenic"){
        
        ERs_w_annotation_gr_tissue_filtered_rand_regions <- 
          get_random_length_matched_regions(query_gr = ERs_w_annotation_gr_tissue_intron_inter_filtered, 
                                            subject_gr = ensembl_grch38_v92_genomic_state_intergenic_space,
                                            n_permuations = 100, n_aggregates = 100,
                                            min_rand_subject_size = 10000)
        
      }
      
      save(ERs_w_annotation_gr_tissue_filtered_rand_regions, 
           file = str_c(results_dir_tissue_intron_inter, "/", tissue_to_filter, "_", intron_inter_to_filter, "_random_regions", permus, ".rda"))
      
    }
    
  }
  
}

