library(tidyverse)
library(stringr)
library(derfinder)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

# Functions -------------------------------------------------------------------------------------------

##### First level #####

source("/home/dzhang/projects/constraint_conservation_wd/constraint_conservation/conservation/conservation_general/conservation_general_functions_from_bw.R")
source("/home/dzhang/projects/constraint_conservation_wd/constraint_conservation/constraint/constraint_general/constraint_general_functions.R")

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

summarise_constraint_conserv <- function(ERs_w_annotation_gr_tissue_filtered_rand_regions_w_conserv_constraint, iter){
  
  ERs_w_annotation_gr_tissue_filtered_rand_regions_w_conserv_constraint_mdata <- 
  elementMetadata(ERs_w_annotation_gr_tissue_filtered_rand_regions_w_conserv_constraint) %>% 
    as.data.frame() %>% 
    as_tibble()
  
  ERs_w_annotation_gr_tissue_filtered_rand_regions_w_conserv_constraint_phast <- 
    ERs_w_annotation_gr_tissue_filtered_rand_regions_w_conserv_constraint_mdata %>% 
    group_by(theRegion, query_aggregate) %>% 
    summarise(mean_phast_cons = mean(mean_phastCons20way, na.rm = T)) %>% 
    ungroup()
  
  ERs_w_annotation_gr_tissue_filtered_rand_regions_w_conserv_constraint_CDTS <- 
    ERs_w_annotation_gr_tissue_filtered_rand_regions_w_conserv_constraint_mdata %>% 
    filter(percent_CDTS_bin_coverage >= 80) %>% 
    group_by(theRegion, query_aggregate) %>% 
    summarise(mean_CDTS = mean(mean_CDTS, na.rm = T)) %>% 
    ungroup()
  
  ERs_w_annotation_gr_tissue_filtered_rand_regions_w_conserv_constraint_phast_CDTS <- 
    ERs_w_annotation_gr_tissue_filtered_rand_regions_w_conserv_constraint_phast %>% 
    left_join(ERs_w_annotation_gr_tissue_filtered_rand_regions_w_conserv_constraint_CDTS) %>% 
    mutate(iteration = iter)
  
  return(ERs_w_annotation_gr_tissue_filtered_rand_regions_w_conserv_constraint_phast_CDTS)
  
}

# Main ------------------------------------------------------------------------------------------------

results_dir <- make_results_dir(results_path = "results/generate_randomised_intron_inter_regions/", 
                                folder_name = "intron_inter_randomised_regions_all_tissues_constraint_conserv")

intron_inter_randomised_regions_all_tissues_paths <- 
  list.files("results/generate_randomised_intron_inter_regions/intron_inter_randomised_regions_all_tissues_v92/", 
             full.names = T)

for(i in seq_along(intron_inter_randomised_regions_all_tissues_paths)){
  
  intron_inter_randomised_regions_all_tissues_path <- intron_inter_randomised_regions_all_tissues_paths[[i]]
  
  print(str_c(Sys.time(), " - ", i, " - ", intron_inter_randomised_regions_all_tissues_path))
  
  tissue_to_filter <- intron_inter_randomised_regions_all_tissues_path %>% str_replace(".*/", "")
  
  intron_inter_randomised_regions_per_tissues_paths <- 
    list.files(path = intron_inter_randomised_regions_all_tissues_path, 
              full.names = T, recursive = T)
  
  ERs_w_annotation_gr_conserv_constraint_summarised_tissue <- data_frame()
  
  for(j in seq_along(intron_inter_randomised_regions_per_tissues_paths)){
    
    intron_inter_randomised_regions_per_tissues_path <- intron_inter_randomised_regions_per_tissues_paths[j]
    
    print(str_c(Sys.time(), " - ", j, " - ", intron_inter_randomised_regions_per_tissues_paths[j]))
    
    load(intron_inter_randomised_regions_per_tissues_path)
    
    iter <- 
      (intron_inter_randomised_regions_per_tissues_path %>% 
      str_extract("random_regions.*\\.rda") %>% 
      str_replace("random_regions", "") %>% 
      str_replace("\\.rda", "") %>% 
      str_split("-") %>% 
      unlist() %>% 
      .[2] %>% 
      as.numeric())/100
    
    suppressWarnings(expr =  
                  
    ERs_w_annotation_gr_tissue_filtered_rand_regions_w_conserv_constraint <- 
      ERs_w_annotation_gr_tissue_filtered_rand_regions[[1]] %>% 
      get_conservation_score_for_regions_bw(bw_path = "/data/conservation/phastCons/hg38.phastCons20way.bw", 
                                            gr = ., summaryFun = "mean") %>% 
      get_constraint_score_for_regions(region_details_gr = ., constraint_gr = CDTS_percentile_N7794_unrelated_all_chrs_gr, get_score_percentile = F)
    
    )
    
    ERs_w_annotation_gr_conserv_constraint_summarised_permu <- 
      summarise_constraint_conserv(ERs_w_annotation_gr_tissue_filtered_rand_regions_w_conserv_constraint, iter)
    
    ERs_w_annotation_gr_conserv_constraint_summarised_tissue <- 
      ERs_w_annotation_gr_conserv_constraint_summarised_tissue %>% 
      bind_rows(ERs_w_annotation_gr_conserv_constraint_summarised_permu)
  
  }
  
  write_delim(x = ERs_w_annotation_gr_conserv_constraint_summarised_tissue %>% mutate(tissue = tissue_to_filter), 
              path = str_c("results/generate_randomised_intron_inter_regions/intron_inter_randomised_regions_all_tissues_constraint_conserv/", 
                           tissue_to_filter, "_rand_conserv_constraint_summarised.csv"), delim = ",")
  
}

