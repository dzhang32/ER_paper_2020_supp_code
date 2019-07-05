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

# Functions -------------------------------------------------------------------------------------------

merge_ERs_diff_chrs_per_tissue <- function(tissues_cut_off_chrs_paths){
  
  for(k in seq_along(tissues_cut_off_chrs_paths)){
    
    tissues_cut_off_chrs_path <- tissues_cut_off_chrs_paths[k]
    
    load(tissues_cut_off_chrs_path)
    
    if(k == 1){
      
      ERs_tissue_single_cut_off_all_chrs <- ERs_tissue_varying_cut_off
      
    }else{
      
      ERs_tissue_single_cut_off_all_chrs <- 
        c(ERs_tissue_single_cut_off_all_chrs, ERs_tissue_varying_cut_off)
      
    }

  }
  
  return(ERs_tissue_single_cut_off_all_chrs)
  
}

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

tissues_paths <- list.files("/data/recount/GTEx_SRP012682/gtex_ERs_varying_cut_offs/by_chr/", full.names = T)

for(i in seq_along(tissues_paths)){
  
  tissues_path <- tissues_paths[i]
  tissue <- tissues_path %>% str_replace("/.*/", "")
  results_path_tissue <- make_results_dir(results_path = "/data/recount/GTEx_SRP012682/gtex_ERs_varying_cut_offs/all_cutoffs_all_chr_merged/", 
                                          folder_name = tissue)
  
  print(str_c(Sys.time(), " - ", i, " - ", "loading ERs by chromsome for: ", tissue))
  
  # arrange by lowest to highest cutoff for ease 
  tissues_cut_off_paths_df <- 
    data_frame(tissues_cut_off_paths = list.files(tissues_path, full.names = T), 
               cutoff = tissues_cut_off_paths %>% 
                 str_replace("/.*/", "") %>% 
                 str_replace("cutoff_", "") %>% 
                 as.double()) %>% 
    arrange(cutoff)
  
  ERs_tissue_all_cut_offs_all_chrs <- list()
  
  for(j in 1:nrow(tissues_cut_off_paths_df)){
    
    tissues_cut_off_path <- tissues_cut_off_paths_df$tissues_cut_off_paths[j]
    cutoff <- tissues_cut_off_paths_df$cutoff[j]
    
    print(str_c(Sys.time(), " - ", j, " - ", "loading cutoff: ", cutoff))
    
    tissues_cut_off_chrs_paths <- list.files(tissues_cut_off_path, full.names = T)
    
    ERs_tissue_single_cut_off_all_chrs <- merge_ERs_diff_chrs_per_tissue(tissues_cut_off_chrs_paths)
    
    ERs_tissue_single_cut_off_all_chrs_sorted <- ERs_tissue_single_cut_off_all_chrs %>% sortSeqlevels() %>% sort()
    
    ERs_tissue_all_cut_offs_all_chrs[[j]] <- ERs_tissue_single_cut_off_all_chrs_sorted
    
  }

  names(ERs_tissue_all_cut_offs_all_chrs) <- as.character(tissues_cut_off_paths_df$cutoff)
  
  save(ERs_tissue_all_cut_offs_all_chrs, 
       file = str_c(results_path_tissue, "/", tissue, "_ERs_", "all_chrs_cutoff_", 
                    min(tissues_cut_off_paths_df$cutoff), "_to_", 
                    max(tissues_cut_off_paths_df$cutoff), "_by_", 
                    tissues_cut_off_paths_df$cutoff[2] - tissues_cut_off_paths_df$cutoff[1], ".rda"))
  
}
