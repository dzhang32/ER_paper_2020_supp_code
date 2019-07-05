library(derfinder)
library(rtracklayer)
library(GenomicRanges)
library(recount)
library(tidyverse)
library(stringr)
library(IRanges)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Functions -------------------------------------------------------------------------------------------

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

# generate df with ERs varying cut off paths, split by Leo or self-generated
# leo's contains 11 brain regions (all minus the crbl)
tissue_varied_cut_off_paths_leo <- list.files("/data/recount/", recursive = T, full.names = T, pattern = "region_cuts.Rdata")
# leo's contains remaining 42 tissues (- Cells - Leukemia cell line CML - no samples passing the smarze)
tissue_varied_cut_off_paths_gtex_tissues <- list.files("/data/recount/GTEx_SRP012682/gtex_ERs_varying_cut_offs/all_cutoffs_all_chr_merged/", 
                                                       recursive = T, full.names = T, pattern = ".rda")

tissue_varied_cut_off_paths_df <- 
  data_frame(tissue_varied_cut_off_paths = c(tissue_varied_cut_off_paths_leo, tissue_varied_cut_off_paths_gtex_tissues), 
             leo_or_dz = c(rep("leo", length(tissue_varied_cut_off_paths_leo)), rep("dz", length(tissue_varied_cut_off_paths_gtex_tissues))), 
             tissue = tissue_varied_cut_off_paths %>% 
               str_replace("/region_cuts.Rdata", "") %>% 
               str_replace("_ERs_all_chrs_cutoff.*.rda", "") %>% 
               str_replace("/.*/", "")) %>% 
  filter(tissue != "brain_hippocampus")

maxgaps_to_test <- seq(0, 100, by = 10)

for(i in 1:nrow(tissue_varied_cut_off_paths_df)){
  
  tissue_varied_cut_off_path <- tissue_varied_cut_off_paths_df$tissue_varied_cut_off_paths[i]
  tissue <- tissue_varied_cut_off_paths_df$tissue[i]
  leo_or_dz <- tissue_varied_cut_off_paths_df$leo_or_dz[i]
  
  print(str_c(Sys.time(), " - loading - ", i, " - ", tissue))
  
  load(tissue_varied_cut_off_path)
  
  if(leo_or_dz == "leo"){
    
    ERs_tissue_all_cut_offs_all_chrs <- region_cuts
    remove(region_cuts)
    
  }
  
  results_path_tissue <- 
    make_results_dir(results_path = "/data/recount/GTEx_SRP012682/gtex_ERs_varying_maxgaps/", folder_name = tissue)
  
  cut_offs <- names(ERs_tissue_all_cut_offs_all_chrs)
  
  # going for only 1 to 10 and by 0.2
  cut_offs <- cut_offs[cut_offs %in% as.character(seq(1, 10, by = 0.2))]
  
  ERs_tissue_all_cut_offs_all_maxgaps <- list()
  
  for(j in seq_along(cut_offs)){
    
    cut_off_to_filter <- cut_offs[j]
    
    print(str_c(Sys.time(), " - loading - ", j, " - ", cut_off_to_filter))
    
    ERs_tissue_all_cut_offs_all_chrs_no_scaffold <- 
      ERs_tissue_all_cut_offs_all_chrs[[cut_off_to_filter]] %>% 
      keepSeqlevels(str_c("chr", c(1:22, "X", "Y", "M")), pruning.mode = "coarse") %>% 
      sort()
    
    ERs_tissue_one_cut_off_all_maxgaps <- list()
    
    for(k in seq_along(maxgaps_to_test)){
      
      maxgap_to_test <- maxgaps_to_test[k]
      
      print(str_c(Sys.time(), " - loading - ", k, " - ", maxgap_to_test))
      
      if(maxgap_to_test == 0){
        
        ERs_tissue_one_cut_off_all_maxgaps[[k]] <- ERs_tissue_all_cut_offs_all_chrs_no_scaffold
        
      }else{
        
        ERs_tissue_one_cut_off_all_maxgaps[[k]] <- 
          GenomicRanges::reduce(ERs_tissue_all_cut_offs_all_chrs_no_scaffold, min.gapwidth = maxgap_to_test)
      
      }
    
    }
    
    names(ERs_tissue_one_cut_off_all_maxgaps) <- 
      str_c(tissue, " - cutoff:", cut_off_to_filter, " - maxgap:", maxgaps_to_test)
    
    ERs_tissue_all_cut_offs_all_maxgaps[[j]] <- ERs_tissue_one_cut_off_all_maxgaps 
    
  }
  
  names(ERs_tissue_all_cut_offs_all_maxgaps) <- as.character(cut_offs)
  
  save(ERs_tissue_all_cut_offs_all_maxgaps, 
       file = str_c(
         results_path_tissue, "/",
         tissue, "_ERs_", "all_chrs_cutoff_", 
         min(as.numeric(cut_offs)), "_to_", 
         max(as.numeric(cut_offs)), "_by_", 
         as.numeric(cut_offs)[2] - as.numeric(cut_offs)[1], "_maxgap_0_to_100.rda"))
  
  rm(ERs_tissue_all_cut_offs_all_maxgaps)
  rm(ERs_tissue_all_cut_offs_all_chrs)
  rm(ERs_tissue_all_cut_offs_all_chrs_no_scaffold)
  
}

