library(derfinder)
library(rtracklayer)
library(GenomicRanges)
library(recount)
library(tidyverse)
library(stringr)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Functions -------------------------------------------------------------------------------------------

download_hipp_cerebellum_fibroblast_bw_urls <- function(){
  
  # to speed things up initially, whilst all samples are downloading, first download the cerebellum and fibroblasts to complete optimisation of cut-offs
  GTEx_sample_bw_urls <- download_study(project = "SRP012682", type = "samples", download = F)
  
  gtex_metadata <-  all_metadata('gtex')
  
  gtex_metadata_hipp_crbl_fib <- 
    subset(gtex_metadata, 
           smtsd %in% c("Brain - Cerebellar Hemisphere", 
                        "Brain - Cerebellum",
                        "Brain - Hippocampus", 
                        "Cells - Transformed fibroblasts"))
  
  gtex_metadata_hipp_crbl_fib_paths <- 
    gtex_metadata_hipp_crbl_fib$bigwig_file %>% 
    str_c("http://duffel.rail.bio/recount/SRP012682/bw/", .)
  
  for(i in seq_along(gtex_metadata_hipp_crbl_fib_paths)){
    
    print(str_c(i, " - ", gtex_metadata_hipp_crbl_fib$bigwig_file[i]))
    
    download.file(url = gtex_metadata_hipp_crbl_fib_paths[i], 
                  destfile = str_c("/data/recount/GTEx_SRP012682/gtex_bigWigs/hipp_crbl_fib_tmp//", gtex_metadata_hipp_crbl_fib$bigwig_file[i]), 
                  method = "wget", quiet = T)
    
  }
  
}

download_all_gtex_tissue_bw <- function(){
  
  # to speed things up initially, whilst all samples are downloading, first download the cerebellum and fibroblasts to complete optimisation of cut-offs
  GTEx_sample_bw_urls <- download_study(project = "SRP012682", type = "samples", download = F)
  
  gtex_metadata <-  all_metadata('gtex')
  
  gtex_metadata_bw_paths <- 
    gtex_metadata$bigwig_file %>% 
    str_c("http://duffel.rail.bio/recount/SRP012682/bw/", .)
  
  for(i in seq_along(gtex_metadata_bw_paths)){
    
    print(str_c(i, " - ", gtex_metadata$bigwig_file[i]))
    
    download.file(url = gtex_metadata_bw_paths[i], 
                  destfile = str_c("/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/", gtex_metadata$bigwig_file[i]), 
                  method = "wget", quiet = T)

  }
  
}

# Main ------------------------------------------------------------------------------------------------

##### Download bw to calculate the mean coverage from #####

# download all the sample bws for the GTEx project from recount, these will take up ~ 2TB and around 5 days to complete - should be deleted after use
# download_study(project = "SRP012682", type = "samples", download = T,
#                outdir = "/data/recount/SRP012682_GTEx/gtex_tissues_ERs_diff_cut_offs/")

# only run once - test if the method works to download urls
# download_hipp_cerebellum_fibroblast_bw_urls()

download_all_gtex_tissue_bw()
