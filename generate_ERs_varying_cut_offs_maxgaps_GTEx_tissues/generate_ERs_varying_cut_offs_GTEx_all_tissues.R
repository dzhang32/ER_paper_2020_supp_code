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

get_bws_no_brain_fib <- function(){
  
  gtex_metadata <-  all_metadata('gtex')
  
  gtex_metadata_all_tissues_no_brain_fib <- 
    gtex_metadata %>% 
    as.data.frame() %>% 
    filter(!str_detect(smtsd, "Brain"), smtsd != "Cells - Transformed fibroblasts")
  
  bigWig_paths <- 
    str_c("/data/recount/GTEx_SRP012682/gtex_bigWigs/all_gtex_tissues_raw_bigWigs/", gtex_metadata_all_tissues_no_brain_fib$bigwig_file)
  
  mapped_read_count <- 
    gtex_metadata_all_tissues_no_brain_fib$mapped_read_count
  
  stopifnot(length(bigWig_paths) == nrow(gtex_metadata_all_tissues_no_brain_fib))
  
  gtex_metadata_all_tissues_no_brain_fib_df <- 
    data_frame(bigWig_paths = bigWig_paths, 
               tissue = bigWig_paths %>% 
                 str_replace_all("/.*_", "") %>% 
                 str_replace("\\.bw", "") %>% 
                 str_replace_all("\\.", "_"), 
               mapped_read_count = mapped_read_count, 
               auc = gtex_metadata_all_tissues_no_brain_fib$auc,
               smafrze = gtex_metadata_all_tissues_no_brain_fib$smafrze) %>% 
    arrange(tissue)
  
  return(gtex_metadata_all_tissues_no_brain_fib_df)
  
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

##### Generate the ERs at varying cut-offs for each tissue ##### 

gtex_metadata_all_tissues_no_brain_fib_df <- get_bws_no_brain_fib()

gtex_tissues <- 
  gtex_metadata_all_tissues_no_brain_fib_df$tissue %>% unique() 
# gtex_tissues <- gtex_tissues[!str_detect(gtex_tissues, "brain_hippocampus")]

# get chr lengths to input into the loadCoverage function()
# chrominfo <- fetchExtendedChromInfoFromUCSC("hg38")
# chrominfo %>% 
#   as_tibble() %>% 
#   write_delim("/home/dzhang/data/UCSC_chr_info/UCSC_chr_info_hg38.csv", delim = ",")
chrominfo <- read_delim("/home/dzhang/data/UCSC_chr_info/UCSC_chr_info_hg38.csv", delim = ",")
chrominfo_no_scaffold <- 
  chrominfo %>% 
  filter(UCSC_seqlevel %in% str_c("chr", c(1:22, "X", "Y", "M")))

# get cut offs to generate ERs 
cut_offs <- seq(1, 10, 0.2)

for(i in seq_along(gtex_tissues)){
  
  gtex_tissue_to_filter <- gtex_tissues[i]
  
  print(str_c(Sys.time(), " - ", i, " - ", gtex_tissue_to_filter))
  
  gtex_metadata_all_tissues_no_brain_fib_df_tissue_filtered <-
    gtex_metadata_all_tissues_no_brain_fib_df %>%
    filter(tissue == gtex_tissue_to_filter, 
           smafrze == "USE ME") 
  
  # this causes an error in generating the mean coverage if no samples are left
  if(nrow(gtex_metadata_all_tissues_no_brain_fib_df_tissue_filtered) == 0){
    
    print(str_c("Warning: no samples found left after smafrze filtering for ", gtex_tissue_to_filter))
    
    next
    
  }
  
  for(j in 1:nrow(chrominfo_no_scaffold)){
    
    chr_to_load <- chrominfo_no_scaffold[["UCSC_seqlevel"]][j]
    chr_length <- chrominfo_no_scaffold[["UCSC_seqlength"]][j]
    
    print(str_c(Sys.time(), " - loading chromosome: ", chr_to_load))
    
    # normalise regions to 40 mill reads - could also do in two steps, with load coverage then filter
    tissue_coverage_w_mean_normalised <-
      loadCoverage(files = gtex_metadata_all_tissues_no_brain_fib_df_tissue_filtered$bigWig_paths, 
                   totalMapped = gtex_metadata_all_tissues_no_brain_fib_df_tissue_filtered$auc, # normalise by auc here as for bws, more accurate since Rail-RNA clips reads
                   targetSize = (40000000 * 100), # target 40 million coverage with 100 bp length reads
                   chr = chr_to_load, 
                   chrlen = chr_length, 
                   inputType = "BigWig",
                   returnMean = T,
                   returnCoverage = F,
                   verbose = F, 
                   cutoff = NULL) # setting cutoff as null here and instead will be applied downstream in findRegions()
    
    results_path_mean_cov_tissue <- 
      make_results_dir(results_path = "/data/recount/GTEx_SRP012682/gtex_mean_coverage/by_tissue_smfrze_use_me/", folder_name = gtex_tissue_to_filter)
    
    save(tissue_coverage_w_mean_normalised, 
         file = str_c(results_path_mean_cov_tissue, "/gtex_", gtex_tissue_to_filter, "_", chr_to_load, "_mean_cov.rda"))
    
    results_path_tissue <- 
      make_results_dir(results_path = "/data/recount/GTEx_SRP012682/gtex_ERs_varying_cut_offs/by_chr/", folder_name = gtex_tissue_to_filter)
    
    for(k in seq_along(cut_offs)){
      
      cut_off_to_filter <- cut_offs[k]
      
      results_path_tissue_cut_off <- make_results_dir(results_path = results_path_tissue, folder_name = str_c("cutoff_", cut_off_to_filter))
      
      # get ERs for a particular cutoff
      # position here is usually returned by loadCoverage/filterData and is a logical vector that indicates which positions pass the cutoff,
      # however as cutoff left as NULL before, this is just a logical vector of length = T
      ERs_tissue_varying_cut_off <- 
        findRegions(position = Rle(TRUE, length(tissue_coverage_w_mean_normalised$meanCoverage)), 
                    fstats = tissue_coverage_w_mean_normalised$meanCoverage, 
                    chr = chr_to_load, 
                    cutoff = cut_off_to_filter, 
                    maxRegionGap = 0L,
                    maxClusterGap = chr_length, # setting this to chr length to reduce run time as we dont use the clusters and checked this doesn't change ERs
                    verbose = T)
      
      names(ERs_tissue_varying_cut_off) <- rep(chr_to_load, length(ERs_tissue_varying_cut_off))
      
      save(ERs_tissue_varying_cut_off, 
           file = str_c(results_path_tissue_cut_off, "/", gtex_tissue_to_filter, "_ERs_", chr_to_load, "_cutoff", cut_off_to_filter, "_regiongap0.rda"))
      
    }
    
    # to save memory after each run
    rm(tissue_coverage_w_mean_normalised)
    
  }
  
}


