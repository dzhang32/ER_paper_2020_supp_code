library(tidyverse)
library(stringr)
library(recount)
library(derfinder)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Functions -------------------------------------------------------------------------------------------

# SRP051844

get_meta_data_SRP051844 <- function(){
  
  all_metadata_sra <- all_metadata("sra")
  
  all_metadata_sra_SRP051844_controls <- 
    all_metadata_sra %>% 
    as_data_frame() %>% 
    filter(project == "SRP051844") %>% 
    filter(str_detect(characteristics, "Neurologically normal"))

  return(all_metadata_sra_SRP051844_controls)
  
}

download_bws_SRP051844 <- function(all_metadata_sra_SRP051844_controls){
  
  all_metadata_sra_SRP051844_controls_bigwig_paths <-
    all_metadata_sra_SRP051844_controls$bigwig_file %>%
    str_c("http://duffel.rail.bio/recount/SRP051844/bw/", .)
  
  for(i in seq_along(all_metadata_sra_SRP051844_controls_bigwig_paths)){
    
    print(str_c(i, " - ", all_metadata_sra_SRP051844_controls$bigwig_file[i]))
    
    download.file(url = all_metadata_sra_SRP051844_controls_bigwig_paths[i],
                  destfile = str_c("/data/recount/SRP051844/sample_big_wigs/", all_metadata_sra_SRP051844_controls$bigwig_file[i]),
                  method = "wget", quiet = T)
    
  }
  
}

generate_mean_cov_per_chr_SRP051844 <- function(all_metadata_sra_SRP051844_controls){
  
  # removed one sample with auc below the 40 million reads
  all_metadata_sra_SRP051844_controls_auc_ab_40_mill <- 
  all_metadata_sra_SRP051844_controls %>% filter(auc >= (40000000 * 100))
  
  all_metadata_sra_SRP051844_controls_bigwig_paths_mrserver_paths <- 
    all_metadata_sra_SRP051844_controls_auc_ab_40_mill$bigwig_file %>%
    str_c("/data/recount/SRP051844/sample_big_wigs/", .)
  
  chrominfo <- read_delim("/home/dzhang/data/UCSC_chr_info/UCSC_chr_info_hg38.csv", delim = ",")
  chrominfo_no_scaffold <- 
    chrominfo %>% 
    filter(UCSC_seqlevel %in% str_c("chr", c(1:22, "X", "Y", "M")))

  for(i in 1:nrow(chrominfo_no_scaffold)){
    
    chr_to_load <- chrominfo_no_scaffold[["UCSC_seqlevel"]][i]
    chr_length <- chrominfo_no_scaffold[["UCSC_seqlength"]][i]
    
    print(str_c(Sys.time(), " - loading chromosome: ", chr_to_load))
    
    SRP051844_control_mean_cov_normalised <-
      loadCoverage(files = all_metadata_sra_SRP051844_controls_bigwig_paths_mrserver_paths, 
                   totalMapped = all_metadata_sra_SRP051844_controls_auc_ab_40_mill$auc, # normalise by auc here as for bws, more accurate since Rail-RNA clips reads
                   targetSize = (40000000 * 100), # target 40 million coverage with 100 bp length reads
                   chr = chr_to_load, 
                   chrlen = chr_length, 
                   inputType = "BigWig",
                   returnMean = T,
                   returnCoverage = F,
                   verbose = F, 
                   cutoff = NULL) # setting cutoff as null here and instead will be applied downstream in findRegions()
    
    save(SRP051844_control_mean_cov_normalised, 
         file = str_c("/data/recount/SRP051844/mean_cov/controls_by_chr_normalised_40mill/SRP051844_control_mean_cov_40mill_", chr_to_load, ".rda"))
    
  }
  
}

# SRP058181

get_meta_data_SRP058181 <- function(){
  
  all_metadata_sra <- all_metadata("sra")
  
  all_metadata_sra_SRP058181_controls <- 
    all_metadata_sra %>% 
    as_data_frame() %>% 
    filter(project == "SRP058181") %>% 
    filter(str_detect(title, "C"))
  
  return(all_metadata_sra_SRP058181_controls)
  
}

download_bws_SRP058181 <- function(all_metadata_sra_SRP058181_controls){
  
  all_metadata_sra_SRP058181_controls_bigwig_paths <-
    all_metadata_sra_SRP058181_controls$bigwig_file %>%
    str_c("http://duffel.rail.bio/recount/SRP058181/bw/", .)
  
  for(i in seq_along(all_metadata_sra_SRP058181_controls_bigwig_paths)){
    
    print(str_c(i, " - ", all_metadata_sra_SRP058181_controls$bigwig_file[i]))
    
    download.file(url = all_metadata_sra_SRP058181_controls_bigwig_paths[i],
                  destfile = str_c("/data/recount/SRP058181/sample_big_wigs/", all_metadata_sra_SRP058181_controls$bigwig_file[i]),
                  method = "wget", quiet = T)
    
  }
  
}

generate_mean_cov_per_chr_SRP058181 <- function(all_metadata_sra_SRP058181_controls){
  
  # removed one sample with auc below the 40 million reads
  all_metadata_sra_SRP058181_controls_auc_ab_40_mill <- 
    all_metadata_sra_SRP058181_controls %>% filter(auc >= (40000000 * 100))
  
  all_metadata_sra_SRP058181_controls_bigwig_paths_mrserver_paths <- 
    all_metadata_sra_SRP058181_controls_auc_ab_40_mill$bigwig_file %>%
    str_c("/data/recount/SRP058181/sample_big_wigs/", .)
  
  chrominfo <- read_delim("/home/dzhang/data/UCSC_chr_info/UCSC_chr_info_hg38.csv", delim = ",")
  chrominfo_no_scaffold_2 <- 
    chrominfo %>% 
    filter(UCSC_seqlevel %in% str_c("chr", c(1:22, "X", "Y", "M")))
  
  for(i in 1:nrow(chrominfo_no_scaffold_2)){
    
    chr_to_load <- chrominfo_no_scaffold_2[["UCSC_seqlevel"]][i]
    chr_length <- chrominfo_no_scaffold_2[["UCSC_seqlength"]][i]
    
    print(str_c(Sys.time(), " - loading chromosome: ", chr_to_load))
    
    SRP058181_control_mean_cov_normalised <-
      loadCoverage(files = all_metadata_sra_SRP058181_controls_bigwig_paths_mrserver_paths, 
                   totalMapped = all_metadata_sra_SRP058181_controls_auc_ab_40_mill$auc, # normalise by auc here as for bws, more accurate since Rail-RNA clips reads
                   targetSize = (40000000 * 100), # target 40 million coverage with 100 bp length reads
                   chr = chr_to_load, 
                   chrlen = chr_length, 
                   inputType = "BigWig",
                   returnMean = T,
                   returnCoverage = F,
                   verbose = F, 
                   cutoff = NULL) # setting cutoff as null here and instead will be applied downstream in findRegions()
    
    save(SRP058181_control_mean_cov_normalised, 
         file = str_c("/data/recount/SRP058181/mean_cov/controls_by_chr_normalised_40mill/SRP058181_control_mean_cov_40mill_", chr_to_load, ".rda"))
    
  }
  
}

# Main ------------------------------------------------------------------------------------------------

all_metadata_sra_SRP051844_controls <- get_meta_data_SRP051844()
download_bws_SRP051844(all_metadata_sra_SRP051844_controls)
generate_mean_cov_per_chr_SRP051844(all_metadata_sra_SRP051844_controls)

all_metadata_sra_SRP058181_controls <- get_meta_data_SRP058181()
download_bws_SRP058181(all_metadata_sra_SRP058181_controls)
generate_mean_cov_per_chr_SRP058181(all_metadata_sra_SRP058181_controls)


