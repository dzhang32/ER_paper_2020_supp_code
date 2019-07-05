library(tidyverse)
library(stringr)
library(annotatER)
library(GenomicRanges)
library(SummarizedExperiment)
library(derfinder)
library(rtracklayer)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Functions -------------------------------------------------------------------------------------------

##### First level #####

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_txDb_from_gtf.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_genomic_state.R")

get_non_overlapping_exons <- function(ensembl_grch38_v87_TxDb){
  
  ensembl_grch38_v87_exons_gr <- 
    ensembl_grch38_v87_TxDb %>% exons(columns = c("EXONNAME", "GENEID", "TXNAME"))
  
  ensembl_grch38_v87_exons_gr_marked_overlapping <- 
    mark_overlapping_genes_gr(gr_1 = ensembl_grch38_v87_exons_gr, gr_2 = ensembl_grch38_v87_exons_gr, identical_gr = T, 
                              maxgap = -1L, minoverlap = 1L)
  
  ensembl_grch38_v87_exons_gr_non_overlapping <- 
    ensembl_grch38_v87_exons_gr_marked_overlapping[ensembl_grch38_v87_exons_gr_marked_overlapping$overlap_gr2 == F]
  
  return(ensembl_grch38_v87_exons_gr_non_overlapping)
  
}

get_ER_annot_details_per_cut_off <- function(ERs_cut_off_gr_no_schaffold, genomic_state_ensembl_grch38_v87){
  
  # annotate ERs and get number of each ER for each 
  ERs_cut_off_filtered_annotation <- 
    annotateRegions(regions = ERs_cut_off_gr_no_schaffold, 
                    genomicState = genomic_state_ensembl_grch38_v87$fullGenome, 
                    annotate = F)
  
  ERs_cut_off_filtered_annotation_w_region_annot <- 
    convert_annot_count_table_to_region_annot(ERs_cut_off_filtered_annotation[["countTable"]])
  
  ERs_cut_off_filtered_annotation_w_region_annot_summarised <- 
    ERs_cut_off_filtered_annotation_w_region_annot %>% 
    group_by(region_annot) %>% 
    summarise(n = n()) %>% 
    spread(key = "region_annot", value = "n") 
  
  colnames(ERs_cut_off_filtered_annotation_w_region_annot_summarised) <- 
    colnames(ERs_cut_off_filtered_annotation_w_region_annot_summarised) %>% str_c("num: ", .)
  
  return(ERs_cut_off_filtered_annotation_w_region_annot_summarised)
  
}

get_exon_delta_details_per_cut_off <- function(ERs_cut_off_gr_no_schaffold, ensembl_grch38_v87_exons_gr_non_overlapping){
  
  ERs_cut_off_gr_no_schaffold_overlapping_exon_hits <- 
    findOverlaps(query = ERs_cut_off_gr_no_schaffold, subject = ensembl_grch38_v87_exons_gr_non_overlapping) %>% 
    as.data.frame() %>% 
    as_tibble()
  
  # get those hits where multiple exons overlap 1 ER 
  num_exon_hits_per_ER_hit <- 
    ERs_cut_off_gr_no_schaffold_overlapping_exon_hits %>% 
    group_by(queryHits) %>% 
    summarise(n_distinct_exons = n_distinct(subjectHits))
  
  # get num ERs whereby multiple ERs overlap 1 exon
  num_ER_hits_per_exon_hit <- 
    ERs_cut_off_gr_no_schaffold_overlapping_exon_hits %>% 
    group_by(subjectHits) %>% 
    summarise(n_distinct_ERs = n_distinct(queryHits))
  
  exon_indexes_with_multiple_overlapping_ERs <- 
    num_ER_hits_per_exon_hit %>% 
    filter(n_distinct_ERs > 1) %>% 
    .[["subjectHits"]]
  
  multiple_ERs_overlapping_1_exon <- 
    ERs_cut_off_gr_no_schaffold_overlapping_exon_hits %>% 
    filter(subjectHits %in% exon_indexes_with_multiple_overlapping_ERs)
  
  # remove hits where multiple exons overlap 1 ER 
  ERs_overlapping_ab_1_exon_index <- 
    num_exon_hits_per_ER_hit %>% 
    filter(n_distinct_exons > 1) %>% 
    .[["queryHits"]]
  
  ERs_hits_not_overlapping_multiple_exons <- 
    ERs_cut_off_gr_no_schaffold_overlapping_exon_hits %>% 
    filter(!queryHits %in% ERs_overlapping_ab_1_exon_index)
  
  # merge ER hits with exon hits and calculate the delta exon
  ER_hits_ranges <- ERs_cut_off_gr_no_schaffold[ERs_hits_not_overlapping_multiple_exons$queryHits] %>% ranges() %>% as.data.frame()
  exon_hits_ranges <- ensembl_grch38_v87_exons_gr_non_overlapping[ERs_hits_not_overlapping_multiple_exons$subjectHits] %>% ranges() %>% as.data.frame()
  
  exon_delta_df <- 
    bind_cols(ER_hits_ranges, exon_hits_ranges) %>% 
    mutate(start_diff = start - start1, 
           end_diff = end - end1, 
           exon_delta = abs(start_diff) + abs(end_diff))
  
  # get details of delta exon including the number of ERs that overlap 1/multiple exons
  exon_delta_details <- 
    data_frame(exon_delta_sum = exon_delta_df$exon_delta %>% sum(), 
               exon_delta_mean = exon_delta_df$exon_delta %>% mean(), 
               exon_delta_sd = exon_delta_df$exon_delta %>% sd(), 
               exon_delta_median = exon_delta_df$exon_delta %>% median(),
               num_exon_delta_eq_0 = sum(exon_delta_df$exon_delta == 0), 
               num_ERs_overlapping_all_exons = nrow(num_exon_hits_per_ER_hit), 
               num_ERs_overlapping_ab_1_exon = length(ERs_overlapping_ab_1_exon_index), 
               num_ERs_for_ab_1_ER_overlapping_1_exon = nrow(multiple_ERs_overlapping_1_exon)) %>% 
    mutate(num_ERs_tested_for_exon_delta = num_ERs_overlapping_all_exons - num_ERs_overlapping_ab_1_exon)
  
  return(exon_delta_details)
  
}

merge_exon_delta_details_diff_tissues <- function(){
  
  exon_details_all_tissues_paths <- list.files("results/optimise_derfinder_cut_off/exon_delta_details_per_tissue/", full.names = T)
  
  exon_details_all_tissues <- data_frame()
  
  for(i in seq_along(exon_details_all_tissues_paths)){
    
    exon_details_one_tissue <- read_delim(exon_details_all_tissues_paths[i], delim = ",")
    
    exon_details_all_tissues <- 
      bind_rows(exon_details_all_tissues, exon_details_one_tissue)
    
  }
  
  return(exon_details_all_tissues)
  
}


##### Second level #####

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/mark_overlapping_genes_gr.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/convert_annot_count_table_to_region_annot.R")

# Main ------------------------------------------------------------------------------------------------

# get all non-overlapping exons from the ensembl TXDb
ensembl_grch38_v87_TxDb <- 
  generate_txDb_from_gtf(gtf_gff3_path = "/data/references/ensembl/gtf_gff3/v87/Homo_sapiens.GRCh38.87.gtf",
                         output_path = "/data/references/ensembl/txdb_sqlite/v87/ensembl_grch38_v87_txdb.sqlite",
                         seq_levels_to_keep = c(1:22, "X", "Y", "MT"), 
                         genome_build = "hg38")

ensembl_grch38_v87_exons_gr_non_overlapping <- 
  get_non_overlapping_exons(ensembl_grch38_v87_TxDb)

# get genomic state to annotate regions
# genomic_state_ensembl_grch38_v87 <- 
#   generate_genomic_state(ensembl_grch38_v87_TxDb, 
#                          output_path = "/data/references/ensembl/genomic_state/v87/ensembl_grch38_v87_genomic_state.rda", 
#                          chrs_to_keep = str_c("chr", c(1:22, "X", "Y", "M")))


tissue_varied_cut_off_maxgap_paths <- 
  list.files("/data/recount/GTEx_SRP012682/gtex_ERs_varying_maxgaps/", full.names = T, recursive = T)

tissue_varied_cut_off_maxgap_paths_df <- 
  data_frame(tissue_varied_cut_off_maxgap_paths = tissue_varied_cut_off_maxgap_paths, 
             tissue = tissue_varied_cut_off_maxgap_paths %>% 
               str_replace("/.*//", "") %>% 
               str_replace("/.*", ""))

# generate df to store all results 
ER_annotation_exon_delta_details_all_tissues <- data_frame()

# loop across brain regions then cutoffs and then grab details of ERs at each cutoff
for(i in 1:nrow(tissue_varied_cut_off_maxgap_paths_df)){
  
  tissue_varied_cut_off_maxgap_path <- tissue_varied_cut_off_maxgap_paths_df$tissue_varied_cut_off_maxgap_paths[i] 
  tissue <- tissue_varied_cut_off_maxgap_paths_df$tissue[i] 
  
  print(str_c("loading - ", i, " - ", tissue))
  
  load(tissue_varied_cut_off_maxgap_path)

  # generate exmpty dataframe to store ER annotation and exon delta details for each brain region 
  ER_annotation_exon_delta_details_per_tissue <- data_frame()
  
  cut_offs <- names(ERs_tissue_all_cut_offs_all_maxgaps)
  
  for(j in seq_along(cut_offs)){
    
    cut_off_to_filter <- cut_offs[j]
    
    print(str_c(j, " - cutoff: ", cut_off_to_filter))
    
    # filter to specific GR for one cutoff
    ERs_one_cut_off_all_maxgaps_gr <- 
      ERs_tissue_all_cut_offs_all_maxgaps[[as.character(cut_off_to_filter)]] 
    
    ER_annotation_exon_delta_details_per_cut_off <- data_frame()
    
    for(k in 1:length(ERs_one_cut_off_all_maxgaps_gr)){
      
      maxgap <- 
      names(ERs_one_cut_off_all_maxgaps_gr[k]) %>% 
        str_replace(".*:", "")
      
      print(str_c(k, " - maxgap: ", maxgap))
      
      ERs_one_cut_off_one_maxgap_gr <- ERs_one_cut_off_all_maxgaps_gr[[k]]
      
      # print("annotating regions by feature...")
      # 
      # ERs_cut_off_filtered_annotation_w_region_annot_summarised <-
      #   get_ER_annot_details_per_cut_off(ERs_one_cut_off_one_maxgap_gr, genomic_state_ensembl_grch38_v87)

      print("getting exon delta details...")
      
      exon_delta_details <- 
        get_exon_delta_details_per_cut_off(ERs_one_cut_off_one_maxgap_gr, ensembl_grch38_v87_exons_gr_non_overlapping) %>% 
        mutate(total_ERs = length(ERs_one_cut_off_one_maxgap_gr), 
               tissue = tissue,
               cut_off = cut_off_to_filter, 
               maxgap = maxgap)
      
      ER_annotation_exon_delta_details_per_cut_off <- 
        ER_annotation_exon_delta_details_per_cut_off %>% 
        bind_rows(exon_delta_details)
      
    }
    
    ER_annotation_exon_delta_details_per_tissue <- 
      ER_annotation_exon_delta_details_per_tissue %>% 
      bind_rows(ER_annotation_exon_delta_details_per_cut_off)
    
  }
  
  write_delim(ER_annotation_exon_delta_details_per_tissue,
              str_c("results/optimise_derfinder_cut_off/exon_delta_details_per_tissue/ER_exon_delta_details_", tissue, "_cutoff_0.2_to_10_by_0.2_maxgap_0_to_100.csv"), delim = ",")
  
  rm(ERs_tissue_all_cut_offs_all_maxgaps)
  rm(ERs_one_cut_off_all_maxgaps_gr)
  rm(ERs_one_cut_off_one_maxgap_gr)
  rm(ER_annotation_exon_delta_details_per_tissue)
  
}

exon_details_all_tissues <- merge_exon_delta_details_diff_tissues()

write_delim(exon_details_all_tissues, 
            "results/optimise_derfinder_cut_off/exon_delta_details_all_tissues_cut_off_1_10_0.2_maxgap_0_100_10.csv", delim = ",")


