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

##### First level #####

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_txDb_from_gtf.R")

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

get_distance_between_adjacent_ranges <- function(gr){
  
  # need to sort for this distance function to work between closest ranges
  gr <- gr %>% sort()
  
  dis_between_adjacent_ERs <- distance(gr[-length(gr)], gr[2:length(gr)]) 
  
  dis_between_adjacent_ERs_no_NA <- dis_between_adjacent_ERs[-is.na(dis_between_adjacent_ERs)]
  
  return(dis_between_adjacent_ERs)
  
}

get_dis_between_ERs_overlapping_1_exon <- function(gr, ensembl_grch38_v87_exons_gr_non_overlapping){
  
  # find the indexes for ERs where more than 1 ER overlaps 1 exon 
  multiple_ERs_overlapping_1_exon <- 
  findOverlaps(gr, ensembl_grch38_v87_exons_gr_non_overlapping) %>% 
    as.data.frame() %>% 
    group_by(subjectHits) %>% 
    mutate(n_dis_ERs = n_distinct(queryHits)) %>% 
    filter(n_dis_ERs > 1)
  
  # get the matching adjacent ERs indexes for ERs overlapping 1 exon 
  multiple_ERs_overlapping_1_exon_ERs_to_get_dis_between <- 
    multiple_ERs_overlapping_1_exon %>% 
    group_by(subjectHits) %>% 
    mutate(queryHits_to_get_dis = lead(queryHits)) %>% 
    filter(!is.na(queryHits_to_get_dis))
  
  # find distance between this adjacent ERs overlapping 1 exon 
  dis_between_ERs_overlapping_1_exon <- 
    distance(gr[multiple_ERs_overlapping_1_exon_ERs_to_get_dis_between$queryHits], 
             gr[multiple_ERs_overlapping_1_exon_ERs_to_get_dis_between$queryHits_to_get_dis])
  
  return(dis_between_ERs_overlapping_1_exon)
  
}

##### Second level #####

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/mark_overlapping_genes_gr.R")

# Main ------------------------------------------------------------------------------------------------

# get all non-overlapping exons from the ensembl TXDb
ensembl_grch38_v87_TxDb <- 
  generate_txDb_from_gtf(gtf_gff3_path = "/data/references/ensembl/gtf_gff3/v87/Homo_sapiens.GRCh38.87.gtf",
                         output_path = "/data/references/ensembl/txdb_sqlite/v87/ensembl_grch38_v87_txdb.sqlite",
                         seq_levels_to_keep = c(1:22, "X", "Y", "MT"), 
                         genome_build = "hg38")

ensembl_grch38_v87_exons_gr_non_overlapping <- 
  get_non_overlapping_exons(ensembl_grch38_v87_TxDb)

# generate df with ERs varying cut off paths, split by Leo or self-generated
# leo's contains 11 brain regions (all minus the crbl)
tissue_varied_cut_off_paths_leo <- list.files("/data/recount/", recursive = T, full.names = T, pattern = "region_cuts.Rdata")
# leo's contains remaining 42 tissues (- Cells - Leukemia cell line CML - no samples passing the smarze)
tissue_varied_cut_off_paths_gtex_tissues <- list.files("/data/recount/GTEx_SRP012682/gtex_ERs_varying_cut_offs/all_cutoffs_all_chr_merged/", 
                                                       recursive = T, full.names = T, pattern = ".rda")
tissue_varied_cut_off_paths_gtex_tissues <- tissue_varied_cut_off_paths_gtex_tissues[-str_detect(tissue_varied_cut_off_paths_gtex_tissues, "hippocampus")]

tissue_varied_cut_off_paths_df <- 
  data_frame(tissue_varied_cut_off_paths = c(tissue_varied_cut_off_paths_leo, tissue_varied_cut_off_paths_gtex_tissues), 
             leo_or_dz = c(rep("leo", length(tissue_varied_cut_off_paths_leo)), rep("dz", length(tissue_varied_cut_off_paths_gtex_tissues))), 
             tissue = tissue_varied_cut_off_paths %>% 
               str_replace("/region_cuts.Rdata", "") %>% 
               str_replace("_ERs_all_chrs_cutoff.*.rda", "") %>% 
               str_replace("/.*/", ""))

# will take a long time to run so filtering to a test group of tissues 
tissue_varied_cut_off_paths_df <- 
  tissue_varied_cut_off_paths_df %>% filter(tissue %in% c("amygdala", "hippocampus", "brain_cerebellum", "whole_blood"))

dis_between_ERs_all_tissues <- data_frame()

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
  
  cut_offs <- names(ERs_tissue_all_cut_offs_all_chrs)
  
  dis_between_ERs_per_tissue <- data_frame()
  
  for(j in seq_along(cut_offs)){
    
    cut_off_to_filter <- cut_offs[j]
    
    print(str_c(Sys.time(), " - loading - ", j, " - ", cut_off_to_filter))
    
    ERs_tissue_all_cut_offs_all_chrs_no_scaffold <- 
      ERs_tissue_all_cut_offs_all_chrs[[cut_off_to_filter]] %>% 
      keepSeqlevels(str_c("chr", c(1:22, "X", "Y", "M")), pruning.mode = "coarse") %>% 
      sort()
    
    # get df for distance between all adjacent ERs and those overlapping 1 exon per cut off
    dis_between_ERs_per_cut_off <- 
      bind_rows(data_frame(dis_between_ERs = get_dis_between_ERs_overlapping_1_exon(ERs_tissue_all_cut_offs_all_chrs_no_scaffold, ensembl_grch38_v87_exons_gr_non_overlapping), 
                           type = "multiple_ERs_overlapping_1_exon"), 
                data_frame(dis_between_ERs = get_distance_between_adjacent_ranges(gr = ERs_tissue_all_cut_offs_all_chrs_no_scaffold), 
                           type = "all")) %>% 
      mutate(cut_off = cut_off_to_filter, 
             tissue = tissue)
    
    dis_between_ERs_per_tissue <- 
      bind_rows(dis_between_ERs_per_tissue, dis_between_ERs_per_cut_off)
    
  }
  
  dis_between_ERs_all_tissues <- 
    bind_rows(dis_between_ERs_all_tissues, dis_between_ERs_per_tissue) 
    
}

write_delim(dis_between_ERs_all_tissues, "results/optimise_derfinder_cut_off/dis_between_ERs_amyg_hipp_crbl_blood.csv", delim = ",")


