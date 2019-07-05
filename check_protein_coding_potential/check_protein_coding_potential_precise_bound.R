# load the packages that will be used in the script
library(tidyverse) # data manipulation package in R 
library(Biostrings) # for manipulating DNA strings
library(UniProt.ws) # blasting your sequences 

# Set Working Directory ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

# ER details 
load(file = "results/annotate_ERs/ERs_optimised_cut_off_max_gap_all_tissues_w_annot_df.rda")

# details for OMIM genes - phenotypes to gene 
OMIM_PM_CS_filtered <- read_delim("results/download_tidy_OMIM_data/OMIM_PM_CS_filtered_2018_05_29.csv", delim = ",")

gtex_tissue_name_formatting <- read_delim("raw_data/gtex_tissue_name_formatting/OMIM_gtex_tissue_name_formatting.csv", delim = ",")

# Functions -------------------------------------------------------------------------------------------

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_txDb_from_gtf.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/get_DNA_sequence.R")

get_gtex_split_read_table_mean_cov_n_samples_df <- function(gtex_tissue_name_formatting){
  
  gtex_split_read_table_annotated_paths <- 
    list.files("/data/recount/GTEx_SRP012682/gtex_split_read_table_annotated_rda/", full.names = T)
  
  gtex_split_read_table_df <- 
    data_frame(gtex_split_read_table_annotated_paths = gtex_split_read_table_annotated_paths, 
               tissue = 
                 gtex_split_read_table_annotated_paths %>% 
                 str_replace("/.*/", "") %>% 
                 str_replace("_split_read_table_annotated.rda", "") %>% 
                 str_replace("brain", "") %>% 
                 str_replace_all("_", ""))
  
  gtex_mean_cov_df <- 
    data_frame(gtex_mean_cov_paths = 
                 list.files("/data/recount/GTEx_SRP012682/gtex_mean_coverage/by_tissue_smfrze_use_me/", full.names = T), 
               tissue = 
                 list.files("/data/recount/GTEx_SRP012682/gtex_mean_coverage/by_tissue_smfrze_use_me/", full.names = F) %>% 
                 str_replace("brain", "") %>% 
                 str_replace_all("_", ""))
  
  gtex_split_read_table_mean_cov_df <- 
    gtex_split_read_table_df %>% 
    inner_join(gtex_mean_cov_df) %>% 
    left_join(gtex_tissue_name_formatting, by = c("tissue" = "gtex_tissue_name_simplified"))
  
  return(gtex_split_read_table_mean_cov_df)
  
}

get_split_read_percent_samples_threshold <- function(ERs_w_annotation_all_tissues_intron_inter, gtex_tissue_name_formatting, percent_samples){
  
  ERs_w_annotation_all_tissues_intron_inter_w_n_samples_total <- 
    ERs_w_annotation_all_tissues_intron_inter %>% 
    left_join(gtex_tissue_name_formatting %>% dplyr::select(OMIM_gtex_name, n_all), by = c("tissue" = "OMIM_gtex_name"))
  
  for(i in 1:nrow(ERs_w_annotation_all_tissues_intron_inter_w_n_samples_total)){
    
    total_sample_size <- ERs_w_annotation_all_tissues_intron_inter_w_n_samples_total$n_all[i] %>% as.integer()
    p_annot_junc_count_samples <- ERs_w_annotation_all_tissues_intron_inter_w_n_samples_total$p_annot_junc_count_samples_split_read_annot[i]
    
    p_annot_junc_count_samples_unlist <- 
      p_annot_junc_count_samples %>% 
      str_split(";") %>% 
      unlist() %>% 
      as.integer()
    
    p_annot_propor_samples <- p_annot_junc_count_samples_unlist/total_sample_size
    
    ERs_w_annotation_all_tissues_intron_inter_w_n_samples_total$p_annot_junc_propor_samples[i] <- str_c(p_annot_propor_samples, collapse = ";")
    ERs_w_annotation_all_tissues_intron_inter_w_n_samples_total$max_p_annot_junc_propor_samples[i] <- max(p_annot_propor_samples)
    
  }
  
  return(ERs_w_annotation_all_tissues_intron_inter_w_n_samples_total)
  
}

collapse_ERs_precise_boundaries <- function(gr){
  
  # get non-redundant hits without self - leaving only groups that are 
  hits <- findOverlaps(gr, type = "equal", drop.redundant = T, drop.self = T)
  
  print("getting ERs with no other ER with a precise boundary in any other tissue...")
  
  # "tissue specific" ERs that do not need to be merged
  gr_no_hits <- gr[-c(queryHits(hits), subjectHits(hits)) %>% unique()]
  
  print(str_c(round((length(gr)-length(gr_no_hits))/length(gr), digits = 4) * 100, "% of ERs have a matching ER with a precise boundary in at least 1 other tissue"))
  
  print("getting groups of ERs to be collapsed across diff tissues...")
  
  uniq_groups <- 
    hits %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    filter(!queryHits %in% subjectHits(hits))
  
  uniq_groups_w_self <- 
    uniq_groups %>% 
    bind_rows(data_frame(queryHits = uniq_groups$queryHits %>% unique(), 
                         subjectHits = uniq_groups$queryHits %>% unique())) %>% 
    arrange(queryHits)
  
  # here we add the groups 
  gr_grouped <- gr[uniq_groups_w_self$subjectHits]
  elementMetadata(gr_grouped)["group"] <- uniq_groups_w_self$queryHits
  
  # concatonating tissues from the same group together
  gr_grouped_w_tissue_concat <- 
    elementMetadata(gr_grouped) %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    group_by(group) %>%
    mutate(tissue_concat = str_c(tissue, collapse = ", "))
  
  elementMetadata(gr_grouped)["tissue"] <- gr_grouped_w_tissue_concat$tissue_concat
  elementMetadata(gr_grouped)["group"] <- NULL
  
  print("collapsing ranges across tissues... ")
  
  gr_grouped_uniq <- gr_grouped %>% unique()
  
  stopifnot(length(gr_grouped_uniq) == (gr_grouped_w_tissue_concat$group %>% unique() %>% length()))
  
  gr_ranges_collapsed <- c(gr_no_hits, gr_grouped_uniq)
  
  return(gr_ranges_collapsed)
  
}

get_ER_boundaries <- function(tissue_to_filter, split_read_ids, gtex_split_read_table_mean_cov_df, ERs_w_annotation_all_tissues_intron_inter_gene_filtered){
  
  load(
    gtex_split_read_table_mean_cov_df %>% 
      filter(OMIM_gtex_name == tissue_to_filter) %>% 
      .[["gtex_split_read_table_annotated_paths"]]
  )
  
  split_reads <- 
    gtex_split_read_table_annotated_only_junc_coverage %>% 
    filter(junID %in% split_read_ids)
  
  if((split_reads$strand %>% unique() %>% length()) != 1){
    
    print("Strands on split reads do not match...")
    
    return("strand_mismatch")
    
  }
  
  split_reads_gr_to_check <- 
    GRanges(str_c("chr", split_reads$chr, ":", (split_reads$start - 1), "-", (split_reads$stop + 1)))
  
  if((findOverlaps(split_reads_gr_to_check, drop.self = T) %>% length()) != 0){
    
    print("Split reads overlap with each other")
    
    return("split_read_overlap")
    
  }
  
  split_reads_gr <- 
    GRanges(c(str_c("chr", split_reads$chr, ":", split_reads$start - 1, "-", split_reads$start - 1), 
              str_c("chr", split_reads$chr, ":", split_reads$stop + 1, "-", split_reads$stop + 1)),
            strand = unique(split_reads$strand))
  
  ER_gr <- 
    GRanges(str_c(ERs_w_annotation_all_tissues_intron_inter_gene_filtered$seqnames, ":", 
                  ERs_w_annotation_all_tissues_intron_inter_gene_filtered$start, "-", 
                  ERs_w_annotation_all_tissues_intron_inter_gene_filtered$end))
  
  split_reads_ends_ER <- split_reads_gr[queryHits(findOverlaps(split_reads_gr, ER_gr))]

  stopifnot(min(split_reads$stop) == (min(start(split_reads_ends_ER)) - 1))
  stopifnot(max(split_reads$start) == (max(start(split_reads_ends_ER)) + 1))
  
  ER_precise <- 
    GRanges(str_c(unique(seqnames(split_reads_ends_ER)), ":", 
                  min(start(split_reads_ends_ER)), "-", 
                  max(start(split_reads_ends_ER))), 
            strand = strand(split_reads_ends_ER) %>% unique())
  
  
  return(ER_precise)
  
}

# Main ------------------------------------------------------------------------------------------------

ensembl_grch38_v92_genes_txdb_exons <- 
  generate_txDb_from_gtf(gtf_gff3_path = "/data/references/ensembl/gtf_gff3/v92/Homo_sapiens.GRCh38.92.gtf", 
                         output_path = "/data/references/ensembl/txdb_sqlite/v92/ensembl_grch38_v92_txdb.sqlite",
                         seq_levels_to_keep = c(1:22, "X", "Y", "MT"), genome_build = "hg38") %>% 
  exons(c("exon_name", "gene_id"))

gtex_split_read_table_mean_cov_df <- get_gtex_split_read_table_mean_cov_n_samples_df(gtex_tissue_name_formatting)

ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific <-
  ERs_w_annotation_all_tissues %>% # this object contains all ERs from all tissues
  filter(width > 30, # filter only ERs above 30bps long, length required for blast significance
         !str_detect(tissue, "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast"), # remove sex specific and cell lines
         ensembl_grch38_v92_region_annot != "exon, intergenic, intron") # remove the category of ER that overlaps an exon, intron and intergenic region

rm(ERs_w_annotation_all_tissues)

ERs_w_annotation_all_tissues_intron_inter <- 
  ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific %>% 
  filter(ensembl_grch38_v92_region_annot %in% c("intron", "intergenic")) %>% # filter for only intronic and intergenic (novel) ERs
  filter(!is.na(uniq_genes_split_read_annot), !str_detect(uniq_genes_split_read_annot, ","), # filter for ERs that connect to only 1 gene
         (str_count(p_annot_junc_ids_split_read_annot, ";") == 1), # exactly 2 split reads (partially annotated, so they must connect to a known exon boundary) intersect the novel ER
         is.na(unannot_junc_ids_split_read_annot)) %>% # no unannotated split reads (where neither end falls onto a known exon boundary) connect to the ER 
  get_split_read_percent_samples_threshold(gtex_tissue_name_formatting, percent_samples) 

rm(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific)

ER_precise_w_seq_all_genes <- data_frame()

for(i in 1:nrow(ERs_w_annotation_all_tissues_intron_inter)){
  
  ensembl_id_to_filter <- ERs_w_annotation_all_tissues_intron_inter$uniq_genes_split_read_annot[i]
  tissue_to_filter <- ERs_w_annotation_all_tissues_intron_inter$tissue[i]
  ERs_w_annotation_all_tissues_intron_inter_gene_filtered <-  ERs_w_annotation_all_tissues_intron_inter[i,]
  
  print(str_c(i , " - getting ER protein sequence from ", ensembl_id_to_filter))
  
  split_read_ids <- 
    ERs_w_annotation_all_tissues_intron_inter$p_annot_junc_ids_split_read_annot[i] %>% 
    str_split(";") %>% 
    unlist()
  
  # check only two split reads connect to the ER
  stopifnot(length(split_read_ids) == 2)
  
  ER_precise <- get_ER_boundaries(tissue_to_filter, split_read_ids, gtex_split_read_table_mean_cov_df, 
                                  ERs_w_annotation_all_tissues_intron_inter_gene_filtered)
  
  if(class(ER_precise) != "GRanges"){
    
    ER_precise_w_seq_all_genes <- 
      bind_rows(ER_precise_w_seq_all_genes, 
        GRanges(str_c(ERs_w_annotation_all_tissues_intron_inter_gene_filtered$seqnames, ":", 
                      ERs_w_annotation_all_tissues_intron_inter_gene_filtered$start, "-", 
                      ERs_w_annotation_all_tissues_intron_inter_gene_filtered$end), 
                ensembl_id = ensembl_id_to_filter, 
                ER_tissue = tissue_to_filter, 
                split_read_ids = str_c(split_read_ids, collapse = ", "),
                DNA_seq = ER_precise, 
                orig_ER_chr = ERs_w_annotation_all_tissues_intron_inter_gene_filtered$seqnames, 
                orig_ER_start = ERs_w_annotation_all_tissues_intron_inter_gene_filtered$start, 
                orig_ER_end = ERs_w_annotation_all_tissues_intron_inter_gene_filtered$end) %>% as.data.frame() 
      )
    
  }else{
    
    ER_precise$ensembl_id <- ensembl_id_to_filter
    ER_precise$ER_tissue <- tissue_to_filter
    ER_precise$split_read_ids <- str_c(split_read_ids, collapse = ", ")
    
    ER_precise_w_seq <- get_DNA_sequence(GRCh38.dna, gr = ER_precise)
    
    ER_precise_w_seq$prot_seq_1 <- Biostrings::translate(DNAString(ER_precise_w_seq$DNA_seq)) %>% as.character()
    ER_precise_w_seq$prot_seq_2 <- Biostrings::translate(DNAString(ER_precise_w_seq$DNA_seq %>% str_replace("^.", ""))) %>% as.character()
    ER_precise_w_seq$prot_seq_3 <- Biostrings::translate(DNAString(ER_precise_w_seq$DNA_seq %>% str_replace("^..", ""))) %>% as.character()
    
    ER_precise_w_seq$orig_ER_chr = ERs_w_annotation_all_tissues_intron_inter_gene_filtered$seqnames
    ER_precise_w_seq$orig_ER_start = ERs_w_annotation_all_tissues_intron_inter_gene_filtered$start 
    ER_precise_w_seq$orig_ER_end = ERs_w_annotation_all_tissues_intron_inter_gene_filtered$end
    
    ER_precise_w_seq_all_genes <- bind_rows(ER_precise_w_seq_all_genes, ER_precise_w_seq %>% as.data.frame())
    
  }
  
}

stopifnot(all(str_count(ER_precise_w_seq_all_genes$split_read_ids[1], ",") == 1))

save(ER_precise_w_seq_all_genes, file = str_c("results/check_protein_coding_potential/ERs_precise_w_prot_seq.rda"))

# x <- list.files("results_tmp/protein_coding_potential/", full.names = T)
# 
# for(i in c(5200, 10400, 15600, 20800, 25466)){
#   
#   print(i)
#   
#   ensembl_id_to_filter <- ERs_w_annotation_all_tissues_intron_inter$uniq_genes_split_read_annot[i]
#   load(x[str_detect(x, ensembl_id_to_filter)])
#   
#   ER_precise_w_seq_all_genes_rbind <- bind_rows(ER_precise_w_seq_all_genes_rbind, ER_precise_w_seq_all_genes)
#   
# }





