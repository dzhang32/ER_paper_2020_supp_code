library(tidyverse) # data manipulation package in R 
library(Biostrings) # for manipulating DNA strings
library(ggpubr)
library(rtracklayer)

# Set Working Directory ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data ----------------------------------------------------------------------------------------------------------

load(file = "results/check_protein_coding_potential/ER_2_split_reads_precise_w_seq_tidy.rda")

load(file = "results/annotate_ERs/ERs_optimised_cut_off_max_gap_all_tissues_w_annot_df.rda")

gtex_tissue_name_formatting <- read_delim("raw_data/gtex_tissue_name_formatting/OMIM_gtex_tissue_name_formatting.csv", delim = ",")

Homo_sapiens.GRCh38.95 <- import("/data/references/ensembl/gtf_gff3/v95/Homo_sapiens.GRCh38.95.gtf")

# Functions ----------------------------------------------------------------------------------------------------------

##### First level #####

get_gtex_split_read_table_mean_cov_n_samples_df <- function(gtex_tissue_name_formatting){
  
  gtex_split_read_table_annotated_paths <- 
    list.files("/data/recount/GTEx_SRP012682/gtex_split_read_table_annotated_rda/", full.names = T)
  
  gtex_split_read_table_df <- 
    data_frame(gtex_split_read_table_annotated_paths = gtex_split_read_table_annotated_paths, 
               tissue = 
                 gtex_split_read_table_annotated_paths %>% 
                 str_replace(".*/", "") %>% 
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

get_ER_table_to_display <- function(ERs_w_annotation_all_tissues_intron_inter_1_split_read_1_gene){
  
  ERs_w_annotation_df_to_display <- 
    ERs_w_annotation_all_tissues_intron_inter_1_split_read_1_gene %>% 
    mutate(misannot_type = ifelse(!is.na(uniq_genes_split_read_annot), "split_read",
                                  ifelse(!is.na(overlap_any_gene_v92_name), "overlap", 
                                         ifelse(nearest_any_gene_v92_distance <= 10000, "within_10Kb", "none"))), 
           associated_gene = ifelse(misannot_type == "split_read", uniq_genes_split_read_annot, 
                                    ifelse(misannot_type == "overlap", overlap_any_gene_v92_name, 
                                           ifelse(misannot_type == "within_10Kb", nearest_any_gene_v92_name, NA))))
  
  ERs_w_annotation_df_to_display_w_split_read_data <- 
    ERs_w_annotation_df_to_display %>% 
    dplyr::select(ER_chr = seqnames, ER_start = start, ER_end = end, tissue, mean_coverage = value, ensembl_grch38_v92_region_annot, misannot_type, associated_gene,
                  mean_CDTS_percentile, mean_phast_cons_7)
  
  return(ERs_w_annotation_df_to_display_w_split_read_data)
  
}

##### Second level #####

get_split_read_info_for_ER_details_df <- function(ERs_w_annotation_df_to_display_tissue_filtered, tissue_n, gtex_split_read_table_annotated_only_junc_coverage){
  
  ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details <- 
    ERs_w_annotation_df_to_display_tissue_filtered %>% mutate(split_read_count_samp = as.integer(NA), 
                                                              split_read_propor_samp = as.numeric(NA),
                                                              split_read_starts = as.integer(NA), 
                                                              split_read_ends = as.integer(NA))
  
  gtex_split_read_table_annotated_only_junc_coverage_int_jun_id <- 
    gtex_split_read_table_annotated_only_junc_coverage %>% as_tibble() %>% mutate(junID = junID %>% as.integer()) %>% arrange(junID)
  
  for(j in 1:nrow(ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details)){
    
    p_annot_junc_ids_split_reads_int <- ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details$p_annot_junc_ids_split_read_annot[j] %>% 
      str_split(";") %>% unlist() %>% as.integer() %>% sort()
    
    if(length(p_annot_junc_ids_split_reads_int) == 0){
      
      next
      
    }else{
      
      gtex_split_read_table_annotated_only_junc_coverage_int_jun_id_filtered <- 
        gtex_split_read_table_annotated_only_junc_coverage_int_jun_id %>% 
        filter(junID %in% p_annot_junc_ids_split_reads_int)
      
      ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details$split_read_count_samp[j] <- 
        gtex_split_read_table_annotated_only_junc_coverage_int_jun_id_filtered$countsSamples %>% as.integer() %>% str_c(collapse = ";")
      
      ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details$split_read_propor_samp[j] <- 
        ((gtex_split_read_table_annotated_only_junc_coverage_int_jun_id_filtered$countsSamples %>% as.integer())/tissue_n) %>% round(digits = 3) %>% str_c(collapse = ";")
      
      ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details$split_read_starts[j] <- 
        gtex_split_read_table_annotated_only_junc_coverage_int_jun_id_filtered$start %>% str_c(collapse = ";")
      
      ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details$split_read_ends[j] <- 
        gtex_split_read_table_annotated_only_junc_coverage_int_jun_id_filtered$stop %>% str_c(collapse = ";")
      
    }
    
  }
  
  return(ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details)
  
}


# Main ---------------------------------------------------------------------------------------------------------------

gtex_split_read_table_mean_cov_df <- get_gtex_split_read_table_mean_cov_n_samples_df(gtex_tissue_name_formatting)

ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific <-
  ERs_w_annotation_all_tissues %>%
  filter(width > 3, !str_detect(tissue, "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast"),
         ensembl_grch38_v92_region_annot != "exon, intergenic, intron")

ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_tidy <- 
  ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific %>% 
  mutate(assoc_to_gene_type = ifelse(!is.na(uniq_genes_split_read_annot), "split_read",
                                ifelse(!is.na(overlap_any_gene_v92_name), "overlap", 
                                       ifelse(nearest_any_gene_v92_distance <= 10000, "within_10Kb", "none"))), 
         associated_gene = ifelse(assoc_to_gene_type == "split_read", uniq_genes_split_read_annot, 
                                  ifelse(assoc_to_gene_type == "overlap", overlap_any_gene_v92_name, 
                                         ifelse(assoc_to_gene_type == "within_10Kb", nearest_any_gene_v92_name, NA)))) %>% 
    dplyr::select(ER_chr = seqnames, ER_start = start, ER_end = end, tissue, 
                  ensembl_grch38_v92_region_annot, assoc_to_gene_type, associated_gene, 
                  mean_CDTS_percentile, mean_phast_cons_7)

ERs_w_annotation_all_tissues_intron_inter_1_split_read_1_gene <- 
  ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific %>% 
  filter(ensembl_grch38_v92_region_annot %in% c("intron", "intergenic")) %>% 
  filter(!is.na(uniq_genes_split_read_annot), !str_detect(uniq_genes_split_read_annot, ","), 
         (str_count(p_annot_junc_ids_split_read_annot, ";") == 0), 
         is.na(unannot_junc_ids_split_read_annot)) 

ERs_w_annotation_all_tissues_intron_inter_1_split_read_1_gene_tidy <- get_ER_table_to_display(ERs_w_annotation_all_tissues_intron_inter_1_split_read_1_gene)

ERs_w_annotation_all_tissues_intron_inter_2_split_read_1_gene_tidy <- 
  ER_precise_gr_w_conserv_constraint_tidy %>% 
  filter(!duplicated(ER_index)) %>% 
  mutate(misannot_type = "split_read") %>% 
  dplyr::select(ER_chr = seqnames, ER_start = start, ER_end = end, tissue = ER_tissue, ensembl_grch38_v92_region_annot, assoc_to_gene_type = misannot_type, associated_gene = ensembl_id,
                prot_coding_potential = protein_potential_any, 
                mean_CDTS_percentile, mean_phast_cons_7 = mean_phastCons7way, mean_phastCons20way = mean_phastCons20way)
  
Homo_sapiens.GRCh38.95_gene_name_key <- Homo_sapiens.GRCh38.95 %>% 
  elementMetadata() %>% 
  as.data.frame() %>% 
  dplyr::select(gene_id, gene_name) %>% 
  filter(!duplicated(gene_id))

ERs_prot_coding_potential_w_constraint_conserv <- 
  ERs_w_annotation_all_tissues_intron_inter_2_split_read_1_gene_tidy %>% 
  filter(prot_coding_potential == T) %>% 
  mutate(ER_chr_start_end = str_c(ER_chr, "_", ER_start, "_", ER_end)) %>% 
  group_by(ER_chr_start_end) %>% 
  mutate(all_tissue = str_c(tissue, collapse = ";")) %>% 
  filter(!duplicated(ER_chr_start_end)) %>% 
  left_join(Homo_sapiens.GRCh38.95_gene_name_key, by = c("associated_gene" = "gene_id")) %>% 
  dplyr::select(-ER_start, -ER_end, -ER_chr, -tissue) %>% 
  dplyr::select(ER_chr_start_end, prot_coding_potential, all_tissue, associated_gene_ens = associated_gene, associated_gene_symbol = gene_name, 
                everything())


# Save data ----------------------------------------------------------------------------------------------------------


write_delim(ERs_w_annotation_all_tissues_intron_inter_1_split_read_1_gene_tidy,  
            "results/export_ER_details/ERs_1_split_read_1_gene.csv", delim = ",")

write_delim(ERs_w_annotation_all_tissues_intron_inter_2_split_read_1_gene_tidy,  
            "results/export_ER_details/ERs_2_split_read_1_gene_high_confidence.csv", delim = ",")

write_delim(ERs_prot_coding_potential_w_constraint_conserv,  
            "results/export_ER_details/ERs_prot_coding_potential.csv", delim = ",")

write_delim(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_tidy,  
            "results/export_ER_details/all_ERs.csv", delim = ",")

