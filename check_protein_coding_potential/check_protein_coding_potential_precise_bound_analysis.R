library(tidyverse) # data manipulation package in R 
library(Biostrings) # for manipulating DNA strings
library(ggpubr)
library(rtracklayer)

# Set Working Directory ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data ----------------------------------------------------------------------------------------------------------

load(file = str_c("results/check_protein_coding_potential/ERs_precise_w_prot_seq.rda"))

OMIM_gtex_tissue_name_formatting <- read_delim("raw_data/gtex_tissue_name_formatting/OMIM_gtex_tissue_name_formatting.csv", delim = ",")

gene_reannotated_properties <- read_delim("results_tmp/gene_reannotated_properties.txt", delim = "\t")

Homo_sapiens.GRCh38.92 <- import("/data/references/ensembl/gtf_gff3/v92/Homo_sapiens.GRCh38.92.gtf")

load(file = "results/annotate_ERs/ERs_optimised_cut_off_max_gap_all_tissues_w_annot_df.rda")

# Functions ------------------------------------------------------------------------------------------------

##### First level #####

get_conserv_constraint_prot_coding_potential_ERs <- function(ER_precise_w_seq_all_genes, ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific){
  
  ERs_intron_inter <- 
  ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific %>% 
    filter(ensembl_grch38_v92_region_annot %in% c("intron", "intergenic")) %>% 
    mutate(chr_start_end = str_c(seqnames, "_", start, "_", end)) %>% 
    filter(!duplicated(chr_start_end)) %>% 
    dplyr::select(chr_start_end, ensembl_grch38_v92_region_annot)
  
  ER_precise_w_seq_all_genes_tidy <- 
    ER_precise_w_seq_all_genes %>% 
    mutate(chr_start_end = str_c(seqnames, "_", start, "_", end), 
           orig_chr_start_end = str_c(orig_ER_chr, "_", orig_ER_start, "_", orig_ER_end),
           chr_start_end_tissue = str_c(seqnames, "_", start, "_", end, "_", ER_tissue), 
           ER_index = row_number()) %>% 
    left_join(ERs_intron_inter, by = c("orig_chr_start_end" = "chr_start_end")) %>% 
    filter(!is.na(prot_seq_1)) %>% 
    gather(key = "frame", value = "prot_seq", contains("prot_seq")) %>% 
    mutate(frame = str_replace(frame, "prot_seq_", ""), 
           prot_potential = !str_detect(prot_seq, "\\*")) %>% 
    group_by(ER_index) %>% 
    mutate(protein_potential_any = any(prot_potential)) %>% 
    ungroup()
  
  ER_precise_gr <- 
    GRanges((str_c(ER_precise_w_seq_all_genes_tidy$seqnames, ":", 
                 ER_precise_w_seq_all_genes_tidy$start, "-", 
                 ER_precise_w_seq_all_genes_tidy$end) %>% unique()))
  
  ER_precise_gr_w_conserv <- get_conservation_score_for_regions_bw(bw_path = "/data/conservation/phastCons/hg38.phastCons7way.bw", gr = ER_precise_gr)
  ER_precise_gr_w_conserv <- get_conservation_score_for_regions_bw(bw_path = "/data/conservation/phastCons/hg38.phastCons20way.bw", gr = ER_precise_gr_w_conserv)
  ER_precise_gr_w_conserv_constraint <- get_constraint_score_for_regions(region_details_gr = ER_precise_gr_w_conserv, constraint_gr = CDTS_percentile_N7794_unrelated_all_chrs_gr)
  
  ER_precise_gr_w_conserv_constraint_tidy <- 
    ER_precise_gr_w_conserv_constraint %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(chr_start_end = str_c(seqnames, "_", start, "_", end)) %>% 
    dplyr::select(-seqnames, -start, -end, -width, -strand) %>% 
    left_join(x = ER_precise_w_seq_all_genes_tidy, y = ., by = "chr_start_end")
  
  return(ER_precise_gr_w_conserv_constraint_tidy)
  
}

plot_total_Kb_potential_prot_coding_ERs <- function(ER_precise_w_seq_all_genes_tidy){
  
  x <- ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific %>% 
    filter(connect_to_OMIM_gene_split_read_annot == "ENSG00000131095", ensembl_grch38_v92_region_annot %in% c("intron", "intergenic"))
  
  ER_precise_w_seq_all_genes_tidy_summarised_per_ER_to_plot <- 
    ER_precise_w_seq_all_genes_tidy %>% 
    filter(!duplicated(chr_start_end)) %>% 
    group_by(protein_potential_any) %>% 
    summarise(total_Kb = sum(width)/1000, 
              n_ER = n()) %>% 
    mutate(protein_potential_any = ifelse(protein_potential_any, "Yes", "No") %>% factor() %>% fct_relevel(c("Yes", "No")))
  
  prot_coding_potential_plot <- 
    ggplot(ER_precise_w_seq_all_genes_tidy_summarised_per_ER_to_plot, aes(x = prot_potential, y = total_Kb)) +
    geom_col(aes(fill = prot_potential), colour = "black") +
    scale_x_discrete(name = "Protein coding potential") +
    scale_y_continuous(name = "Total Kb") +
    scale_fill_manual(name = "Protein coding potential", values = c(get_palette("npg", 3)[3], "grey"), guide = F) + 
    theme_pubr() 
  
  return(prot_coding_potential_plot)
    
  
}

plot_num_tissues_vs_ER <- function(ER_precise_w_seq_all_genes){
  
  ER_precise_w_seq_all_genes_tidy <- 
    ER_precise_w_seq_all_genes %>% 
    filter(!is.na(prot_seq_1), 
           !(ER_tissue %in% c("brain_cerebellum", "cortex"))) %>% 
    mutate(ER_index = row_number(), 
           chr_start_end = str_c(seqnames, "_", start, "_", end), 
           brain_derived = ER_tissue %in% (OMIM_gtex_tissue_name_formatting %>% 
                                          filter(gtex_tissue_group == "brain") %>% 
                                            .[["OMIM_gtex_name"]]))

  ER_precise_w_seq_all_genes_no_dup <- 
    ER_precise_w_seq_all_genes %>% 
    mutate(chr_start_end = str_c(seqnames, "_", start, "_", end)) %>% 
    filter(!is.na(prot_seq_1))
  
  ER_precise_w_seq_all_genes_tissue_specificity <- 
    ER_precise_w_seq_all_genes_tidy %>% 
    group_by(chr_start_end) %>% 
    summarise(n_dis_tissues = n_distinct(ER_tissue), 
              brain_derived_any = any(brain_derived))
  
  hist_tissue_specificity <- 
    ggplot(ER_precise_w_seq_all_genes_tissue_specificity, aes(x = n_dis_tissues)) +
    geom_histogram(fill = "grey", colour = "black", bins = 41) + 
    scale_y_continuous(name = "Count") +
    scale_x_continuous(name = "Number of GTEx tissues", breaks = seq(0, 40, 5)) + 
    theme_pubr()
  
  cum_propor_tissue_specificity <- 
  ggplot(ER_precise_w_seq_all_genes_tissue_specificity, aes(x = n_dis_tissues)) +
    stat_ecdf(geom = "step", n = 41) +
    scale_y_continuous(name = "Cumulative proportion") +
    scale_x_continuous(name = "Number of GTEx tissues", breaks = seq(0, 40, 5)) + 
    theme_pubr()
  
  hist_tissue_specificity_brain <- 
    ggplot(ER_precise_w_seq_all_genes_tissue_specificity %>% filter(brain_derived_any == T) , aes(x = n_dis_tissues)) +
    geom_histogram(fill = "grey", colour = "black", bins = 41) + 
    scale_y_continuous(name = "Count") +
    scale_x_continuous(name = "Number of GTEx tissues (CNS only)", breaks = seq(0, 40, 5)) + 
    theme_pubr()
  
  cum_propor_tissue_specificity_brain <- 
    ggplot(ER_precise_w_seq_all_genes_tissue_specificity %>% filter(brain_derived_any == T), aes(x = n_dis_tissues)) +
    stat_ecdf(geom = "step", n = 41) +
    scale_y_continuous(name = "Cumulative proportion") +
    scale_x_continuous(name = "Number of GTEx tissues (CNS only)", breaks = seq(0, 40, 5)) + 
    theme_pubr()
  
  ER_precise_w_seq_all_genes_tissue_specificity %>% 
    filter(brain_derived_any == T) %>% 
    group_by(n_dis_tissues) %>% 
    summarise(propor_dis_tissues = n()/nrow(ER_precise_w_seq_all_genes_tissue_specificity %>% 
                                              filter(brain_derived_any == T))) %>% 
    filter(n_dis_tissues %in% 1:5) %>% 
    .[["propor_dis_tissues"]] %>% sum()
  
  ER_precise_w_seq_all_genes_tissue_specificity %>% 
    group_by(n_dis_tissues) %>% 
    summarise(propor_dis_tissues = n()/nrow(ER_precise_w_seq_all_genes_tissue_specificity)) %>% 
    filter(n_dis_tissues %in% 1:5) %>% 
    .[["propor_dis_tissues"]] %>% sum()
  
  return(ggarrange(plotlist = list(hist_tissue_specificity, cum_propor_tissue_specificity), ncol = 2, nrow = 1, align = "h"))
  
}

##### Second level #####

source("/home/dzhang/projects/constraint_conservation_wd/constraint_conservation/conservation/conservation_general/conservation_general_functions_from_bw.R")
source("/home/dzhang/projects/constraint_conservation_wd/constraint_conservation/constraint/constraint_general/constraint_general_functions.R")

# Main ------------------------------------------------------------------------------------------------

ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific <-
  ERs_w_annotation_all_tissues %>%
  filter(width > 3, !str_detect(tissue, "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast"),
         ensembl_grch38_v92_region_annot != "exon, intergenic, intron")

ER_precise_gr_w_conserv_constraint_tidy <- 
  get_conserv_constraint_prot_coding_potential_ERs(ER_precise_w_seq_all_genes, 
                                                   ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific)

# get constrained + prot potential genes 
ER_precise_w_seq_all_genes_tidy_constrained <- 
  ER_precise_gr_w_conserv_constraint_tidy %>% 
  filter(!duplicated(ER_index)) %>% 
  group_by(chr_start_end) %>% 
  mutate(ER_present_tissues = ER_tissue %>% sort() %>% unique() %>% str_c(collapse = ";")) %>% 
  ungroup() %>% 
  filter(!duplicated(chr_start_end), 
         protein_potential_any == T, 
         !is.na(mean_CDTS_percentile), 
         mean_CDTS_percentile <= 20) 

Homo_sapiens.GRCh38.92_genes_df <- 
  Homo_sapiens.GRCh38.92 %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  filter(type == "gene")

brain_specific_reannotated <- 
  Homo_sapiens.GRCh38.92_genes_df %>% 
  filter(gene_id %in% ER_precise_w_seq_all_genes_tidy_constrained$ensembl_id) %>% 
  dplyr::select(gene_chr = seqnames, gene_start = start, gene_end = end, gene_strand = strand, ensembl_gene_id = gene_id, gene_name, gene_biotype) %>% 
  left_join(gene_reannotated_properties %>% dplyr::select(ensembl_gene_id, brain_specific)) %>% 
  left_join(ER_precise_w_seq_all_genes_tidy_constrained %>% dplyr::select(protein_potential_any_frame = protein_potential_any, ER_chr = seqnames, ER_start = start, ER_end = end, ER_present_tissues,
                                                                          ensembl_gene_id = ensembl_id, 
                                                                          ER_mean_PC7 = mean_phastCons7way, ER_mean_PC20 = mean_phastCons20way, 
                                                                          ER_mean_CDTS = mean_CDTS, ER_mean_CDTS_percentile = mean_CDTS_percentile))


prot_coding_potential_kb <- plot_total_Kb_potential_prot_coding_ERs(ER_precise_w_seq_all_genes, 
                                                                    ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific)

hist_cum_propor_tissue_specificity <- plot_num_tissues_vs_ER(ER_precise_w_seq_all_genes)

# Save data -------------------------------------------------------------------------------------------

ggsave(plot = prot_coding_potential_kb, filename = "total_Kb_split_read_potential_prot_coding.png", 
       path = "OMIM_paper/figures/ERs_functional/", width = 8.27/2, height = (11.69/2), units = "in", dpi = 600)

ggsave(plot = hist_cum_propor_tissue_specificity, filename = "hist_cum_propor_tissue_specificity.png", 
       path = "OMIM_paper/supp_figures/tissue_specificity/", width = 8.27, height = (11.69/3.25), units = "in", dpi = 600)

save(ER_precise_gr_w_conserv_constraint_tidy, file = "results/check_protein_coding_potential/ER_2_split_reads_precise_w_seq_tidy.rda")

write_delim(brain_specific_reannotated, "OMIM_paper/supp_tables//reannotated_constrained_PCP_ERs.txt", delim = "\t")

load("results/check_protein_coding_potential/ER_2_split_reads_precise_w_seq_tidy.rda")
