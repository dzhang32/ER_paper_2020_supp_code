library(tidyverse)
library(stringr)
library(forcats)
library(ggpubr)
library(plotrix)
library(lubridate)
library(regioneR)
library(stargazer)
library(miscbioinfoR)


# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

load(file = "results/annotate_ERs/ERs_optimised_cut_off_max_gap_all_tissues_w_annot_df.rda")

OMIM_PM_CS_filtered <- read_delim("results/download_tidy_OMIM_data/OMIM_PM_CS_filtered_2018_05_29.csv", delim = ",")

pheno_abnorm_groups_gtex_tissue_matched_TPC <- read_delim("raw_data/pheno_abnorm_groups/pheno_abnorm_groups_gtex_tissue_matched_TPC.txt", delim = "\t")

OMIM_gtex_tissue_name_formatting <- read_delim("raw_data/gtex_tissue_name_formatting/OMIM_gtex_tissue_name_formatting.csv", delim = ",")

gtex_gene_tpm_by_tissue <- read_delim("/home/dzhang/data/gtex/gene_tpm/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct", delim = "\t", skip = 2)

funicane_tstats <- read_delim("/home/dzhang/data/finucane_2018/tstats/GTEx.tstat.tsv", delim = "\t")

load("/data/references/STOPGAP/stopgap.bestld.RData")

load("results/check_protein_coding_potential/ER_2_split_reads_precise_w_seq_tidy.rda")

# Functions -------------------------------------------------------------------------------------------

##### First level #####

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_txDb_from_gtf.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_genomic_state.R")

get_list_brain_specific_genes <- function(funicane_tstats, top_propor = 0.1){
  
  funicane_tstats_brain <- 
  funicane_tstats %>% 
    gather(key = "tissue", value = "t-stat", -ENSGID) %>% 
    filter(str_detect(tissue, "Brain"))
  
  funicane_tstats_brain_percent_filtered <- 
    funicane_tstats_brain %>% 
    group_by(tissue) %>% 
    mutate(percentile_tstat = percent_rank(`t-stat`)) %>% 
    ungroup() %>% 
    filter(percentile_tstat >= (1 - top_propor))
  
  finucane_brain_specific_genes <- 
    funicane_tstats_brain_percent_filtered$ENSGID %>% unique()
  
  return(finucane_brain_specific_genes)
  
}

get_gene_properties <- function(ensembl_grch38_v92_genes_txdb, ensembl_grch38_v92_gtf_gr, ensembl_grch38_v92_genes_genomic_state, OMIM_PM_CS_filtered, ER_precise_gr_w_conserv_constraint_tidy,
                                ERs_w_annotation_all_tissues_intron_inter, gtex_gene_tpm_by_tissue, OMIM_gtex_tissue_name_formatting, finucane_brain_specific_genes){
  
  # get all genes 
  ensembl_grch38_v92_genes_txdb_genes <- genes(ensembl_grch38_v92_genes_txdb)
  
  # get all overlapping genes
  ensembl_grch38_v92_genes_txdb_genes_overlapping_genes <- 
    mark_overlapping_genes_gr(gr_1 = ensembl_grch38_v92_genes_txdb_genes, gr_2 = ensembl_grch38_v92_genes_txdb_genes, identical_gr = T, maxgap = -1, minoverlap = 1) 
  
  ensembl_grch38_v92_genes_txdb_genes_overlapping_genes <- 
    ensembl_grch38_v92_genes_txdb_genes_overlapping_genes[ensembl_grch38_v92_genes_txdb_genes_overlapping_genes$overlap_gr2 == T]
  
  # get all reannotated genes
  uniq_genes_split_read_annot_all <- 
    ERs_w_annotation_all_tissues_intron_inter %>% 
    group_by(uniq_genes_split_read_annot) %>% 
    summarise(ensembl_grch38_v92_region_annot = str_c(unique(ensembl_grch38_v92_region_annot), collapse = ", "))
  
  # get biotype and transcript count
  ensembl_grch38_v92_gene_biotype_transcript_count <- 
    ensembl_grch38_v92_gtf_gr %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    group_by(gene_id, gene_biotype) %>% 
    summarise(transcript_count = sum(type == "transcript")) 
  
  ensembl_grch38_v92_gene_biotype_desc_n_genes <- 
    ensembl_grch38_v92_gene_biotype_transcript_count %>% 
    group_by(gene_biotype) %>% 
    summarise(n_genes = n()) %>% 
    arrange(desc(n_genes)) %>% 
    .[1:5,]
  
  # get gene names for all ensembl ids
  ensembl_grch38_v92_gene_name_to_ens_id <- 
    ensembl_grch38_v92_gtf_gr[ensembl_grch38_v92_gtf_gr$type == "gene"] %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    dplyr::select(gene_name, gene_id)
  
  ensembl_grch38_v92_gene_biotype_transcript_count <- 
    ensembl_grch38_v92_gene_biotype_transcript_count %>% 
    mutate(gene_biotype = ifelse(gene_biotype %in% ensembl_grch38_v92_gene_biotype_desc_n_genes$gene_biotype, gene_biotype, "other"))
  
  # get longest intron
  # ensembl_grch38_v92_genes_longest_intron <-
  #   get_longest_intron_length(ensembl_grch38_v92_genes_genomic_state, ensembl_grch38_v92_genes_txdb)
  
  # get average expression across all tissues
  gtex_gene_tpm_averaged <- 
    gtex_gene_tpm_by_tissue %>% 
    gather(key = "tissue", value = "tpm", -gene_id, -Description) 
  
  stopifnot(all(unique(gtex_gene_tpm_averaged$tissue) %in% OMIM_gtex_tissue_name_formatting$gtex_tissues_name))
  
  gtex_gene_tpm_averaged <- 
    gtex_gene_tpm_averaged %>% 
    left_join(OMIM_gtex_tissue_name_formatting %>% dplyr::select(gtex_tissues_name, OMIM_gtex_name), by = c("tissue" = "gtex_tissues_name")) %>% 
    filter(!str_detect(OMIM_gtex_name, "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast"), tpm >= 0.1) %>% 
    mutate(gene_id = str_replace(gene_id, "\\..*", "")) %>% 
    group_by(gene_id) %>% 
    summarise(mean_tpm = mean(tpm))
  
  # get neuro pheno (OMIM)
  OMIM_gene_reannot_neuro_pheno <- 
    OMIM_PM_CS_filtered %>% 
    group_by(ensembl_gene_id) %>% 
    mutate(neuro_pheno = ifelse(!is.na(neurologic), T, F)) %>% 
    summarise(neuro_pheno_by_gene = any(neuro_pheno))
  
  # get date created (OMIM)
  OMIM_gene_PM_data_created <- 
    OMIM_PM_CS_filtered %>% 
    mutate(PM_date_created_formatted = 
             PM_date_created %>% 
             str_sub(13, 17) %>%
             as.integer()) %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(PM_date_created_formatted = min(PM_date_created_formatted))
  
  gene_reannotated_properties <- 
    data_frame(ensembl_gene_id = ensembl_grch38_v92_genes_txdb_genes$gene_id, 
               gene_length = ensembl_grch38_v92_genes_txdb_genes %>% width(),
               reannotated = ensembl_gene_id %in% uniq_genes_split_read_annot_all$uniq_genes_split_read_annot, 
               reannotated_2_split_read = ensembl_gene_id %in% (ER_precise_gr_w_conserv_constraint_tidy$ensembl_id %>% unique()),
               OMIM_gene = ensembl_gene_id %in% OMIM_PM_CS_filtered$ensembl_gene_id, 
               brain_specific = ensembl_gene_id %in% finucane_brain_specific_genes, 
               overlapping_gene = ensembl_gene_id %in% ensembl_grch38_v92_genes_txdb_genes_overlapping_genes$gene_id) %>% 
    left_join(ensembl_grch38_v92_gene_name_to_ens_id, by = c("ensembl_gene_id" = "gene_id")) %>% 
    left_join(uniq_genes_split_read_annot_all, by = c("ensembl_gene_id" = "uniq_genes_split_read_annot")) %>% 
    left_join(ensembl_grch38_v92_gene_biotype_transcript_count, by = c("ensembl_gene_id" = "gene_id")) %>% 
    # left_join(ensembl_grch38_v92_genes_longest_intron, by = c("ensembl_gene_id" = "gene")) %>% 
    left_join(gtex_gene_tpm_averaged, by = c("ensembl_gene_id" = "gene_id")) %>% 
    left_join(OMIM_gene_reannot_neuro_pheno) %>% 
    left_join(OMIM_gene_PM_data_created) %>% 
     mutate(mean_tpm = ifelse(is.na(mean_tpm), 0, mean_tpm), 
            gene_prot_coding = ifelse(gene_biotype == "protein_coding", T, F))
  
  all_genes_reannot_model <- 
    glm(formula = reannotated ~ OMIM_gene + brain_specific + mean_tpm + overlapping_gene + transcript_count + gene_biotype + gene_length, 
        data = gene_reannotated_properties)
  
  return(gene_reannotated_properties)
  
}

analyse_all_genes_reannot <- function(gene_reannotated_properties){
  
  OMIM_vs_non_OMIM_gene_percent_reannot_plot <- 
    ggplot(gene_reannotations %>% group_by(reannotated) %>% summarise(percent_OMIM_gene = mean(OMIM_gene) * 100), 
           aes(x = reannotated, y = percent_OMIM_gene, fill = reannotated)) +
    geom_col(color="black", width = 0.75) +
    scale_x_discrete(name = "Gene reannotated") +
    scale_y_continuous(name = "OMIM gene (%)") +
    scale_fill_manual(values = get_palette("jco", 3)[c(3,1)], guide = F) +
    theme_bw()
  
  all_genes_reannot_percent_protein_coding_plot <- 
    ggplot(gene_reannotations %>% group_by(reannotated) %>% summarise(percent_protein_coding = mean(gene_biotype == "protein_coding") * 100), 
           aes(x = reannotated, y = percent_protein_coding, fill = reannotated)) +
    geom_col(color="black", width = 0.75) +
    scale_x_discrete(name = "Gene reannotated") +
    scale_y_continuous(name = "Protein coding (%)") +
    scale_fill_manual(values = get_palette("jco", 3)[c(3,1)], guide = F) +
    theme_bw()
  
  all_genes_transcript_count_plot <- 
    ggplot(gene_reannotations, 
           aes(x = reannotated, y = transcript_count, fill = reannotated)) +
    geom_bar(stat = "summary", fun.y = mean, color="black", width = 0.75) +
    geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2) +
    scale_x_discrete(name = "Gene reannotated") +
    scale_y_continuous(name = "Transcript count") +
    scale_fill_manual(values = get_palette("jco", 3)[c(3,1)], guide = F) +
    stat_compare_means(comparisons = list(c("TRUE", "FALSE")), label = "p.signif", label.y = 11, tip.length = 0.0075) +
    theme_bw() 
  
  list_all_genes_reannot_plots <- 
    ggarrange(plotlist = list(OMIM_vs_non_OMIM_gene_percent_reannot_plot, gene_reannot_percent_protein_coding_plot, all_genes_transcript_count_plot), 
              ncol = 3, nrow = 1)
  
}

OMIM_genes_reannot_propor_neuro_pheno <- function(OMIM_PM_CS_filtered, ERs_w_annotation_all_tissues_intron_inter){
  
  OMIM_genes_reannot <-
    ERs_w_annotation_all_tissues_intron_inter %>% 
    filter(!is.na(connect_to_OMIM_gene_split_read_annot)) %>% 
    .[["connect_to_OMIM_gene_split_read_annot"]] %>% unique()
  
  OMIM_gene_reannot_neuro_pheno <- 
  OMIM_PM_CS_filtered %>% 
    mutate(reannotated = ifelse(ensembl_gene_id %in% OMIM_genes_reannot, T, F), 
           neuro_pheno = ifelse(!is.na(neurologic), T, F)) %>% 
    group_by(ensembl_gene_id, reannotated) %>% 
    summarise(neuro_pheno_by_gene = any(neuro_pheno)) 
  
  OMIM_gene_reannot_neuro_pheno %>% 
    filter(reannotated == T) %>% 
    .[["neuro_pheno_by_gene"]] %>% 
    mean()
  
  return(list(OMIM_gene_percent_neuro_reannot_plot, OMIM_gene_neuro_reannot_quant_plot))
  
}

analyse_OMIM_gene_reannot <- function(OMIM_PM_CS_filtered, ERs_w_annotation_all_tissues_intron_inter_OMIM_genes, 
                                      ERs_w_annotation_all_tissues_intron_inter){
  
  OMIM_gene_reannot_neuro_pheno <- 
    OMIM_PM_CS_filtered %>% 
    mutate(reannotated = ifelse(ensembl_gene_id %in% ERs_w_annotation_all_tissues_intron_inter_OMIM_genes, T, F)) %>% 
    group_by(ensembl_gene_id, reannotated) %>% 
    mutate(neuro_pheno = ifelse(!is.na(neurologic), T, F)) %>% 
    summarise(neuro_pheno_by_gene = any(neuro_pheno))
  
  OMIM_gene_PM_data_created <- 
    OMIM_PM_CS_filtered %>% 
    mutate(PM_date_created_formatted = 
             PM_date_created %>% 
             str_sub(13, 17) %>%
             as.integer()) %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(PM_date_created_formatted = min(PM_date_created_formatted))
  
  OMIM_gene_biotype_count_length <- 
    query_biomart(mart = 38, attributes = c("ensembl_gene_id", "transcript_count", "gene_biotype", "percentage_gene_gc_content"), 
                  filters = "ensembl_gene_id", values = OMIM_gene_reannot_neuro_pheno$ensembl_gene_id) %>% 
    as_tibble() 
  
  OMIM_gene_annotations <- 
    OMIM_gene_reannot_neuro_pheno %>% 
    left_join(OMIM_gene_PM_data_created) %>% 
    left_join(OMIM_gene_biotype_count_length)
  
  OMIM_gene_reannot_model <- 
    glm(formula = reannotated ~ transcript_count + neuro_pheno_by_gene + PM_date_created_formatted + gene_biotype + percentage_gene_gc_content, 
        data = OMIM_gene_annotations)
    
  print(summary(OMIM_gene_reannot_model))
  
  OMIM_gene_transcript_count_plot <- 
    ggplot(OMIM_gene_annotations, 
           aes(x = reannotated, y = transcript_count, fill = reannotated)) +
    geom_bar(stat = "summary", fun.y = mean, color="black", width = 0.75) +
    geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2) +
    scale_x_discrete(name = "Gene reannotated") +
    scale_y_continuous(name = "Transcript count") +
    scale_fill_manual(values = get_palette("jco", 3)[c(3,1)], guide = F) +
    stat_compare_means(comparisons = list(c("TRUE", "FALSE")), label = "p.signif", label.y = 13, tip.length = 0.0075) +
    theme_bw() 
  
  OMIM_gene_percent_GC_plot <- 
    ggplot(OMIM_gene_annotations, 
           aes(x = reannotated, y = percentage_gene_gc_content, fill = reannotated)) +
    geom_bar(stat = "summary", fun.y = mean, color="black", width = 0.75) +
    geom_errorbar(stat = "summary", fun.data = "mean_se", width = 0.2) +
    scale_x_discrete(name = "Gene reannotated") +
    scale_y_continuous(name = "Percentage GC content") +
    scale_fill_manual(values = get_palette("jco", 3)[c(3,1)], guide = F) +
    stat_compare_means(comparisons = list(c("TRUE", "FALSE")), label = "p.signif", label.y = 50, tip.length = 0.0075) +
    theme_bw() 
  
  list_OMIM_reannot_plots <- 
    ggarrange(plotlist = list(OMIM_gene_transcript_count_plot, OMIM_gene_percent_GC_plot), ncol = 2, nrow = 1)
  
  return(list_OMIM_reannot_plots)
  
}

check_OMIM_gene_reannot_match_gtex_tissue <- 
  function(ERs_w_annotation_all_tissues_intron_inter, OMIM_PM_CS_filtered, pheno_abnorm_groups_gtex_tissue_matched_TPC, OMIM_gtex_tissue_name_formatting){
  
    OMIM_genes_pheno_abnorm_w_gtex_tissue <- 
    match_OMIM_gene_pheno_abnorm_to_GTEx_tissue(OMIM_PM_CS_filtered, pheno_abnorm_groups_gtex_tissue_matched_TPC, OMIM_gtex_tissue_name_formatting)
    
    # check all pheno matched tissues are named the same as in our ER annotation
    pheno_matched_tissues_uniq <- 
      OMIM_genes_pheno_abnorm_w_gtex_tissue %>% 
      filter(!is.na(gtex_tissue_pheno_matched)) %>% 
      .[["gtex_tissue_pheno_matched"]] %>% 
      str_split(", ") %>% 
      unlist() %>% 
      unique()
    
    stopifnot(pheno_matched_tissues_uniq %in% OMIM_gtex_tissue_name_formatting$OMIM_gtex_name)
    
    OMIM_genes_ER_expressed_tissue_pheno_match <- 
      ERs_w_annotation_all_tissues_intron_inter %>% 
      filter(!is.na(connect_to_OMIM_gene_split_read_annot)) %>% 
      group_by(connect_to_OMIM_gene_split_read_annot) %>% 
      summarise(ER_expressed_in = str_c(unique(tissue), collapse = ", ")) %>% 
      left_join(OMIM_genes_pheno_abnorm_w_gtex_tissue, by = c("connect_to_OMIM_gene_split_read_annot" = "ensembl_gene_id")) 
   
    for(i in 1:nrow(OMIM_genes_ER_expressed_tissue_pheno_match)){
      
      ER_expressed_in_tissues_uniq <- 
      OMIM_genes_ER_expressed_tissue_pheno_match$ER_expressed_in[i] %>% 
        str_split(", ") %>% unlist() %>% unique()
      
      gtex_tissue_pheno_matched_uniq <- 
        OMIM_genes_ER_expressed_tissue_pheno_match$gtex_tissue_pheno_matched[i] %>% 
        str_split(", ") %>% unlist() %>% unique()
      
      OMIM_genes_ER_expressed_tissue_pheno_match$percent_ER_tissue_in_matched[[i]] <- mean(ER_expressed_in_tissues_uniq %in% gtex_tissue_pheno_matched_uniq)
      OMIM_genes_ER_expressed_tissue_pheno_match$percent_matched_tissue_in_ER[[i]] <- mean(gtex_tissue_pheno_matched_uniq %in% ER_expressed_in_tissues_uniq)
      OMIM_genes_ER_expressed_tissue_pheno_match$any_matched_tissue_in_ER[[i]] <- any(gtex_tissue_pheno_matched_uniq %in% ER_expressed_in_tissues_uniq)
      
    }
    
    OMIM_genes_ER_expressed_tissue_pheno_match %>% filter(!is.na(gtex_tissue_pheno_matched)) %>% .[["percent_ER_tissue_in_matched"]] %>% mean()
    OMIM_genes_ER_expressed_tissue_pheno_match %>% filter(!is.na(gtex_tissue_pheno_matched)) %>% .[["percent_matched_tissue_in_ER"]] %>% mean()
    OMIM_genes_ER_expressed_tissue_pheno_match %>% filter(!is.na(gtex_tissue_pheno_matched)) %>% .[["any_matched_tissue_in_ER"]] %>% mean()
    
    return(OMIM_genes_ER_expressed_tissue_pheno_match)
    
  }

get_OMIM_gene_logis_regress_table <- function(all_genes_reannot_model){
  
  all_genes_reannot_model_table <- 
  coef(summary(all_genes_reannot_model)) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "gene_property") %>% 
    as_tibble() %>% 
    mutate(`Odds ratio` = exp(Estimate) %>% round(digits = 3), 
           `p-value` = ifelse(`Pr(>|t|)` < 1e-5, "***", round(`Pr(>|t|)`, digits = 3)), 
           `Gene property` = c("Constant", "Brain-specific", "Gene TPM", "Overlapping gene", "Transcript count", 
                                 "Gene biotype - lincRNA", "Gene biotype - other", "Gene biotype - processed pseudogene", "Gene biotype - protein coding", "Gene biotype - unprocessed pseudogene", "Gene length") %>% 
             factor() %>% fct_relevel(c("Brain-specific", "Transcript count", "Gene length", 
                                        "Gene biotype - protein coding", "Gene biotype - lincRNA", "Gene biotype - processed pseudogene", "Gene biotype - unprocessed pseudogene", "Gene biotype - other",
                                        "Gene TPM", "Overlapping gene"))) %>% 
    arrange(`Gene property`) %>% 
    dplyr::select(`Gene property`, `Odds ratio`, `p-value`, everything())
  
  return(all_genes_reannot_model_table)
  
  
}

plot_propor_OMIM_genes_reannotated <- function(gene_reannotated_properties, OMIM_genes_ER_expressed_tissue_pheno_match){
  
  gene_reannotated_properties_w_propor <- 
  gene_reannotated_properties %>% 
    filter(OMIM_gene == T) %>% 
    group_by(reannotated) %>% 
    summarise(n = n()) %>% 
    mutate(propor = n/2898, 
           n_pheno = c(72), 
           reannotated = reannotated %>% factor() %>% fct_relevel(c("TRUE", "FALSE")))
  
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )
  
  gene_tissue_matched_propor <- 
    OMIM_genes_ER_expressed_tissue_pheno_match %>% 
    group_by(any_matched_tissue_in_ER) %>% 
    summarise(n = n(), 
              propor = n/nrow(OMIM_genes_ER_expressed_tissue_pheno_match)) %>% 
    mutate(any_matched_tissue_in_ER = any_matched_tissue_in_ER %>% factor() %>% fct_relevel(c("TRUE", "FALSE")))
  
  OMIM_gene_tissue_matched_propor_plot <- 
    ggplot(gene_tissue_matched_propor, 
         aes(x = "", y = n, fill = any_matched_tissue_in_ER)) +
    geom_col(width = 1, colour = "black") + 
    coord_polar("y", start=0) + 
    scale_fill_manual(name = "Reannotated OMIM genes", values = c("black", "white")) + 
    blank_theme + 
    theme(axis.text.x=element_blank())
  
  OMIM_gene_reannotated_propor_plot <- 
    ggplot(gene_reannotated_properties_w_propor, 
         aes(x = "", y = n, fill = reannotated)) +
    geom_col(width = 1, colour = "black") + 
    coord_polar("y", start=0) + 
    scale_fill_manual(name = "OMIM genes", values = get_palette("jco", 3)[c(1,3)]) + 
    blank_theme + 
    theme(axis.text.x=element_blank())
  
  return(list(OMIM_gene_tissue_matched_propor_plot, OMIM_gene_reannotated_propor_plot))
  
  
}


##### Second level #####

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/query_biomart.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/mark_overlapping_genes_gr.R")

match_OMIM_gene_pheno_abnorm_to_GTEx_tissue <- function(OMIM_PM_CS_filtered, pheno_abnorm_groups_gtex_tissue_matched_TPC, OMIM_gtex_tissue_name_formatting){
  
  gtex_tissues_formatted_uniq <- 
    pheno_abnorm_groups_gtex_tissue_matched_TPC$gtex_tissue_formatted %>% 
    na.omit() %>%
    str_split( ", ") %>% 
    unlist() %>% 
    unique()
  
  stopifnot(all(gtex_tissues_formatted_uniq %in% OMIM_gtex_tissue_name_formatting$OMIM_gtex_name))
  
  OMIM_genes_pheno_abnorm <- data_frame()
  
  print(str_c(Sys.time(), " - getting pheno abnorms that exist for each gene"))
  
  for(i in seq_along(unique(OMIM_PM_CS_filtered$ensembl_gene_id))){
    
    ensembl_gene_id_to_check <- unique(OMIM_PM_CS_filtered$ensembl_gene_id)[i]
    
    print(str_c(Sys.time(), " - ", i, " - ", ensembl_gene_id_to_check))
    
    OMIM_PM_CS_filtered_gene_filtered <- 
      OMIM_PM_CS_filtered %>% 
      filter(ensembl_gene_id == ensembl_gene_id_to_check)
    
    pheno_aborm_cols_unique <- 
      OMIM_PM_CS_filtered_gene_filtered[unique(pheno_abnorm_groups_gtex_tissue_matched_TPC$pheno_abnorm_group_title)] %>% 
      unlist() %>% 
      str_split(" @") %>% 
      unlist() %>% 
      str_replace_all("\n", "") %>% 
      str_replace(":.*", "") %>% 
      unique()
    
    pheno_aborm_cols_unique_no_NA <- pheno_aborm_cols_unique[!is.na(pheno_aborm_cols_unique)]
    pheno_aborm_cols_unique_no_NA_no_blank <- pheno_aborm_cols_unique_no_NA[pheno_aborm_cols_unique_no_NA != ""]
    
    stopifnot(all(pheno_aborm_cols_unique_no_NA_no_blank %in% pheno_abnorm_groups_gtex_tissue_matched_TPC$OMIM_pheno_abnorm))
    
    ensembl_gene_id_pheno_aborm <- 
      data_frame(ensembl_gene_id = ensembl_gene_id_to_check, 
               OMIM_pheno_abnorm = pheno_aborm_cols_unique_no_NA_no_blank)
    
    OMIM_genes_pheno_abnorm <- 
      bind_rows(OMIM_genes_pheno_abnorm, ensembl_gene_id_pheno_aborm)
    
  }
  
  
  OMIM_genes_pheno_abnorm_w_gtex_tissue <- 
    OMIM_genes_pheno_abnorm %>% 
    left_join(pheno_abnorm_groups_gtex_tissue_matched_TPC %>% dplyr::select(OMIM_pheno_abnorm, gtex_tissue_formatted)) %>% 
    group_by(ensembl_gene_id) %>% 
    summarise(gtex_tissue_pheno_matched = remove_NA_uniq(gtex_tissue_formatted)) 
  
  return(OMIM_genes_pheno_abnorm_w_gtex_tissue)
  
}

get_longest_intron_length <- function(ensembl_grch38_v92_genes_genomic_state_full_genome, ensembl_grch38_v92_genes_txdb){
  
  ensembl_grch38_v92_genes_txdb_genes <- genes(ensembl_grch38_v92_genes_txdb)
  
  ensembl_grch38_v92_genes_genomic_state_full_genome <- ensembl_grch38_v92_genes_genomic_state$fullGenome
  
  stopifnot(max(ensembl_grch38_v92_genes_genomic_state_full_genome$gene %>% unlist()) == length(genes(ensembl_grch38_v92_genes_txdb)))
  
  ensembl_grch38_v92_genes_genomic_state_full_genome_intron <- 
    ensembl_grch38_v92_genes_genomic_state_full_genome[ensembl_grch38_v92_genes_genomic_state_full_genome$theRegion == "intron"]
  
  ensembl_grch38_v92_genes_genomic_state_full_genome[ensembl_grch38_v92_genes_genomic_state_full_genome$theRegion == "intron" & 
                                                       ensembl_grch38_v92_genes_genomic_state_full_genome$gene %in% 29684]
  
  all_introns_unlisted_index <- rep(seq_along(ensembl_grch38_v92_genes_genomic_state_full_genome_intron$gene), lapply(ensembl_grch38_v92_genes_genomic_state_full_genome_intron$gene, length) %>% unlist())
  
  ensembl_grch38_v92_genes_genomic_state_full_genome_intron_unlist <- ensembl_grch38_v92_genes_genomic_state_full_genome_intron[all_introns_unlisted_index] 
  
  ensembl_grch38_v92_genes_genomic_state_full_genome_intron_unlist$gene_unlist <- ensembl_grch38_v92_genes_genomic_state_full_genome_intron$gene %>% unlist()
  
  ensembl_grch38_v92_genes_genomic_state_full_genome_intron_unlist$gene_id_unlist <- genes(ensembl_grch38_v92_genes_txdb)$gene_id[ensembl_grch38_v92_genes_genomic_state_full_genome_intron_unlist$gene_unlist]
  
  names(ensembl_grch38_v92_genes_genomic_state_full_genome_intron_unlist) <- NULL
  
  ensembl_grch38_v92_genes_longest_intron <- 
    ensembl_grch38_v92_genes_genomic_state_full_genome_intron_unlist %>% 
    as.data.frame() %>% 
    as_tibble() %>%
    group_by(gene_id_unlist) %>% 
    summarise(longest_intron_length = max(width)) %>% 
    dplyr::rename(gene = gene_id_unlist) 

  return(ensembl_grch38_v92_genes_longest_intron)
  
}

##### Third level ######

remove_NA_uniq <- function(x){
  
  y <- x[!is.na(x)]
  
  if(length(y) == 0){
    
    return(as.character(NA))
    
  }else{
    
    z <- 
      y %>% 
      str_split(", ") %>% 
      unlist() %>% 
      unique() %>% 
      str_c(collapse = ", ")
    
    return(z)
    
  }
  
}

# Main ------------------------------------------------------------------------------------------------

ensembl_grch38_v92_genes_txdb <- 
  generate_txDb_from_gtf(gtf_gff3_path = "/data/references/ensembl/gtf_gff3/v92/Homo_sapiens.GRCh38.92.gtf", 
                         output_path = "/data/references/ensembl/txdb_sqlite/v92/ensembl_grch38_v92_txdb.sqlite",
                         seq_levels_to_keep = c(1:22, "X", "Y", "MT"), genome_build = "hg38")

ensembl_grch38_v92_genes_genomic_state <- 
  generate_genomic_state(ensembl_grch38_v92_genes_txdb, 
                         output_path = "/data/references/ensembl/genomic_state/v92/ensembl_grch38_v92_genomic_state.rda")

ensembl_grch38_v92_gtf_gr <- 
  import("/data/references/ensembl/gtf_gff3/v92/Homo_sapiens.GRCh38.92.gtf")

ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific <-
  ERs_w_annotation_all_tissues %>%
  filter(width > 3, !str_detect(tissue, "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast"),
         ensembl_grch38_v92_region_annot != "exon, intergenic, intron")

ERs_w_annotation_all_tissues_intron_inter <- 
  ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific %>% 
  filter(ensembl_grch38_v92_region_annot %in% c("intron", "intergenic") ,  
         !is.na(uniq_genes_split_read_annot), !str_detect(uniq_genes_split_read_annot, ","))

finucane_brain_specific_genes <- get_list_brain_specific_genes(funicane_tstats, top_propor = 0.10)

gene_reannotated_properties <- 
  get_gene_properties(ensembl_grch38_v92_genes_txdb, ensembl_grch38_v92_gtf_gr, ensembl_grch38_v92_genes_genomic_state, OMIM_PM_CS_filtered, ER_precise_gr_w_conserv_constraint_tidy, 
                      ERs_w_annotation_all_tissues_intron_inter, gtex_gene_tpm_by_tissue, OMIM_gtex_tissue_name_formatting, finucane_brain_specific_genes)

# analyse and print logis regress on all genes 
all_genes_reannot_model <- 
  glm(formula = reannotated ~ brain_specific + mean_tpm + overlapping_gene + transcript_count + gene_biotype + gene_length, 
      data = gene_reannotated_properties)

summary(all_genes_reannot_model)

all_genes_reannot_model_table <- get_OMIM_gene_logis_regress_table(all_genes_reannot_model)

write_delim(all_genes_reannot_model_table, "OMIM_paper/tables/OMIM_gene_logis_regress/logis_regress_all_gene_properties.csv", delim = ",")

OMIM_genes_ER_expressed_tissue_pheno_match <- 
  check_OMIM_gene_reannot_match_gtex_tissue(ERs_w_annotation_all_tissues_intron_inter, OMIM_PM_CS_filtered, pheno_abnorm_groups_gtex_tissue_matched_TPC, OMIM_gtex_tissue_name_formatting)

list_OMIM_reannot_tissue_match_plots <- plot_propor_OMIM_genes_reannotated(gene_reannotated_properties, OMIM_genes_ER_expressed_tissue_pheno_match)

reannotated_genes_w_ppc_ER <- 
  ER_precise_gr_w_conserv_constraint_tidy %>% 
  filter(protein_potential_any == T) %>% 
  .[["ensembl_id"]] %>% 
  unique()

mean(unique(OMIM_PM_CS_filtered$ensembl_gene_id) %in% reannotated_genes_w_ppc_ER)

# gene_reannotated_properties %>% 
#   mutate(ensembl_grch38_v92_region_annot = ifelse(ensembl_grch38_v92_region_annot == "intergenic, intron", 
#                                                   "intron, intergenic", 
#                                                   ensembl_grch38_v92_region_annot)) %>% 
#   write_delim(path = "results_tmp/gene_list_ens_v92.csv", delim = ",")
# filter(reannotated == T, brain_specific == T, gene_prot_coding == T) %>% 
#   .[["ensembl_grch38_v92_region_annot"]] %>% 
#   table()

# Save data -------------------------------------------------------------------------------------------

write_delim(OMIM_genes_ER_expressed_tissue_pheno_match, 
            "results/analyse_ER_annotation/OMIM_gene_analysis/OMIM_genes_ER_expressed_tissue_pheno_match.csv", 
            delim = ",")

write_delim(gene_reannotated_properties, "results/analyse_ER_annotation/OMIM_gene_analysis/gene_reannotated_properties.txt", delim = "\t")

ggsave(plot = list_OMIM_reannot_tissue_match_plots[[1]], 
       path = "OMIM_paper/figures/OMIM_example/", 
       filename = "OMIM_tissue_match.png", width = 8, height = 6, dpi = 600)

ggsave(plot = list_OMIM_reannot_tissue_match_plots[[2]], 
       path = "OMIM_paper/figures/OMIM_example/", 
       filename = "OMIM_reannot_propor.png", width = 8, height = 6, dpi = 600)
