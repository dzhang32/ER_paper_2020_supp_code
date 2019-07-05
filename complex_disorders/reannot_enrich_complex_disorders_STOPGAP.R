library(tidyverse)
library(stringr)
library(readxl)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

gene_reannotated_properties <- read_delim("results/analyse_ER_annotation/OMIM_gene_analysis/gene_reannotated_properties.txt", delim = "\t")

neuro_GWAS_all_STOPGAP_classified <- 
  read_xlsx("results/complex_disorders/neuro_GWAS_all_STOPGAP_MR.xlsx", sheet = 1) %>% 
  dplyr::rename(to_use = `To Use`)

load("/data/references/STOPGAP/stopgap.bestld.RData")

load("results/check_protein_coding_potential/ER_2_split_reads_precise_w_seq_tidy.rda")

# Functions -------------------------------------------------------------------------------------------

tidy_stopgap_data <- function(stopgap.bestld, ER_precise_w_seq_all_genes_tidy, gene_reannotated_properties, neuro_GWAS_all_STOPGAP_classified){
  
  stopgap.bestld_tidy <- 
  stopgap.bestld %>% 
    as_tibble() %>% 
    dplyr::select(gene.best, disease, pvalue, msh, msh.cat) %>% 
    mutate(gene.best_disease = str_c(gene.best, "_", disease)) %>% 
    filter(pvalue <= 5e-8, !duplicated(gene.best_disease)) %>%
    mutate(neuro_gene = ifelse(msh.cat == "Neurological/behavioral", T, F)) %>% # categorise by neuro 
    group_by(gene.best) %>% # use their best gene prediction
    mutate(n_dis_studies = n_distinct(disease)) %>% # add a column with the number of distinct studies 
    ungroup() %>% 
    inner_join(gene_reannotated_properties %>% filter(!duplicated(gene_name)) %>% dplyr::select(gene_name, reannotated, reannotated_2_split_read), by = c("gene.best" = "gene_name")) %>% 
    left_join(neuro_GWAS_all_STOPGAP_classified) %>% 
    mutate(to_use = ifelse(is.na(to_use), "other", to_use))
    
  return(stopgap.bestld_tidy)

}

get_parameters <- function(){
  
  to_use_cats <- c(neuro_GWAS_all_STOPGAP_classified$to_use %>% unique(), "all")
  n_dis_studies <- c("1", "any")
  reannotated <- c("any", "2_split_read")
  
  params_df <- 
    data.frame(to_use_cat = rep(to_use_cats, sum(length(n_dis_studies), length(reannotated))) %>% sort(), 
             n_dis_study = rep(n_dis_studies, length(reannotated)) %>% sort(), 
             reannotated_type = reannotated)
  
  return(params_df)
  
}

# Main ------------------------------------------------------------------------------------------------

to_use_cats <- c(neuro_GWAS_all_STOPGAP_classified$to_use %>% unique(), "all")

params_df <- 
  get_parameters() 

neuro_GWAS_all_STOPGAP_classified <- 
  neuro_GWAS_all_STOPGAP_classified %>% 
  mutate(to_use = as.character(to_use))

# tidy stopgap data
stopgap.bestld_tidy <- 
  tidy_stopgap_data(stopgap.bestld, ER_precise_w_seq_all_genes_tidy, gene_reannotated_properties, neuro_GWAS_all_STOPGAP_classified)

for(i in 1:nrow(params_df)){
  
  to_use_cat_current <- params_df$to_use_cat[i]
  n_dis_study_current <- params_df$n_dis_study[i]
  reannotated_type_current <- params_df$reannotated_type[i]
  
  print(str_c(i, " - performing fishers exact for: ", to_use_cat_current, " msh cats, ", n_dis_study_current, " num studies, ", reannotated_type_current, " reannotations"))
  
  if(to_use_cat_current != "all"){
    
    stopgap.bestld_tidy_w_to_use_cat <- 
      stopgap.bestld_tidy %>% 
      mutate(to_use_current = ifelse(to_use == to_use_cat_current, T, F)) %>% 
      group_by(gene.best) %>% 
      mutate(to_use_any = any(to_use_current)) 
    
  }else{
    
    stopgap.bestld_tidy_w_to_use_cat <- 
      stopgap.bestld_tidy %>% 
      mutate(to_use_current = ifelse(to_use %in% c("1", "2", "3"), T, F)) %>% 
      group_by(gene.best) %>% 
      mutate(to_use_any = any(to_use_current))
  
    
    # stopgap.bestld_tidy_w_to_use_cat %>% filter(to_use_current) %>% .[["to_use"]] %>% unique()
    
  }
  
  if(n_dis_study_current == "1"){
    
    stopgap.bestld_tidy_w_to_use_cat <- 
    stopgap.bestld_tidy_w_to_use_cat %>% 
      filter(n_dis_studies == 1)
    
    # stopgap.bestld_tidy_w_to_use_cat$n_dis_studies %>% unique()
    
  }
  
  if(reannotated_type_current == "any"){
    
    stopgap.bestld_tidy_w_to_use_cat_reannot_labeled <- 
      stopgap.bestld_tidy_w_to_use_cat %>% 
      dplyr::rename(reannotated_current = reannotated) 
    
  }else{
    
    stopgap.bestld_tidy_w_to_use_cat_reannot_labeled <- 
      stopgap.bestld_tidy_w_to_use_cat %>%
      dplyr::rename(reannotated_current = reannotated_2_split_read) 
  }
  
  stopgap.bestld_tidy_w_to_use_cat_n_genes <- 
    stopgap.bestld_tidy_w_to_use_cat_reannot_labeled %>% 
    group_by(reannotated_current, to_use_any) %>% 
    summarise(n_genes = n_distinct(gene.best)) %>% 
    spread(key = "to_use_any", value = "n_genes") %>% 
    ungroup() 
  
  n_neuro_vs_reannot_mat <- 
    stopgap.bestld_tidy_w_to_use_cat_n_genes %>% 
    dplyr::select(non_neuro = `FALSE`, neuro = `TRUE`) %>% 
    as.matrix() 
  
  rownames(n_neuro_vs_reannot_mat) <- c("no_reannot", "reannot")
  
  fisher_neuro_vs_reannot <- fisher.test(n_neuro_vs_reannot_mat)
  chisq_neuro_vs_reannot <- chisq.test(n_neuro_vs_reannot_mat)
  
  params_df$pvalue_fishers[i] <- fisher_neuro_vs_reannot$p.value
  params_df$pvalue_chi[i] <- chisq_neuro_vs_reannot$p.value
  params_df$reannot_neuro[i] <- stopgap.bestld_tidy_w_to_use_cat_n_genes %>% filter(reannotated_current == T) %>% .[["TRUE"]]
  params_df$reannot_non_neuro[i] <- stopgap.bestld_tidy_w_to_use_cat_n_genes %>% filter(reannotated_current == T) %>% .[["FALSE"]]
  params_df$no_reannot_neuro[i] <- stopgap.bestld_tidy_w_to_use_cat_n_genes %>% filter(reannotated_current == F) %>% .[["TRUE"]]
  params_df$no_reannot_non_neuro[i] <- stopgap.bestld_tidy_w_to_use_cat_n_genes %>% filter(reannotated_current == F) %>% .[["FALSE"]]
  
}

params_df_w_propor <- 
  params_df %>% 
  mutate(propor_reannot_neuro = reannot_neuro/(reannot_neuro + reannot_non_neuro), 
         propor_no_reannot_neuro = no_reannot_neuro/(no_reannot_neuro + no_reannot_non_neuro), 
         propor_diff = propor_reannot_neuro - propor_no_reannot_neuro)

# just to check if we split reannot into 2 and 1 split reads
# matrix_reannto <- 
#   data_frame(neuro = c(256, (715 - 256), 294), 
#              non_neuro = c(3232, (9041 - 3232), 5281)) %>% 
#   as.matrix()
# 
# rownames(matrix_reannto) <- c("reannot_2", "reannot_1", "no_reannot")
# 
# chisq.test(matrix_reannto)

reannot_complex_neuro_enrichment_table <- 
  params_df_w_propor %>%
  filter(n_dis_study == "any", reannotated_type == "any", to_use_cat != "all") %>% 
  mutate(subgroup = c("other", "neuropsychiatric", "neurodegenerative", "other neurological conditions"), 
         pvalue = p.adjust(pvalue_fishers, method = "fdr", n = length(pvalue_fishers)), 
         sign_change = ifelse(propor_diff < 0, "-", "+")) %>% 
  dplyr::select(subgroup, pvalue, sign_change)

neuro_GWAS_all_STOPGAP_classified_tidy <- 
  neuro_GWAS_all_STOPGAP_classified %>% 
  left_join(data_frame(to_use = 0:3, subgroup = c("other", "neuropsychiatric", "neurodegenerative", "other neurological conditions"))) %>% 
  filter(to_use != 0) %>% 
  dplyr::select(-to_use)


# Save data -------------------------------------------------------------------------------------------

write_delim(params_df_w_propor, "results/complex_disorders/reannot_vs_complex_neuro.csv", delim = ",")

write_delim(reannot_complex_neuro_enrichment_table, "OMIM_paper/supp_tables/reannot_complex_neuro_enrichment_table.csv", delim = ",")

write_delim(neuro_GWAS_all_STOPGAP_classified_tidy, "OMIM_paper/supp_tables/neuro_GWAS_STOPGAP_subgroup_classification.csv", delim = ",")

