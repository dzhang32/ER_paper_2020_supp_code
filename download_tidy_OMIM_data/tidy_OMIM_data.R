library(biomaRt)
library(tidyverse)
library(stringr)
library(forcats)
library(RColorBrewer)
library(ggpubr)
library(miscbioinfoR)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

current_date <- Sys.Date() %>% str_replace_all("-", "_")

OMIM_clinical_synopsis_prefix_num_sign_raw <- 
  read_delim(str_c("results/download_tidy_OMIM_data/OMIM_clinical_synopsis_prefix_num_sign_raw_", current_date, ".csv"), 
             delim = ",")

OMIM_phenotype_map_prefix_num_sign_raw <- 
  read_delim(str_c("results/download_tidy_OMIM_data//OMIM_phenotype_map_prefix_num_sign_raw_", current_date, ".csv"), 
             delim = ",")

pheno_abnorm_groups_df <-  
  read_delim("raw_data/pheno_abnorm_groups/pheno_abnorm_groups.csv", delim = ",")

download.file(url = "https://www.omim.org/phenotypicSeriesTitle/all?format=tsv", 
              destfile = "raw_data/OMIM_website_data/Phenotypic-Series-Titles-all.txt.tsv", mode = "wget")

pheno_series_titles <- 
  read_delim("raw_data/OMIM_website_data/Phenotypic-Series-Titles-all.txt.tsv", delim = "\t", skip = 4)

# downloaded from: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
# homo_sapiens_gene_info <- 
#     read_delim("raw_data/Homo_sapiens_031017.gene_info", delim = "\t", na = "-")

download.file(url = "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt", 
              destfile = "raw_data/hgnc_complete_gene_details/hgnc_complete_set.txt", mode = "wget")

hgnc_complete_set <- 
  read_delim("raw_data/hgnc_complete_gene_details/hgnc_complete_set.txt", delim = "\t")

download.file(url = "https://omim.org/static/omim/data/mim2gene.txt", 
              destfile = "raw_data/OMIM_website_data/mim2gene.txt", mode = "wget")

mim2gene <- 
  read_delim("raw_data/OMIM_website_data/mim2gene.txt", delim = "\t", skip = 4)

# Functions -------------------------------------------------------------------------------------------

##### First level #####

check_phenotype_map_no_data_exist <- function(OMIM_phenotype_map_prefix_num_sign_raw){
  
  OMIM_phenotype_map_prefix_num_sign_raw_no_details_exist <- 
    OMIM_phenotype_map_prefix_num_sign_raw %>% filter(details_exist == F)
  
  OMIM_clinical_synopsis_prefix_num_sign_raw_no_details_exist <- 
      OMIM_clinical_synopsis_prefix_num_sign_raw %>% 
      filter(details_exist == F)
  
  OMIM_clinical_synopsis_prefix_num_sign_raw_no_details_exist %>% 
    filter(MIM_number_pheno %in% OMIM_phenotype_map_prefix_num_sign_raw_no_details_exist$MIM_number_pheno)
  
}

analyse_pheno_abnorm_nomenclature <- function(OMIM_clinical_synopsis_prefix_num_sign_raw){
  
  clinical_synopsis_pheno_abnorm_cols <- extract_pheno_abnorm_columns(OMIM_clinical_synopsis_prefix_num_sign_raw)
  
  clinical_synopsis_pheno_abnorm_cols_no_old_format <- extract_pheno_abnorm_columns(OMIM_clinical_synopsis_prefix_num_sign_raw %>% 
                                                                                      filter(oldFormatExists == T))
  
  nomenclature_freq_plot <- 
    analyse_nomenclature(clinical_synopsis_pheno_abnorm_cols)
  
  return(nomenclature_freq_plot)
  
}

group_pheno_abnorm_columns <- function(OMIM_clinical_synopsis_prefix_num_sign_raw, pheno_abnorm_groups_df){
  
  pheno_abnorm_groups <- 
    pheno_abnorm_groups_df %>% 
    .[["pheno_abnorm_group_title"]] %>% 
    unique()
  
  OMIM_clinical_synopsis_pheno_grouped <- 
    OMIM_clinical_synopsis_prefix_num_sign_raw %>% dplyr::select(MIM_number_pheno)
  
  for(i in seq_along(pheno_abnorm_groups)){
    
    pheno_abnorm_group_to_filter <- pheno_abnorm_groups[i]
    
    pheno_abnorm_cols <- 
      pheno_abnorm_groups_df %>% 
      filter(pheno_abnorm_group_title == pheno_abnorm_group_to_filter) %>% 
      .[["OMIM_pheno_abnorm"]]
    
    pheno_abnorm_to_combine <- 
      OMIM_clinical_synopsis_prefix_num_sign_raw %>% 
      dplyr::select(one_of(pheno_abnorm_cols))
    
    if(length(pheno_abnorm_to_combine > 1)){
      
      print(str_c(i,  " - Combining ", pheno_abnorm_group_to_filter))
      
      pheno_abnorm_to_combine_tidy <- append_col_name(pheno_abnorm_to_combine)
      
      pheno_abnorm_combined <- combine_cols(pheno_abnorm_to_combine_tidy, pheno_abnorm_group_to_filter)
      
      OMIM_clinical_synopsis_pheno_grouped <- bind_cols(OMIM_clinical_synopsis_pheno_grouped, pheno_abnorm_combined)
      
    }else{
      
      colnames(pheno_abnorm_to_combine) <- pheno_abnorm_group_to_filter
      
      OMIM_clinical_synopsis_pheno_grouped <- bind_cols(OMIM_clinical_synopsis_pheno_grouped, pheno_abnorm_to_combine)
      
    }
    
  }
  
  return(OMIM_clinical_synopsis_pheno_grouped)
  
}

append_clinical_synopsis_keys <- function(OMIM_clinical_synopsis_prefix_num_sign_raw, OMIM_clinical_synopsis_pheno_grouped){
  
  clinical_synopsis_keys <- 
    OMIM_clinical_synopsis_prefix_num_sign_raw %>% 
    mutate(oldFormatExists = ifelse(is.na(oldFormatExists), FALSE, oldFormatExists)) %>% 
    dplyr::select(MIM_number_pheno, omim_pheno_preferred_title, CS_details_exist = details_exist, CS_old_format_exists = oldFormatExists, 
           CS_date_created = dateCreated, CS_date_updated = dateUpdated, PM_date_created, PM_date_updated)
  
  stopifnot(nrow(clinical_synopsis_keys) == nrow(OMIM_clinical_synopsis_pheno_grouped))
  
  OMIM_clinical_synopsis_tidy <- 
    inner_join(clinical_synopsis_keys, OMIM_clinical_synopsis_pheno_grouped)
  
  stopifnot(nrow(clinical_synopsis_keys) == nrow(OMIM_clinical_synopsis_tidy))
  
  return(OMIM_clinical_synopsis_tidy)
  
}

get_pheno_series_title <- function(OMIM_phenotype_map_prefix_num_sign_raw, pheno_series_titles){
  
  OMIM_phenotype_map_with_pheno_series_title <- 
    OMIM_phenotype_map_prefix_num_sign_raw %>% 
    mutate(pheno_series_title = NA)
  
  for(i in seq_along(OMIM_phenotype_map_prefix_num_sign_raw$phenotypicSeriesNumber)){
    
    pheno_series_num <- OMIM_phenotype_map_prefix_num_sign_raw$phenotypicSeriesNumber[i]
    
    if(is.na(pheno_series_num)){
      
      OMIM_phenotype_map_with_pheno_series_title$pheno_series_title[i] <-  NA
      
    }else{
      
      pheno_series_num_list <- str_split(pheno_series_num, ",")
      pheno_series_title_all <- character()
      
      for(j in seq_along(pheno_series_num_list[[1]])){
        
        pheno_series_num_to_search <- pheno_series_num_list[[1]][j]
        pheno_series_title <- 
          pheno_series_titles %>% 
          filter(`Phenotypic Series number` == pheno_series_num_to_search) %>% 
          .[["Phenotypic Series Title"]]
        
        
        pheno_series_title_all[j] <- pheno_series_title
        
      }
      
      OMIM_phenotype_map_with_pheno_series_title$pheno_series_title[i] <- 
        str_c(pheno_series_title_all, collapse = " & ")
      
    }
    
  }
  
  return(OMIM_phenotype_map_with_pheno_series_title)
  
}

extract_rename_phenotype_map_cols <- function(OMIM_phenotype_map_with_pheno_series_title){
  
  OMIM_phenotype_map_tidy <- 
    OMIM_phenotype_map_with_pheno_series_title %>% 
    mutate(gene_symbol = str_replace(geneSymbols, ",.*", "")) %>% 
    dplyr::select(MIM_number_pheno, omim_pheno_preferred_title, PM_details_exist = details_exist, pheno_annotated = phenotype, 
           phenotypic_series_num = phenotypicSeriesNumber, pheno_series_title, phenotype_map_key = phenotypeMappingKey, pheno_inheritence = phenotypeInheritance, 
           MIM_number_gene = mimNumber, gene_symbols = geneSymbols, gene_symbol, gene_inheritance = geneInheritance, transcript, chromosome, start = chromosomeLocationStart, 
           end = chromosomeLocationEnd, cyto_location = cytoLocation)
  
  return(OMIM_phenotype_map_tidy)
  
}

join_PM_CS <- function(OMIM_phenotype_map_tidy, OMIM_clinical_synopsis_tidy){
  
  OMIM_phenotype_map_clinical_synopsis_joined <- 
  OMIM_phenotype_map_tidy %>% 
    left_join(OMIM_clinical_synopsis_tidy) 
  
  return(OMIM_phenotype_map_clinical_synopsis_joined)
  
}

get_approved_symbols <- function(OMIM_phenotype_map_clinical_synopsis_joined, mim2gene, hgnc_complete_set){
  
  stopifnot(!any(mim2gene$`# MIM Number` %>% duplicated()))
  
  MIM_number_gene_only_PM_details_exist <- 
    OMIM_phenotype_map_clinical_synopsis_joined %>% 
    filter(PM_details_exist == T) %>% 
    .[["MIM_number_gene"]]
  
  # check if all MIM numbers downloaded are found within the mim2gene.txt database
  # seems as though mim2gene.txt is not updated as frequently as the entries themselves so have omitted this
  # stopifnot(all(MIM_number_gene_only_PM_details_exist %in% mim2gene$`# MIM Number`))
  
  # generate dataframe with column to store output gene symbol
  OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs <- 
    OMIM_phenotype_map_clinical_synopsis_joined %>% 
    left_join(mim2gene %>% 
                dplyr::select(-`MIM Entry Type (see FAQ 1.3 at https://omim.org/help/faq)`), 
              by = c("MIM_number_gene" = "# MIM Number"))
  
  indexes_to_search_for_ens_ID <- 
    OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs %>% 
    rowid_to_column() %>% 
    filter(PM_details_exist == T , is.na(`Ensembl Gene ID (Ensembl)`)) %>% 
    .[["rowid"]]
  
  for(i in indexes_to_search_for_ens_ID){
    
    gene_symbols_to_search <- 
      OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs$gene_symbols[i] %>% 
      str_split(", ") %>% 
      unlist()
    
    print(str_c(i, " - finding ENSG ID for: ", str_c(gene_symbols_to_search, collapse = ", ")))
    
    for(j in seq_along(hgnc_complete_set$symbol)){
      
      symbol_to_match <- 
        hgnc_complete_set$symbol[j] 
      
      if(any(symbol_to_match %in% gene_symbols_to_search)){
        
        OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs$`Approved Gene Symbol (HGNC)`[i] <-
          symbol_to_match
        
        ens_ID_to_store <- 
          hgnc_complete_set$ensembl_gene_id[j]
        
        if(!is.na(ens_ID_to_store)){
          
          print(str_c("matched! - ", ens_ID_to_store))
          
          OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs$`Ensembl Gene ID (Ensembl)`[i] <- 
            ens_ID_to_store
          
        }else{
          
          print(str_c("matched! - but no ens ID found"))
          
          next
          
        }
        
        break
      
      }
    }
  }
   
  return(OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs)
    
}

tidy_gene_info <- function(OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs){
  
  OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs_tidy <- 
    OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs %>% 
    dplyr::select(MIM_number_pheno, omim_pheno_preferred_title, PM_details_exist, CS_details_exist, CS_old_format_exists, CS_date_created, CS_date_updated, PM_date_created, PM_date_updated, 
           pheno_annotated, phenotypic_series_num, pheno_series_title, phenotype_map_key, pheno_inheritence, MIM_number_gene, 
           gene_symbols, gene_symbol, approved_hgnc_symbol = `Approved Gene Symbol (HGNC)`, entrez_gene_id = `Entrez Gene ID (NCBI)`, ensembl_gene_id = `Ensembl Gene ID (Ensembl)`, 
           everything())
  
  return(OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs_tidy)
  
}

filter_OMIM_pheno_gene <- function(OMIM_phenotype_map_clinical_synopsis_joined_tidy, current_date){
  
  OMIM_phenotype_map_clinical_synopsis_joined_tidy_marked <- 
    OMIM_phenotype_map_clinical_synopsis_joined_tidy %>% 
    mutate(CSNF_PM_exist = ifelse(PM_details_exist == T & CS_details_exist == T & CS_old_format_exists == F, T, F),  # CS new format genes only
           provisional = ifelse(str_detect(pheno_annotated, "\\?"), T, F), # mark provisional associations 
           nondisease = ifelse(str_detect(pheno_annotated, "\\["), T, F), # mark nondiseases
           susceptibility = ifelse(str_detect(pheno_annotated, "\\{"), T, F)) # mark susceptibility phenotypes
  
  OMIM_PM_CS_filtered <- 
    OMIM_phenotype_map_clinical_synopsis_joined_tidy_marked %>% 
    filter(CSNF_PM_exist == T, 
           provisional == F, 
           nondisease == F, 
           susceptibility == F) 
  
  print(str_c(OMIM_PM_CS_filtered$approved_hgnc_symbol %>% 
                unique() %>% 
                length(), 
              " unique genes remaining..."))
  
  OMIM_PM_CS_filtered_no_ens_id <- 
    OMIM_PM_CS_filtered %>% 
    filter(is.na(ensembl_gene_id)) %>% 
    dplyr::select(MIM_number_pheno, omim_pheno_preferred_title, MIM_number_gene, gene_symbols)
  
  write_delim(OMIM_PM_CS_filtered_no_ens_id, 
              str_c("results/download_tidy_OMIM_data/OMIM_entries_no_ens_ids_", current_date, ".csv"), delim = ",")
  
  print(str_c(length(OMIM_PM_CS_filtered_no_ens_id$MIM_number_gene), 
              " genes with no ens IDs removed from OMIM_PM_CS_filtered and saved to raw_data_tidy/CS_PM_filtered/OMIM_entries_no_ens_ids.csv"))
  
  OMIM_PM_CS_filtered_no_na_ens_ids <- 
    OMIM_PM_CS_filtered %>% 
    filter(!is.na(ensembl_gene_id))
  
  return(OMIM_PM_CS_filtered_no_na_ens_ids)
  
}

##### Second level #####

extract_pheno_abnorm_columns <- function(OMIM_clinical_synopsis_prefix_num_sign_raw){
  
  index_cols_contain_UMLS <- vector(mode = "logical", length = length(OMIM_clinical_synopsis_prefix_num_sign_raw))
  
  for(i in seq_along(OMIM_clinical_synopsis_prefix_num_sign_raw)){
    
    clinical_synopsis_col <- OMIM_clinical_synopsis_prefix_num_sign_raw[[i]]
    
    if(is.character(clinical_synopsis_col)){
      
      index_cols_contain_UMLS[i] <- clinical_synopsis_col %>% na.omit() %>% str_detect("\\{UMLS") %>% any()
      
    }else{
      
      index_cols_contain_UMLS[i] <- FALSE
      
    }
  }
   
  clinical_synopsis_pheno_abnorm_cols <- 
    OMIM_clinical_synopsis_prefix_num_sign_raw[which(index_cols_contain_UMLS)]
 
  return(clinical_synopsis_pheno_abnorm_cols)
   
}

analyse_nomenclature <- function(clinical_synopsis_pheno_abnorm_cols){
  
  pheno_abnorm_all_df <- 
  data_frame(pheno_abnorm = 
               unlist(clinical_synopsis_pheno_abnorm_cols) %>% 
               na.omit() %>% 
               str_split("\n") %>% 
               unlist()) 
  
  pheno_aborm_no_UMLS <- 
    pheno_abnorm_all_df %>% 
    filter(!str_detect(pheno_abnorm, "UMLS"))
  
  nomenclatures <- 
    c("UMLS", "SNOMEDCT", "HP:", "ICD")
  
  nomenclatures_to_search <- character()
  
  for(i in seq_along(nomenclatures)){
    
    nomenclature <- nomenclatures[i]
    
    nomenclatures_to_search <- c(nomenclatures_to_search, nomenclature)
    
    if(i < length(nomenclatures)){
      
      for(j in (i+1):length(nomenclatures)){
        
        nomenclatures_to_search <- c(nomenclatures_to_search, str_c(nomenclature,"|", nomenclatures[j]))
        
      }
    }
  }
    
  nomenclature_freq <- 
    data_frame(nomenclature = nomenclatures_to_search, 
               freq = 0)
  
  for(i in seq_along(nomenclature_freq$nomenclature)){
    
    nomenclature_to_search <- nomenclature_freq$nomenclature[i]
    
    nomenclature_freq$freq[i] <- 
      pheno_abnorm_all_df$pheno_abnorm %>% str_detect(nomenclature_to_search) %>% sum()
    
  }
  
  nomenclature_freq_proportion <- 
    nomenclature_freq %>% 
    mutate(proportion = freq/length(pheno_abnorm_all_df$pheno_abnorm))
  
  nomenclature_freq_plot <- 
    ggplot(nomenclature_freq_proportion, aes(x = fct_reorder(nomenclature, proportion, .desc = T), y = freq, label = round(proportion, digits = 2))) +
    geom_col(aes(fill = nomenclature)) +
    geom_text(nudge_y = 1500, size = 3) +
    scale_x_discrete(name = "Nomenclature") +
    scale_y_continuous(name = "Number of Phenotypic Abonormalities Labelled") +
    scale_fill_manual(values = brewer.pal(length(nomenclature_freq_proportion$nomenclature), "Spectral"), guide = FALSE) +
    ggtitle("Frequency of Different Nomenclatures within OMIM Phenotypes") +
    theme_pubr(x.text.angle = 90) + 
    theme(axis.text.x = element_text(vjust = 0.5, hjust = 1))
  
  return(nomenclature_freq_plot)
  
}

append_col_name <- function(pheno_abnorm_to_combine){
  
  pheno_abnorm_to_combine_tidy <- data_frame()
  
  for(c in seq_along(pheno_abnorm_to_combine)){
    
    pheno_abnorm_to_combine_col <- pheno_abnorm_to_combine[c]
    pheno_abnorm_to_combine_colname <- pheno_abnorm_to_combine_col %>% colnames()
   
    for(r in seq_along(pheno_abnorm_to_combine_col[[1]])){
      
      pheno_abnorm_to_combine_col_row <- pheno_abnorm_to_combine_col[[1]][r]
        
      pheno_abnorm_to_combine_col_row <- 
      str_c(pheno_abnorm_to_combine_colname, ": ", pheno_abnorm_to_combine_col_row, " @")
      
      pheno_abnorm_to_combine_col[[1]][r] <- pheno_abnorm_to_combine_col_row
      
    }
    
    if(c == 1){
      
      pheno_abnorm_to_combine_tidy <- pheno_abnorm_to_combine_col
      
    }else if (c >1){
      
      pheno_abnorm_to_combine_tidy <- bind_cols(pheno_abnorm_to_combine_tidy, pheno_abnorm_to_combine_col)
      
    }

  }

  return(pheno_abnorm_to_combine_tidy)
  
}

combine_cols <- function(pheno_abnorm_to_combine_tidy, pheno_abnorm_group){
  
  pheno_abnorm_combined_vec <- character()
  
  for(i in seq_along(pheno_abnorm_to_combine_tidy[[1]])){
    
    pheno_abnorm_to_combine_tidy_row_vec <- pheno_abnorm_to_combine_tidy[i,] %>% unlist()
    pheno_abnorm_to_combine_tidy_row_vec_no_NA <- 
      pheno_abnorm_to_combine_tidy_row_vec[!is.na(pheno_abnorm_to_combine_tidy_row_vec)]
    
    if(length(pheno_abnorm_to_combine_tidy_row_vec_no_NA) == 0){
      
      pheno_abnorm_combined_vec[i] <- NA
      
    }else if(length(pheno_abnorm_to_combine_tidy_row_vec_no_NA) == 1){
      
      pheno_abnorm_combined_vec[i] <- pheno_abnorm_to_combine_tidy_row_vec_no_NA
      
    }else if(length(pheno_abnorm_to_combine_tidy_row_vec_no_NA) > 1){
    
      pheno_abnorm_combined_vec[i] <- str_c(pheno_abnorm_to_combine_tidy_row_vec_no_NA, collapse = "\n")
      
    }
  }
  
  pheno_abnorm_combined_df <- data_frame(x = pheno_abnorm_combined_vec)
  colnames(pheno_abnorm_combined_df) <- pheno_abnorm_group
  
  return(pheno_abnorm_combined_df)
    
}

# Main ------------------------------------------------------------------------------------------------

##### Analyse nomenclature of pheno abnorm terms #####

nomenclature_freq_plot <- 
  analyse_pheno_abnorm_nomenclature(OMIM_clinical_synopsis_prefix_num_sign_raw)

##### Tidy Clinical Synopsis data #####

OMIM_clinical_synopsis_pheno_grouped <- 
  group_pheno_abnorm_columns(OMIM_clinical_synopsis_prefix_num_sign_raw, pheno_abnorm_groups_df)

OMIM_clinical_synopsis_tidy <- append_clinical_synopsis_keys(OMIM_clinical_synopsis_prefix_num_sign_raw, OMIM_clinical_synopsis_pheno_grouped)

##### Tidy Phenotype Map data #####

OMIM_phenotype_map_with_pheno_series_title <- 
  get_pheno_series_title(OMIM_phenotype_map_prefix_num_sign_raw, pheno_series_titles)

OMIM_phenotype_map_tidy <- extract_rename_phenotype_map_cols(OMIM_phenotype_map_with_pheno_series_title)

##### Join Phenotype Map and Clinical Synopsis #####

OMIM_phenotype_map_clinical_synopsis_joined <- 
  join_PM_CS(OMIM_phenotype_map_tidy, OMIM_clinical_synopsis_tidy)

OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs <- 
  get_approved_symbols(OMIM_phenotype_map_clinical_synopsis_joined, mim2gene, hgnc_complete_set)

OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs_tidy <- 
  tidy_gene_info(OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs)

##### Filter OMIM genes #####

OMIM_PM_CS_filtered <-
  filter_OMIM_pheno_gene(OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs_tidy, current_date)

##### Generate GR of OMIM genes #####

ensembl_grch38_v87_TxDb <- 
  generate_txDb_from_gtf(gtf_gff3_path = "/data/references/ensembl/gtf_gff3/v87/Homo_sapiens.GRCh38.87.gtf",
                         output_path = "/data/references/ensembl/txdb_sqlite/v87/ensembl_grch38_v87_txdb.sqlite",
                         seq_levels_to_keep = c(1:22, "X", "Y", "MT"), 
                         genome_build = "hg38")

ensembl_grch38_v87_genes_gr <- 
  genes(ensembl_grch38_v87_TxDb)

OMIM_gene_start_stop_gr <- 
  ensembl_grch38_v87_genes_gr[ensembl_grch38_v87_genes_gr$gene_id %in% unique(OMIM_PM_CS_filtered[["ensembl_gene_id"]])]

OMIM_gene_start_stop_gr_marked_overlapping_genes <- 
  mark_overlapping_genes_gr(gr_1 = OMIM_gene_start_stop_gr, gr_2 = ensembl_grch38_v87_genes_gr, identical_gr = T)
OMIM_gene_start_stop_gr_marked_overlapping_genes$overlap_any_other_gene_v92 <- 
  OMIM_gene_start_stop_gr_marked_overlapping_genes$overlap_gr2 %>% as.vector()
OMIM_gene_start_stop_gr_marked_overlapping_genes$overlap_gr2 <- NULL

# Save data -------------------------------------------------------------------------------------------

ggsave(plot = nomenclature_freq_plot, filename = "pheno_nomenclature_freq_plot.png", 
       path = "results/download_tidy_OMIM_data/", width = 8, height = 6)

write_delim(OMIM_clinical_synopsis_tidy, 
            str_c("results/download_tidy_OMIM_data/OMIM_clinical_synopsis_tidy_", current_date, ".csv"), delim = ",")

write_delim(OMIM_phenotype_map_tidy, 
            str_c("results/download_tidy_OMIM_data/OMIM_phenotype_map_tidy_", current_date, ".csv"), delim = ",")

write_delim(OMIM_phenotype_map_clinical_synopsis_joined_w_gene_IDs_tidy, 
            str_c("results/download_tidy_OMIM_data/OMIM_phenotype_map_clinical_synopsis_joined_tidy_", current_date, ".csv"), delim = ",")

write_delim(OMIM_PM_CS_filtered, 
            str_c("results/download_tidy_OMIM_data/OMIM_PM_CS_filtered_", current_date, ".csv"), delim = ",")

save(list = c("OMIM_gene_start_stop_gr_marked_overlapping_genes"), 
     file = str_c("results/download_tidy_OMIM_data/OMIM_gene_start_stop_gr_marked_overlapping_genes_ens_v92", current_date, ".Rdata"))
