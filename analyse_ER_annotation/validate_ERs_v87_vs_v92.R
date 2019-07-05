library(tidyverse)
library(stringr)
library(forcats)
library(ggpubr)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

load(file = "results/annotate_ERs/ERs_optimised_cut_off_max_gap_all_tissues_w_annot_df.rda")

gtex_tissue_name_formatting <- read_delim("raw_data/gtex_tissue_name_formatting/OMIM_gtex_tissue_name_formatting.csv", delim = ",")

# Functions -------------------------------------------------------------------------------------------

##### First level ##### 

get_ensembl_v87_to_v92_tidy <- function(ERs_w_annotation_all_tissues, feature_type_order_colors, gtex_tissue_name_formatting){
  
  gtex_tissue_name_formatting_w_face <- 
    gtex_tissue_name_formatting %>% 
    mutate(face = ifelse(gtex_tissue_group == "brain", "bold", "plain"))
  
  ERs_v87_vs_v92 <- 
    ERs_w_annotation_all_tissues %>% 
    mutate(ensembl_v87_v92_annot = 
             ifelse(ensembl_grch38_v87_region_annot == "intron" & ensembl_grch38_v92_region_annot %in% c("exon"), 
                    yes = "intron in v87", 
                    no = ifelse(ensembl_grch38_v87_region_annot == "intergenic" & ensembl_grch38_v92_region_annot %in% c("exon"), 
                                yes = "intergenic in v87", 
                                no = ifelse(ensembl_grch38_v87_region_annot == "exon" & ensembl_grch38_v92_region_annot %in% c("intron", "intergenic"), 
                                            yes = "exon in v87", 
                                            no = ifelse(ensembl_grch38_v87_region_annot == "exon, intron" & ensembl_grch38_v92_region_annot %in% c("intron", "intergenic"), 
                                                        yes = "exon, intron in v87", 
                                                        no = ifelse(ensembl_grch38_v87_region_annot == "exon, intergenic" & ensembl_grch38_v92_region_annot %in% c("intron", "intergenic"), 
                                                                    yes = "exon, intergenic in v87", 
                                                                    no = NA))))), 
           ensembl_v87_v92_annot = 
             factor(ensembl_v87_v92_annot) %>% 
             fct_relevel(c("exon in v87", "exon, intron in v87", "exon, intergenic in v87", "intergenic in v87", "intron in v87")))
  
  ERs_v87_vs_v92_to_plot <- 
    ERs_v87_vs_v92 %>% 
    filter(!is.na(ensembl_v87_v92_annot)) %>% 
    group_by(tissue, ensembl_v87_v92_annot) %>% 
    summarise(total_kb_ERs = sum(width)/1000) %>% 
    ungroup() %>%
    mutate(total_kb_ERs = ifelse(ensembl_v87_v92_annot %in% c("exon in v87", "exon, intergenic in v87", "exon, intron in v87"), -total_kb_ERs, total_kb_ERs)) %>% 
    left_join(gtex_tissue_name_formatting_w_face %>% dplyr::select(OMIM_gtex_name, gtex_tissues_name_to_plot, gtex_tissue_group, tissue_color_hex, face), by = c("tissue" = "OMIM_gtex_name")) %>%
    mutate(gtex_tissues_name_to_plot = gtex_tissues_name_to_plot %>% str_replace(" $", "") %>% str_c("   ")) 
  
  return(ERs_v87_vs_v92_to_plot)
  
}

plot_num_unannot_vs_annot_v87_v92 <- function(ERs_v87_vs_v92_to_plot, feature_type_order_colors){ 
  
  ERs_v87_vs_v92_to_plot_order <- 
    ERs_v87_vs_v92_to_plot %>% 
    filter(!ensembl_v87_v92_annot %in% c("exon in v87", "exon, intergenic in v87", "exon, intron in v87")) %>% 
    group_by(gtex_tissues_name_to_plot, face, tissue_color_hex, tissue) %>% 
    summarise(intron_inter_total_kb_ERs = sum(total_kb_ERs)) %>% 
    arrange(desc(intron_inter_total_kb_ERs)) %>% 
    ungroup()
  
  ERs_v87_vs_v92_total_kb_ERs <- 
    ggplot(ERs_v87_vs_v92_to_plot %>% 
             mutate(gtex_tissues_name_to_plot = gtex_tissues_name_to_plot %>% 
                      factor() %>% 
                      fct_relevel(ERs_v87_vs_v92_to_plot_order$gtex_tissues_name_to_plot)), 
           aes(x = gtex_tissues_name_to_plot, y = total_kb_ERs)) +
    geom_col(aes(fill = ensembl_v87_v92_annot), colour = "black", width = 0.75) +
    scale_x_discrete(name = "Tissue") +
    scale_y_continuous(name = "Total expressed region length (Kb)", breaks = seq(-25, 125, 25)) + 
    scale_fill_manual(name = "", 
                      values = feature_type_order_colors$feature_type_colour[-4]) +
    theme_pubr(legend = "top") + 
    theme(axis.text.x = 
            element_text(angle = 90, 
                         vjust = 0.5, 
                         hjust = 1, 
                         face = ERs_v87_vs_v92_to_plot_order$face)) +
    geom_point(data = ERs_v87_vs_v92_to_plot_order, 
               aes(x = gtex_tissues_name_to_plot, y = -50), 
               colour = "black", fill = ERs_v87_vs_v92_to_plot_order$tissue_color_hex  %>% str_c("#", .), shape = 22, size = 3.5) +
    coord_cartesian(ylim = c(-30, 130), 
                    clip = "off")
  
  return(ERs_v87_vs_v92_total_kb_ERs)
  
}  

get_legend_v87_v92 <- function(feature_type_order_colors){
  
  feature_type_order_colors_alpha <- 
    feature_type_order_colors %>% 
    mutate(feature_type_order = str_c(feature_type_order, " - alpha"), 
           feature_type_colour = alpha(feature_type_colour, 0.4))
  
  feature_type_order_colors_w_alpha <- 
    bind_rows(feature_type_order_colors, 
              feature_type_order_colors_alpha) %>% 
    mutate(feature_type_order = feature_type_order %>% factor() %>% fct_relevel(feature_type_order))
  
  v87_v92_legend <- 
    (ggplot(feature_type_order_colors_w_alpha, aes(x = feature_type_order, fill = feature_type_order)) +
       geom_bar(colour = "black") +
       scale_fill_manual(name = "", values = feature_type_order_colors_w_alpha$feature_type_colour) +
       guides(fill = guide_legend(ncol = 2)) +
       theme(legend.text = element_blank())) %>% 
    get_legend()
  
  return(v87_v92_legend)
  
}

compare_Mb_total_to_random_inter_intron_regions <- function(ERs_v87_vs_v92_to_plot, intron_inter_randomised_regions_all_tissues_v87_to_v92_paths){
  
  ERs_v87_vs_v92_to_plot_rand <- 
    ERs_v87_vs_v92_to_plot %>% 
    filter(!str_detect(ensembl_v87_v92_annot, "exon")) %>% 
    mutate(ensembl_v87_v92_annot = ensembl_v87_v92_annot %>% as.character() %>% str_replace(" in v87", "")) 
  
  intron_inter_randomised_regions_all_tissues_v87_to_v92 <- data_frame()
  
  for(i in seq_along(intron_inter_randomised_regions_all_tissues_v87_to_v92_paths)){
    
    intron_inter_randomised_regions_all_tissues_v87_to_v92_path <- intron_inter_randomised_regions_all_tissues_v87_to_v92_paths[i]
    
    tissue_to_filter <- intron_inter_randomised_regions_all_tissues_v87_to_v92_path %>% 
      str_replace(".*/", "") %>% 
      str_replace("_intron_inter_exon_in_v92.csv", "")
    
    print(str_c(i, " - ", tissue_to_filter))
    
    suppressMessages(
      
    intron_inter_randomised_regions_per_tissue_v87_to_v92 <- 
      read_delim(intron_inter_randomised_regions_all_tissues_v87_to_v92_path, delim = ",") %>% 
      mutate(iter = rep(sort(rep(1:(10000/100), 100)), 2)) %>% 
      mutate(tissue = tissue_to_filter)
    
    )
    
    # check_distinct_aggregates_per_iter <- 
    #   intron_inter_randomised_regions_per_tissue_v87_to_v92 %>% 
    #   group_by(iter) %>% 
    #   summarise(n_aggregates = n_distinct(query_aggregate))
    
    intron_inter_randomised_regions_all_tissues_v87_to_v92 <- 
      intron_inter_randomised_regions_all_tissues_v87_to_v92 %>% 
      bind_rows(intron_inter_randomised_regions_per_tissue_v87_to_v92)
    
  }
  
  stopifnot(all(ERs_v87_vs_v92_to_plot_rand$tissue %in% unique(intron_inter_randomised_regions_all_tissues_v87_to_v92$tissue)))
  
  tissues_to_filter <- ERs_v87_vs_v92_to_plot_rand$tissue %>% unique()
  
  ens_v87_93_wilcox_test_pval_all_tissues <- data_frame()
  
  for(i in seq_along(tissues_to_filter)){
    
    tissue_to_filter <- tissues_to_filter[i]
    
    print(str_c(i, " - ", tissue_to_filter))
    
    ens_v87_93_wilcox_test_pval <- 
      data_frame(tissue = tissue_to_filter, 
                 ensembl_v87_v92_annot = c("intergenic", "intron"), 
                 wilcox_pvalue = as.numeric(NA))
    
    for(j in 1:nrow(ens_v87_93_wilcox_test_pval)){
      
      ensembl_v87_v92_annot_to_filter <- ens_v87_93_wilcox_test_pval$ensembl_v87_v92_annot[j]
      
      total_kb_ERs_mu <- 
        ERs_v87_vs_v92_to_plot_rand %>% filter(tissue == tissue_to_filter, 
                                             ensembl_v87_v92_annot == ensembl_v87_v92_annot_to_filter) %>% 
        .[["total_kb_ERs"]]
      
      total_kb_random_regions <- 
        intron_inter_randomised_regions_all_tissues_v87_to_v92 %>% 
         filter(tissue == tissue_to_filter, theRegion == ensembl_v87_v92_annot_to_filter) %>% 
      .[["total_kb"]]
      
      wilcox_results <- wilcox.test(x = total_kb_random_regions, mu = total_kb_ERs_mu)
      
      ens_v87_93_wilcox_test_pval$wilcox_pvalue[j] <- wilcox_results$p.value
        
    }
    
    ens_v87_93_wilcox_test_pval_all_tissues <- 
      ens_v87_93_wilcox_test_pval_all_tissues %>% 
      bind_rows(ens_v87_93_wilcox_test_pval)
    
  }
    
    
  
  intron_inter_randomised_regions_all_tissues_v87_to_v92_w_total_kb_ERs <- 
    intron_inter_randomised_regions_all_tissues_v87_to_v92 %>% 
    left_join(ERs_v87_vs_v92_to_plot_rand, by = c("tissue" = "tissue", "theRegion" = "ensembl_v87_v92_annot")) %>% 
    left_join(ens_v87_93_wilcox_test_pval_all_tissues, by = c("tissue" = "tissue", "theRegion" = "ensembl_v87_v92_annot"))
  
  return(intron_inter_randomised_regions_all_tissues_v87_to_v92_w_total_kb_ERs)
  
  
}

plot_num_unannot_vs_annot_v87_v92_vs_rand_regions <- function(ERs_v87_vs_v92_to_plot, v87_to_v92_random_vs_ERs_tidy){
  
  ERs_v87_vs_v92_to_plot_order <- 
    ERs_v87_vs_v92_to_plot %>% 
    filter(!ensembl_v87_v92_annot %in% c("exon in v87", "exon, intergenic in v87", "exon, intron in v87")) %>% 
    group_by(gtex_tissues_name_to_plot, face, tissue_color_hex, tissue) %>% 
    summarise(intron_inter_total_kb_ERs = sum(total_kb_ERs)) %>% 
    arrange(desc(intron_inter_total_kb_ERs)) %>% 
    ungroup()
  
  ERs_rand_v87_vs_v92_all_tissues <- 
    ggplot(v87_to_v92_random_vs_ERs_tidy %>% 
           mutate(gtex_tissues_name_to_plot = gtex_tissues_name_to_plot %>% factor() %>% fct_relevel(ERs_v87_vs_v92_to_plot_order$gtex_tissues_name_to_plot)), 
         aes(x = gtex_tissues_name_to_plot, y = total_kb)) +
    geom_boxplot(size = rel(0.25)) +
    geom_point(aes(y = total_kb_ERs), colour = "#CD534CFF", shape = 4, size = rel(1.5)) +
    scale_y_continuous("Total Kb becoming to exonic in v92") + 
    scale_x_discrete(name = "Tissue") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, 
                                     face = ERs_v87_vs_v92_to_plot_order$face)) +
    facet_wrap(theRegion ~ ., nrow = 2, ncol = 1) 
  
  ggsave(plot = ERs_rand_v87_vs_v92_all_tissues, filename = "ERs_rand_v87_vs_v92_all_tissues.png", path = "OMIM_paper/supp_figures/rand_vs_ERs_v87_to_v92/", 
         width = 8.27, height = (11.69/1.3), dpi = 600, units = "in")

  return(ERs_rand_v87_vs_v92_all_tissues)
  

}

# Main ------------------------------------------------------------------------------------------------

intron_inter_randomised_regions_all_tissues_v87_to_v92_paths <- 
  list.files("results/generate_randomised_intron_inter_regions/intron_inter_randomised_regions_all_tissues_v87_to_v92/", full.names = T)

feature_type_order_colors <- 
  data_frame(feature_type_order = c("exon", "exon, intron", "exon, intergenic", 
                                    "exon, intergenic, intron", 
                                    "intergenic", "intron"), 
             feature_type_colour = get_palette("jco", 9)[c(6, 1, 5, 3, 9, 4)])


ERs_w_annotation_all_tissues <- 
  ERs_w_annotation_all_tissues %>%
  filter(width > 3, !str_detect(tissue, "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast"),
         ensembl_grch38_v92_region_annot != "exon, intergenic, intron")

ERs_v87_vs_v92_to_plot <- get_ensembl_v87_to_v92_tidy(ERs_w_annotation_all_tissues, feature_type_order_colors, gtex_tissue_name_formatting)

ERs_v87_vs_v92_n_total_kb_ERs <- 
  plot_num_unannot_vs_annot_v87_v92(ERs_v87_vs_v92_to_plot, feature_type_order_colors)

legend_v87_v92 <- get_legend_v87_v92(feature_type_order_colors %>% filter(feature_type_order != c("exon, intergenic, intron")))

v87_to_v92_random_vs_ERs_tidy <- 
  compare_Mb_total_to_random_inter_intron_regions(ERs_v87_vs_v92_to_plot, intron_inter_randomised_regions_all_tissues_v87_to_v92_paths)

ERs_rand_v87_vs_v92_all_tissues <- 
  plot_num_unannot_vs_annot_v87_v92_vs_rand_regions(ERs_v87_vs_v92_to_plot, v87_to_v92_random_vs_ERs_tidy)

# Save data -------------------------------------------------------------------------------------------

ggsave(plot = ERs_v87_vs_v92_n_total_kb_ERs, filename = "ERs_v87_vs_v92_kb.png", path = "OMIM_paper/figures/validation_of_ERs/",
       width = 8.27, height = (11.69/2), dpi = 600, units = "in")

ggsave(plot = legend_v87_v92, filename = "legend_v87_v92.png", path = "OMIM_paper/figures/validation_of_ERs/", 
       width = 1.27, height = (11.69/3))

ggsave(plot = ERs_rand_v87_vs_v92_all_tissues, filename = "ERs_rand_v87_vs_v92_all_tissues.png", path = "OMIM_paper/supp_figures/rand_vs_ERs_v87_to_v92/", 
       width = 8.27, height = (11.69/1.3), dpi = 600, units = "in")

