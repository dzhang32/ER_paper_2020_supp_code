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

plot_num_propor_total_Mb_ERs_by_annot_features <- 
  function(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific, refer_annot_name = "ensembl_grch38_v92", 
           feature_type_order_colors, gtex_tissue_name_formatting){
  
  gtex_tissue_name_formatting_w_colours <- 
  gtex_tissue_name_formatting %>% 
    mutate(face = ifelse(gtex_tissue_group == "brain", "bold", "plain"))
  
  # tidy data to plot 
  ERs_w_annotation_all_tissues_tidy  <- 
    ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific %>% 
    dplyr::select(contains("region_annot"), width, tissue, contains("overlap_")) %>% 
    group_by(tissue) %>% 
    mutate(n_total_ERs = n()) %>% 
    gather(key = "refer_annot_ver", value = "feature_type", contains("ensembl")) %>% 
    left_join(gtex_tissue_name_formatting %>% dplyr::select(OMIM_gtex_name, gtex_tissues_name_to_plot, gtex_tissue_group), 
              by = c("tissue" = "OMIM_gtex_name"))  %>% 
    mutate(refer_annot_ver = refer_annot_ver %>% str_replace("_region_annot", ""), 
           gtex_tissues_name_to_plot_w_space = gtex_tissues_name_to_plot %>% str_c(., "   "))
  
  ERs_w_annotation_all_tissues_tidy_to_plot <- 
    ERs_w_annotation_all_tissues_tidy %>% 
    group_by(tissue, gtex_tissues_name_to_plot, gtex_tissues_name_to_plot_w_space, n_total_ERs, refer_annot_ver, feature_type) %>% 
    summarise(n_feature_type = n(), 
              total_Mb = sum(width)/1000000, 
              mean_width = mean(width)) %>% 
    ungroup() %>% 
    mutate(propor_feature_type = n_feature_type/n_total_ERs, 
           feature_type = feature_type %>% factor() %>% fct_relevel(feature_type_order_colors$feature_type_order))
  
  # get the order of tissues in terms of descending number of intron/intergenic regions 
  # desc_intron_inter_tissue_order_propor_intron_inter <- 
  #   ERs_w_annotation_all_tissues_tidy_to_plot %>% 
  #   filter(refer_annot_ver == refer_annot_name, 
  #          feature_type == "intron" | feature_type == "intergenic") %>% 
  #   group_by(gtex_tissues_name_to_plot, refer_annot_ver) %>% 
  #   summarise(sum_propor_intron_inter = sum(propor_feature_type)) %>% 
  #   arrange(desc(sum_propor_intron_inter))
  # 
  # propor_ERs_feature_type_plot <- 
  #   ggplot(ERs_w_annotation_all_tissues_tidy_to_plot %>% 
  #          filter(refer_annot_ver == refer_annot_name) %>% 
  #          mutate(gtex_tissues_name_to_plot = gtex_tissues_name_to_plot %>% 
  #                   factor() %>% fct_relevel(desc_intron_inter_tissue_order_propor_intron_inter[["gtex_tissues_name_to_plot"]])), 
  #        aes(x = gtex_tissues_name_to_plot, y = propor_feature_type)) +
  #   geom_col(aes(fill = feature_type), colour = "black") +
  #   scale_x_discrete(name = "Tissue") +
  #   scale_y_continuous(name = "Proportion of expressed regions") + 
  #   scale_fill_manual(name = "Annotation feature", 
  #                     values = feature_type_order_colors$feature_type_colour, 
  #                     guide = F) +
  #   theme_pubr(legend = "right") + 
  #   theme(axis.text.x = 
  #           element_text(angle = 90, 
  #                        vjust = 0.5, 
  #                        hjust=1))
  
  # get the order of tissues in terms of descending number of intron/intergenic regions 
  desc_intron_inter_tissue_order_total_Mb <- 
    ERs_w_annotation_all_tissues_tidy_to_plot %>% 
    filter(refer_annot_ver == refer_annot_name,  
           feature_type == "intron" | feature_type == "intergenic") %>% 
    group_by(gtex_tissues_name_to_plot, gtex_tissues_name_to_plot_w_space, refer_annot_ver) %>% 
    summarise(sum_total_intron_intergenic_Mb = sum(total_Mb)) %>% 
    arrange(desc(sum_total_intron_intergenic_Mb)) %>% 
    left_join(gtex_tissue_name_formatting_w_colours %>% dplyr::select(gtex_tissues_name_to_plot, tissue_color_hex, face))
  
  total_Mb_ERs_feature_type_plot <- 
    ggplot(ERs_w_annotation_all_tissues_tidy_to_plot %>% 
           filter(refer_annot_ver == refer_annot_name) %>% 
           mutate(gtex_tissues_name_to_plot_w_space = gtex_tissues_name_to_plot_w_space %>% factor() %>% fct_relevel(desc_intron_inter_tissue_order_total_Mb[["gtex_tissues_name_to_plot_w_space"]])), 
         aes(x = gtex_tissues_name_to_plot_w_space, y = total_Mb)) +
    geom_col(aes(fill = feature_type), colour = "black") +
    scale_x_discrete(name = "Tissue") +
    scale_y_continuous(name = "Total expressed region length (Mb)") + 
    scale_fill_manual(name = "Annotation feature:", 
                      values = feature_type_order_colors$feature_type_colour) +
    # facet_wrap(~feature_type, scales = "free") +
    theme_pubr(legend = "right") + 
    theme(axis.text.x = 
            element_text(angle = 90, 
                         vjust = 0.5, 
                         hjust=1, 
                         face = desc_intron_inter_tissue_order_total_Mb$face, 
                         size = rel(0.6)), 
          legend.position = "top")
  
  total_Mb_ERs_all_feature_type_plot_list <- list()
  
  for(i in seq_along(feature_type_order_colors$feature_type_order)){
    
    feature_type_to_filter <- feature_type_order_colors$feature_type_order[i]
    feature_type_color <- feature_type_order_colors$feature_type_colour[i]
    
    total_Mb_ERs_one_feature_type_plot <- 
      ggplot(ERs_w_annotation_all_tissues_tidy_to_plot %>% 
               filter(refer_annot_ver == refer_annot_name, 
                      feature_type == feature_type_to_filter) %>% 
               mutate(gtex_tissues_name_to_plot_w_space = gtex_tissues_name_to_plot_w_space %>% factor() %>% fct_relevel(desc_intron_inter_tissue_order_total_Mb[["gtex_tissues_name_to_plot_w_space"]])), 
             aes(x = gtex_tissues_name_to_plot_w_space, y = total_Mb)) +
      geom_col(fill = feature_type_color, colour = "black") +
      scale_x_discrete(name = "Tissue") +
      scale_y_continuous(name = "Total expressed region length (Mb)") + 
      theme_pubr(legend = "right") + 
      theme(axis.text.x = 
              element_text(angle = 90, 
                           vjust = 0.5, 
                           hjust=1, 
                           face = desc_intron_inter_tissue_order_total_Mb$face, 
                           size = rel(0.6))) 
    
    if(feature_type_to_filter == "intergenic"){
      
      total_Mb_ERs_one_feature_type_plot <- 
        total_Mb_ERs_one_feature_type_plot +
        geom_point(data = desc_intron_inter_tissue_order_total_Mb, 
                   aes(x = gtex_tissues_name_to_plot_w_space, y = -0.9), 
                   colour = "black", fill = desc_intron_inter_tissue_order_total_Mb$tissue_color_hex  %>% str_c("#", .), shape = 22, size = 2) +
        coord_cartesian(ylim = c(0, 7), 
                        clip = "off")
        
    }else if(feature_type_to_filter == "intron"){
      
      total_Mb_ERs_one_feature_type_plot <- 
        total_Mb_ERs_one_feature_type_plot +
        geom_point(data = desc_intron_inter_tissue_order_total_Mb, 
                   aes(x = gtex_tissues_name_to_plot_w_space, y = -1.8), 
                   colour = "black", fill = desc_intron_inter_tissue_order_total_Mb$tissue_color_hex  %>% str_c("#", .), shape = 22, size = 2) +
        coord_cartesian(ylim = c(0, 14), 
                        clip = "off")
      
    }
    
    total_Mb_ERs_all_feature_type_plot_list[[i]] <- total_Mb_ERs_one_feature_type_plot
    
  }
  
  total_Mb_ERs_all_feature_type_plot_list[[i+1]] <- total_Mb_ERs_feature_type_plot
    
  names(total_Mb_ERs_all_feature_type_plot_list) <- 
    c(feature_type_order_colors$feature_type_order, "all")

  # generate table with the propor and total Mb across feature types including sum of intron and intergenic
  # total_Mb_across_diff_feature_types <- 
  #   ERs_w_annotation_all_tissues_tidy_to_plot %>% 
  #   filter(refer_annot_ver == refer_annot_name, 
  #          feature_type %in% c("exon", "intron", "intergenic")) %>% 
  #   group_by(tissue, refer_annot_ver, feature_type) %>% 
  #   summarise(sum_total_Mb= sum(total_Mb)) %>% 
  #   ungroup() %>% 
  #   mutate(feature_type = feature_type %>% str_c(., "_total_Mb")) %>% 
  #   spread(key = "feature_type", value = "sum_total_Mb") %>% 
  #   mutate(intron_intergenic_total_Mb = intergenic_total_Mb + intron_total_Mb, 
  #          intron_intergenic_ratio = intron_total_Mb/intergenic_total_Mb) %>% 
  #   dplyr::select(-refer_annot_ver) %>% 
  #   arrange(desc(intron_intergenic_total_Mb))
  # 
  # total_Mb_across_diff_feature_types_ggtable <- 
  #   bind_cols(total_Mb_across_diff_feature_types[1], total_Mb_across_diff_feature_types[2:length(total_Mb_across_diff_feature_types)] %>% round(digits = 2)) %>% 
  #   ggtexttable(rows = NULL, theme = ttheme("mBlue"))
  
  return(total_Mb_ERs_all_feature_type_plot_list)
  
}

plot_brain_vs_non_brain_inter_intron_Mb <- function(){
  
  ERs_w_annotation_all_tissues_tidy_to_plot %>% 
    filter(refer_annot_ver == "ensembl_grch38_v87", feature_type %in% c("intergenic", "intron")) 
    
  ERs_w_annotation_all_tissues_tidy_to_plot <- 
    ERs_w_annotation_all_tissues_tidy %>% 
    group_by(gtex_tissue_group, tissue, gtex_tissues_name_to_plot, n_total_ERs, refer_annot_ver, feature_type) %>% 
    summarise(n_feature_type = n(), 
              total_Mb = sum(width)/1000000, 
              mean_width = mean(width)) %>% 
    ungroup() %>% 
    mutate(propor_feature_type = n_feature_type/n_total_ERs, 
           feature_type = feature_type %>% factor() %>% fct_relevel(feature_type_order_colors$feature_type_order))
  
  ERs_w_annotation_all_tissues_tidy  <- 
    ERs_w_annotation_all_tissues %>% 
    dplyr::select(contains("region_annot"), width, tissue, contains("overlap_"), 
                  annotationType_split_read_annot) %>% 
    group_by(tissue) %>% 
    mutate(n_total_ERs = n()) %>% 
    gather(key = "refer_annot_ver", value = "feature_type", contains("ensembl")) %>% 
    mutate(refer_annot_ver = refer_annot_ver %>% str_replace("_region_annot", "")) %>% 
    left_join(gtex_tissue_name_formatting %>% dplyr::select(OMIM_gtex_name, gtex_tissues_name_to_plot, gtex_tissue_group), by = c("tissue" = "OMIM_gtex_name"))
  
  total_Mb_annotated_by_split_reads_across_feature_type <- 
    ERs_w_annotation_all_tissues_tidy %>% 
    filter(refer_annot_ver == refer_annot_name, !is.na(annotationType_split_read_annot)) %>% 
    group_by(gtex_tissue_group, tissue, gtex_tissues_name_to_plot, refer_annot_ver, n_total_ERs, feature_type) %>% 
    summarise(total_Mb = sum(width)/1000000, 
              n_ERs = n()) %>% 
    ungroup() %>% 
    mutate(propor_ERs = (n_ERs/n_total_ERs))
  
    desc_intron_inter_tissue_order_total_Mb <- 
      total_Mb_annotated_by_split_reads_across_feature_type %>% 
      group_by(gtex_tissues_name_to_plot, gtex_tissue_group, refer_annot_ver) %>% 
      mutate(total_tissue_Mb = sum(total_Mb)) %>% 
      ungroup() %>% 
      filter(refer_annot_ver == refer_annot_name,  
             feature_type %in% c("exon, intergenic") ) %>% 
      group_by(gtex_tissues_name_to_plot, gtex_tissue_group, refer_annot_ver) %>% 
      summarise(sum_total_intron_intergenic_Mb = sum(total_Mb), 
                propor_total_intron_intergenic_Mb = (sum(total_Mb)/(total_tissue_Mb)))
    
    desc_intron_inter_tissue_order_total_Mb %>% 
      filter(gtex_tissue_group == "brain", gtex_tissues_name_to_plot != "Cerebellum", 
             gtex_tissues_name_to_plot != "Cortex") %>% .[["propor_total_intron_intergenic_Mb"]] %>% mean()
    
    desc_intron_inter_tissue_order_total_Mb %>% 
      filter(gtex_tissue_group != "brain") %>% .[["propor_total_intron_intergenic_Mb"]] %>% mean()
    
    brain_inter_intron_Mb <- 
      desc_intron_inter_tissue_order_total_Mb %>% 
      filter(gtex_tissue_group == "brain", gtex_tissues_name_to_plot != "Cerebellum", 
             gtex_tissues_name_to_plot != "Cortex") %>% .[["propor_total_intron_intergenic_Mb"]]
    
    non_brain_inter_intron_Mb <- 
      desc_intron_inter_tissue_order_total_Mb %>% 
      filter(gtex_tissue_group != "brain") %>% .[["propor_total_intron_intergenic_Mb"]]
    
    t.test(brain_inter_intron_Mb, non_brain_inter_intron_Mb, alternative="greater")
  
}

compare_standardised_variance_exon_intron_intergenic <- function(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific){
  
  ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_total_Mb <- 
  ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific %>% 
    group_by(tissue, ensembl_grch38_v92_region_annot) %>%  
    summarise(total_Mb = sum(width)/1000000)
  
  ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_total_Mb %>% 
    group_by(ensembl_grch38_v92_region_annot) %>% 
    summarise(standardised_var = sd(total_Mb)/mean(total_Mb))
  
  
}

get_novel_exon_intron <- function(ERs_w_annotation_all_tissues){
  
  ERs_w_annotation_all_tissues_exon_intron <- 
  ERs_w_annotation_all_tissues %>% 
    filter(ensembl_grch38_v87_region_annot == "exon, intron") 
  
  ERs_w_annotation_all_tissues_exon_intron_gr <- 
    GRanges(str_c(ERs_w_annotation_all_tissues_exon_intron$seqnames, ":", 
                  ERs_w_annotation_all_tissues_exon_intron$start, "-", 
                  ERs_w_annotation_all_tissues_exon_intron$end), 
            tissue = ERs_w_annotation_all_tissues_exon_intron$tissue)
  
  genomic_state_ensembl_grch38_v87 <-
    generate_genomic_state(ensembl_grch38_v87_TxDb,
                           output_path = "/data/references/ensembl/genomic_state/v87/ensembl_grch38_v87_genomic_state.rda",
                           chrs_to_keep = str_c("chr", c(1:22, "X", "Y", "M")))
  
  genomic_state_ensembl_grch38_v87_introns <- 
    genomic_state_ensembl_grch38_v87$fullGenome[genomic_state_ensembl_grch38_v87$fullGenome$theRegion == "intron"]
  
  exon_intron_novel_tissues_all <- GRanges()
  
  for(i in ERs_w_annotation_all_tissues_exon_intron_gr$tissue %>% unique()){
    
    tissue <- i
    print(tissue)
    ERs_w_annotation_all_tissues_exon_intron_gr_tissue_filtered <- ERs_w_annotation_all_tissues_exon_intron_gr[ERs_w_annotation_all_tissues_exon_intron_gr$tissue == tissue]
    
    
    exon_intron_novel_tissue_filtered <- 
      intersect(ERs_w_annotation_all_tissues_exon_intron_gr_tissue_filtered, 
                genomic_state_ensembl_grch38_v87_introns, ignore.strand=T)
    exon_intron_novel_tissue_filtered$tissue <- tissue
    
    exon_intron_novel_tissues_all <- c(exon_intron_novel_tissues_all, exon_intron_novel_tissue_filtered)
    
  }
  
  to_plot <- 
    exon_intron_novel_tissues_all %>% as.data.frame() %>% 
    left_join(gtex_tissue_name_formatting %>% dplyr::select(OMIM_gtex_name, gtex_tissues_name_to_plot, gtex_tissue_group), by = c("tissue" = "OMIM_gtex_name")) %>% 
    group_by(gtex_tissues_name_to_plot) %>% 
    summarise(total_Mb = sum(width)/1000000)
  
  total_Mb_ERs_feature_type_plot <- 
    ggplot(to_plot%>% 
             mutate(gtex_tissues_name_to_plot = gtex_tissues_name_to_plot %>% factor() %>% fct_relevel(desc_intron_inter_tissue_order_total_Mb[["gtex_tissues_name_to_plot"]])), 
           aes(x = gtex_tissues_name_to_plot, y = total_Mb)) +
    geom_col(fill = "#0073C2FF", colour = "black") +
    scale_x_discrete(name = "Tissue") +
    scale_y_continuous(name = "Total expressed regions length (Mb)") + 
    scale_fill_manual(name = "Annotation feature", 
                      values = feature_type_order_colors$feature_type_colour) +
    # facet_wrap(~feature_type, scales = "free") +
    theme_pubr(legend = "right") + 
    theme(axis.text.x = 
            element_text(angle = 90, 
                         vjust = 0.5, 
                         hjust=1, 
                         face = desc_intron_inter_tissue_order_total_Mb$face)) 
  
  
}

plot_num_propor_total_Mb_ERs_annotated_by_split_reads <- function(ERs_w_annotation_all_tissues, refer_annot_name = "ensembl_grch38_v87", feature_type_order_colors){
  
  gtex_tissue_name_formatting_w_colours <- 
    gtex_tissue_name_formatting %>% 
    mutate(face = ifelse(gtex_tissue_group == "brain", "bold", "plain"))
  
  ERs_w_annotation_all_tissues_tidy  <- 
    ERs_w_annotation_all_tissues %>% 
    dplyr::select(contains("region_annot"), width, tissue, contains("overlap_"), 
                  annotationType_split_read_annot) %>% 
    group_by(tissue) %>% 
    mutate(n_total_ERs = n()) %>% 
    gather(key = "refer_annot_ver", value = "feature_type", contains("ensembl")) %>% 
    mutate(refer_annot_ver = refer_annot_ver %>% str_replace("_region_annot", "")) %>% 
    left_join(gtex_tissue_name_formatting %>% dplyr::select(OMIM_gtex_name, gtex_tissues_name_to_plot, gtex_tissue_group), by = c("tissue" = "OMIM_gtex_name"))
  
  total_Mb_annotated_by_split_reads_across_feature_type <- 
    ERs_w_annotation_all_tissues_tidy %>% 
    filter(refer_annot_ver == refer_annot_name, !is.na(annotationType_split_read_annot)) %>% 
    group_by(tissue, gtex_tissues_name_to_plot, refer_annot_ver, n_total_ERs, feature_type) %>% 
    summarise(total_Mb = sum(width)/1000000, 
              n_ERs = n()) %>% 
    ungroup() %>% 
    mutate(propor_ERs = (n_ERs/n_total_ERs))

  
  # get the order of tissues in terms of descending number of intron/intergenic regions 
  desc_intron_inter_tissue_order_total_Mb <- 
    ERs_w_annotation_all_tissues_tidy_to_plot %>% 
    filter(refer_annot_ver == refer_annot_name,  
           feature_type == "intron" | feature_type == "intergenic") %>% 
    group_by(gtex_tissues_name_to_plot, refer_annot_ver) %>% 
    summarise(sum_total_intron_intergenic_Mb = sum(total_Mb)) %>% 
    arrange(desc(sum_total_intron_intergenic_Mb)) %>% 
    left_join(gtex_tissue_name_formatting_w_colours %>% dplyr::select(gtex_tissues_name_to_plot, face))
  
  total_Mb_ERs_feature_type_plot <- 
    ggplot(total_Mb_annotated_by_split_reads_across_feature_type %>% 
             filter(refer_annot_ver == refer_annot_name, !str_detect(tissue, "cells|testis|vagina|ovary|uterus|prostate"))  %>% 
             mutate(gtex_tissues_name_to_plot = gtex_tissues_name_to_plot %>% factor() %>% fct_relevel(desc_intron_inter_tissue_order_total_Mb[["gtex_tissues_name_to_plot"]]), 
                    feature_type = feature_type %>% factor() %>% fct_relevel(total_Mb_annotated_by_split_reads_across_feature_type$feature_type %>% unique() %>% .[c(1, 4, 2, 3, 5, 6)])), 
           aes(x = gtex_tissues_name_to_plot, y = total_Mb)) +
    geom_col(aes(fill = feature_type), colour = "black") +
    scale_x_discrete(name = "Tissue") +
    scale_y_continuous(name = "Total expressed regions length (Mb)") + 
    scale_fill_manual(name = "Annotation feature", 
                      values = feature_type_order_colors$feature_type_colour) +
    facet_wrap(~feature_type, scales = "free") +
    theme_pubr(legend = "right") + 
    theme(axis.text.x = 
            element_text(angle = 90, 
                         vjust = 0.5, 
                         hjust=1, 
                         face = desc_intron_inter_tissue_order_total_Mb$face)) 
  
  ggsave(plot = total_Mb_ERs_feature_type_plot, filename = "total_Mb_ERs_feature_type_plot_facet_feature_type.png", 
         path = "ESHG_2018/", dpi = 600, width = 28, height = 14)
  
  propor_total_Mb_ERs_feature_type_arranged <- 
    ggarrange(plotlist = list(propor_ERs_split_read_feature_type_plot, total_Mb_ERs_split_read_feature_type_plot), 
              ncol = 2, nrow = 1, legend = "right", widths = c(1.6, 2)) 
  

  
  # generate table with the propor and total Mb across feature types including sum of intron and intergenic
  total_Mb_across_diff_feature_types <- 
    total_Mb_annotated_by_split_reads_across_feature_type %>% 
    filter(refer_annot_ver == refer_annot_name, 
           feature_type %in% c("exon", "intron", "intergenic")) %>%
    group_by(tissue, refer_annot_ver, feature_type) %>% 
    summarise(sum_total_Mb= sum(total_Mb)) %>% 
    ungroup() %>% 
    mutate(feature_type = feature_type %>% str_c(., "_total_Mb")) %>% 
    spread(key = "feature_type", value = "sum_total_Mb") %>% 
    mutate(intron_intergenic_total_Mb = intergenic_total_Mb + intron_total_Mb, 
           intron_intergenic_ratio = intron_total_Mb/intergenic_total_Mb) %>% 
    dplyr::select(-refer_annot_ver) %>% 
    arrange(desc(intron_intergenic_total_Mb))
  
  total_Mb_across_diff_feature_types_ggtable <- 
    bind_cols(total_Mb_across_diff_feature_types[1], total_Mb_across_diff_feature_types[2:length(total_Mb_across_diff_feature_types)] %>% round(digits = 2)) %>% 
    ggtexttable(rows = NULL, theme = ttheme("mBlue"))
  
  return(list(propor_total_Mb_ERs_feature_type_arranged, total_Mb_across_diff_feature_types_ggtable))
  
  
}

plot_width_mean_cov_for_validated_or_split_read_supported_ERs <- function(ERs_w_annotation_all_tissues){
  
  hipp_regions_seb_optimal <- region_cuts[["10.1"]]
  hipp_regions_dz_optimal <- 
    list_ERs_w_annotation[["2 - hippocampus" ]] %>% 
    keepSeqlevels(hipp_regions_seb_optimal %>% seqlevels(), pruning.mode = "coarse")
  
  overlap_eq <- findOverlaps(hipp_regions_dz_optimal, hipp_regions_seb_optimal, type = "equal")
  overlap_any <- findOverlaps(hipp_regions_dz_optimal, hipp_regions_seb_optimal, type = "any", maxgap = -1, minoverlap = 1)

  hipp_regions_dz_optimal_w_seb_valid <- 
    hipp_regions_dz_optimal %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(validated_in_hipp_seb_equal_any = ifelse(row_number() %in% queryHits(overlap_eq), "validated - precise boundaries", 
                                                    ifelse(row_number() %in% queryHits(overlap_any), "validated - any", "not validated"))) %>% 
    mutate(annotationType_split_read_annot_fct = 
             ifelse(is.na(annotationType_split_read_annot), "none", 
                    ifelse(annotationType_split_read_annot == 1, "(partially) annotated split read", "unannotated split read")) %>% 
             factor() %>% 
             fct_relevel("(partially) annotated split read", "unannotated split read", "none"),
           ensembl_v87_v92_annot = 
             ifelse(ensembl_grch38_v87_region_annot == "intron" & ensembl_grch38_v92_region_annot %in% c("exon", "exon, intron", "exon, intergenic"), 
                    yes = "intron - annotated in v92", 
                    no = ifelse(ensembl_grch38_v87_region_annot == "intergenic" & ensembl_grch38_v92_region_annot %in% c("exon", "exon, intron", "exon, intergenic"), 
                                yes = "intergenic - annotated in v92", 
                                no = ifelse(ensembl_grch38_v87_region_annot == "exon" & ensembl_grch38_v92_region_annot %in% c("intron", "intergenic"), 
                                            yes = "exon - unannotated in v92", 
                                            no = ifelse(ensembl_grch38_v87_region_annot == "exon, intron" & ensembl_grch38_v92_region_annot %in% c("intron", "intergenic"), 
                                                        yes = "exon, intron - unannotated in v92", 
                                                        no = ifelse(ensembl_grch38_v87_region_annot == "exon, intergenic" & ensembl_grch38_v92_region_annot %in% c("intron", "intergenic"), 
                                                                    yes = "exon, intergenic - unannotated in v92", 
                                                                    no = NA))))),
           ensembl_v87_v92_annot = 
             factor(ensembl_v87_v92_annot) %>% 
             fct_relevel(c("exon - unannotated in v92", "exon, intron - unannotated in v92", "exon, intergenic - unannotated in v92", 
                           "intergenic - annotated in v92", "intron - annotated in v92")))
  
  ggplot(hipp_regions_dz_optimal_w_seb_valid %>% 
           mutate(ensembl_grch38_v87_region_annot = ensembl_grch38_v87_region_annot %>% factor() %>% fct_relevel(feature_type_order_colors$feature_type_order)), 
         aes(x = width, y = value)) + 
    geom_point(aes(colour = ensembl_grch38_v87_region_annot), size = 1) +
    scale_x_log10(name = "Width", breaks = c(0, 10, 100, 1000, 10000)) +
    scale_y_log10(name = "Coverage", breaks = c(0, 10, 100, 1000, 10000)) + 
    scale_colour_manual(name = "Annotation feature", 
                      values = feature_type_order_colors$feature_type_colour) +
    facet_grid(validated_in_hipp_seb_equal_any ~ ensembl_grch38_v87_region_annot) +
    theme_bw()
  
  ggplot(hipp_regions_dz_optimal_w_seb_valid %>% 
           mutate(ensembl_grch38_v87_region_annot = ensembl_grch38_v87_region_annot %>% factor() %>% fct_relevel(feature_type_order_colors$feature_type_order)), 
         aes(x = width, y = value)) + 
    geom_point(aes(colour = ensembl_grch38_v87_region_annot), size = 1) +
    scale_x_log10(name = "Width", breaks = c(0, 10, 100, 1000, 10000)) +
    scale_y_log10(name = "Coverage", breaks = c(0, 10, 100, 1000, 10000)) + 
    scale_colour_manual(name = "Annotation feature", 
                        values = feature_type_order_colors$feature_type_colour) +
    facet_grid(annotationType_split_read_annot_fct ~ ensembl_grch38_v87_region_annot) +
    theme_bw()
  
  ggplot(hipp_regions_dz_optimal_w_seb_valid %>% 
           mutate(ensembl_grch38_v87_region_annot = ensembl_grch38_v87_region_annot %>% factor() %>% fct_relevel(feature_type_order_colors$feature_type_order)), 
         aes(x = width, y = value)) + 
    geom_point(aes(colour = ensembl_grch38_v87_region_annot), size = 1) +
    scale_x_log10(name = "Width", breaks = c(0, 10, 100, 1000, 10000)) +
    scale_y_log10(name = "Coverage", breaks = c(0, 10, 100, 1000, 10000)) + 
    scale_colour_manual(name = "Annotation feature", 
                        values = feature_type_order_colors$feature_type_colour) +
    facet_grid(annotationType_split_read_annot_fct ~ ensembl_v87_v92_annot) +
    theme_bw()
 
  
}

##### Second level ##### 

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_txDb_from_gtf.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_genomic_state.R")


# Main ------------------------------------------------------------------------------------------------

feature_type_order_colors <- 
  data_frame(feature_type_order = c("exon", "exon, intron", "exon, intergenic", 
                                    "exon, intergenic, intron", 
                                    "intergenic", "intron"), 
             feature_type_colour = get_palette("jco", 9)[c(6, 1, 5, 3, 9, 4)])

ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific <- 
  ERs_w_annotation_all_tissues %>% 
  filter(width > 3, !str_detect(tissue, "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast"), 
         ensembl_grch38_v92_region_annot != "exon, intergenic, intron")

propor_total_Mb_ERs_by_annot_features <- 
  plot_num_propor_total_Mb_ERs_by_annot_features(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific, 
                                                 refer_annot_name = "ensembl_grch38_v92", 
                                                 feature_type_order_colors %>% filter(feature_type_order != "exon, intergenic, intron"), 
                                                 gtex_tissue_name_formatting)

propor_total_Mb_ERs_by_annot_features_arranged <- 
  ggarrange(
  
  ggarrange(propor_total_Mb_ERs_by_annot_features$all + rremove("xlab") + rremove("x.text") + rremove("x.ticks") + rremove("ylab") + rremove("legend"), 
            propor_total_Mb_ERs_by_annot_features$`exon, intron` + rremove("xlab") + rremove("x.text") + rremove("x.ticks"), 
            propor_total_Mb_ERs_by_annot_features$intergenic + rremove("ylab"), 
            ncol = 1, nrow = 3, align = "v", heights = c(0.6, 0.6, 1)),
  
  ggarrange(propor_total_Mb_ERs_by_annot_features$exon + rremove("xlab") + rremove("x.text") + rremove("x.ticks") + rremove("ylab") ,
            propor_total_Mb_ERs_by_annot_features$`exon, intergenic` + rremove("xlab") + rremove("x.text") + rremove("x.ticks"),
            propor_total_Mb_ERs_by_annot_features$intron + rremove("ylab"),
            ncol = 1, nrow = 3, align = "v", heights = c(0.6, 0.6, 1)), 
  ncol = 2, nrow = 1
  
)

propor_total_Mb_ERs_by_annot_features_legend <- 
  propor_total_Mb_ERs_by_annot_features$all %>% get_legend()

# Save data -------------------------------------------------------------------------------------------

ggsave(plot = propor_total_Mb_ERs_by_annot_features_legend, filename = "propor_total_Mb_ERs_by_annot_features_legend.png", 
       path = "OMIM_paper/figures/annotation_features/",
       width = 8, height = 2, dpi = 600, units = "in")

ggsave(plot = propor_total_Mb_ERs_by_annot_features_arranged, filename = "propor_total_Mb_ERs_by_annot_features_arranged.png", 
       path = "OMIM_paper/figures/annotation_features/",
       width = 10, height = (11.69/1.5), dpi = 600, units = "in")

