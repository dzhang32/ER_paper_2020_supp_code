library(tidyverse)
library(stringr)
library(forcats)
library(ggpubr)
library(foreach)
library(derfinder)
library(GenomicRanges)
library(genomation)
library(rtracklayer)
library(gridExtra)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

load(file = "results/annotate_ERs/ERs_optimised_cut_off_max_gap_all_tissues_w_annot_df.rda")

tissue_optimal_cut_off_max_gap_df <- read_delim("results/optimise_derfinder_cut_off/exon_delta_details_optimised_maxgap_cutoff.csv", delim = ",")

scz_ERs_jaffe_gr <- readBed(file = "raw_data/scz_ERs_jaffe_for_validation/szControl_ERs.bed")

gtex_tissue_name_formatting <- read_delim("raw_data/gtex_tissue_name_formatting/OMIM_gtex_tissue_name_formatting.csv", delim = ",")

# Functions -------------------------------------------------------------------------------------------

##### First level ##### 

get_GTEx_hipp_ERs_quantification_Seb <- function(ERs_w_annotation_all_tissues, mean_cov_hipp_seb_paths){
  
  ERs_w_annotation_all_tissues_hipp <- 
    ERs_w_annotation_all_tissues %>% 
    filter(tissue == "hippocampus") 
  
  ERs_w_annotation_all_tissues_hipp_gr <- 
    GRanges(str_c(ERs_w_annotation_all_tissues_hipp$seqnames, ":", 
                  ERs_w_annotation_all_tissues_hipp$start, "-",
                  ERs_w_annotation_all_tissues_hipp$end)) 
    
  DZ_hipp_ERs_seb_quant_all_chrs <- GRanges()
  
  for(i in seq_along(mean_cov_hipp_seb_paths)){
    
    mean_cov_hipp_seb_path <- mean_cov_hipp_seb_paths[i]
    chr_to_filter <- mean_cov_hipp_seb_path %>% str_replace("/.*/", "") %>% str_replace("\\.rda", "") %>% str_c("chr", .)
    
    print(str_c(Sys.time(), " - getting coverage for Seb's hipp regions: ", chr_to_filter))
    
    load(mean_cov_hipp_seb_path)
    
    ERs_w_annotation_all_tissues_hipp_gr_chr_filtered <- 
      ERs_w_annotation_all_tissues_hipp_gr %>% keepSeqlevels(value = chr_to_filter, pruning.mode = "coarse")
    
    DZ_hipp_ERs_seb_quant_one_chr <- 
      get_region_coverage(regions = ERs_w_annotation_all_tissues_hipp_gr_chr_filtered, 
                          fstats = meanCov[[chr_to_filter]][["meanCoverage"]], 
                          position = Rle(TRUE, length(meanCov[[chr_to_filter]][["meanCoverage"]])), 
                          maxRegionGap = 0L, chr = chr_to_filter)
      
    DZ_hipp_ERs_seb_quant_all_chrs <- c(DZ_hipp_ERs_seb_quant_all_chrs, DZ_hipp_ERs_seb_quant_one_chr)
    
  }
  
  DZ_hipp_ERs_seb_quant_all_chrs <- 
    DZ_hipp_ERs_seb_quant_all_chrs %>% sortSeqlevels()

  return(DZ_hipp_ERs_seb_quant_all_chrs) 
  
}

check_GTEx_hipp_ER_validation <- 
  function(ERs_w_annotation_all_tissues, DZ_hipp_ERs_seb_quant_all_chrs, seb_hipp_optimised_cutoff = 3.3, feature_type_order_colors){
  
  ERs_w_annotation_all_tissues_hipp <- 
      ERs_w_annotation_all_tissues %>% 
      filter(tissue == "hippocampus")   
    
  ERs_w_annotation_all_tissues_hipp_chr_filtered <- 
  ERs_w_annotation_all_tissues_hipp %>% 
    filter(seqnames %in% str_c("chr", c(1:22, "X"))) %>% 
    mutate(seqnames = seqnames %>% factor() %>% fct_relevel(str_c("chr", c(1:22, "X")))) %>% 
    arrange(seqnames, start)
    
  ERs_w_annotation_all_tissues_hipp_chr_filtered_gr <- 
    GRanges(str_c(ERs_w_annotation_all_tissues_hipp_chr_filtered$seqnames, ":", 
                  ERs_w_annotation_all_tissues_hipp_chr_filtered$start, "-",
                  ERs_w_annotation_all_tissues_hipp_chr_filtered$end))
  
  stopifnot(identical(ranges(ERs_w_annotation_all_tissues_hipp_chr_filtered_gr), ranges(DZ_hipp_ERs_seb_quant_all_chrs)))
  
  ERs_w_annotation_all_tissues_hipp_chr_filtered_w_valid_seb_hipp <- 
    ERs_w_annotation_all_tissues_hipp_chr_filtered %>% 
    mutate(quant_seb_hipp = DZ_hipp_ERs_seb_quant_all_chrs$value, 
           validated_seb_hipp = quant_seb_hipp >= seb_hipp_optimised_cutoff, 
           annotationType_split_read_annot_fct = 
             ifelse(is.na(annotationType_split_read_annot), "none", 
                    ifelse(annotationType_split_read_annot == 1, "(partially) annotated split read", "unannotated split read")) %>% 
             factor() %>% 
             fct_relevel("(partially) annotated split read", "unannotated split read", "none")) 
  
  ERs_w_annotation_all_tissues_hipp_chr_filtered_w_valid_seb_hipp_to_plot <- 
    ERs_w_annotation_all_tissues_hipp_chr_filtered_w_valid_seb_hipp %>% 
    filter(width > 3) %>% 
    group_by(ensembl_grch38_v92_region_annot, annotationType_split_read_annot_fct) %>% 
    summarise(n_validated_seb_hipp = sum(validated_seb_hipp), 
              n_annotationType_split_read = n()) %>% 
    mutate(percent_validated_seb_hipp = (n_validated_seb_hipp/n_annotationType_split_read) * 100) 
  
  ERs_w_annotation_all_tissues_hipp_chr_filtered_w_valid_seb_hipp_plot <- 
    ggplot(ERs_w_annotation_all_tissues_hipp_chr_filtered_w_valid_seb_hipp_to_plot, 
         aes(x = annotationType_split_read_annot_fct, y = percent_validated_seb_hipp)) + 
    geom_col(aes(fill = annotationType_split_read_annot_fct), colour = "black") +
    geom_text(aes(x = annotationType_split_read_annot_fct, y = percent_validated_seb_hipp, 
                  label = round(percent_validated_seb_hipp, digits = 2) %>% str_c(., "%")), 
              vjust = -1, size = 3) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Percentage validated - quanitfied", limits = c(0, 100)) + 
    facet_wrap(~ensembl_grch38_v92_region_annot, ncol = 6) +
    scale_fill_manual(name = "Split read support",
                      values = get_palette("jco", 3)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ERs_w_annotation_all_tissues_hipp_chr_filtered_w_valid_seb_hipp_to_plot_intron_inter <- 
    ERs_w_annotation_all_tissues_hipp_chr_filtered_w_valid_seb_hipp_to_plot %>% 
    filter(ensembl_grch38_v92_region_annot %in% c("intergenic", "intron")) %>% 
    group_by(ensembl_grch38_v92_region_annot) %>% 
    summarise(n_annotationType_split_read = sum(n_annotationType_split_read), 
              n_validated_seb_hipp = sum(n_validated_seb_hipp)) %>% 
    mutate(percent_validated_seb_hipp = n_validated_seb_hipp/n_annotationType_split_read * 100)
  
  return(ERs_w_annotation_all_tissues_hipp_chr_filtered_w_valid_seb_hipp_plot)
  
}

get_percent_validation_ERs <- 
  function(ERs_w_annotation_all_tissues, seed_tissue, tissue_to_quantify = NULL, optimised_cutoff_manual = NULL, 
           mean_cov_tissue_dir, gtex_SRP051844_SRP058181, tissue_optimal_cut_off_max_gap_df, feature_type_order_colors){
  
  ERs_quant_all_chrs <- 
    quantify_ERs_all_chrs(ERs_w_annotation_all_tissues = ERs_w_annotation_all_tissues, seed_tissue = seed_tissue, 
                          mean_cov_tissue_dir = mean_cov_tissue_dir, gtex_SRP051844_SRP058181 = gtex_SRP051844_SRP058181)
  
  ERs_validated_plot <- 
    plot_ERs_percentage_validation(ERs_w_annotation_all_tissues = ERs_w_annotation_all_tissues, seed_tissue = seed_tissue, 
                                   gtex_SRP051844_SRP058181 = gtex_SRP051844_SRP058181, 
                                   tissue_to_quantify = tissue_to_quantify, ERs_quant_all_chrs = ERs_quant_all_chrs, 
                                   tissue_optimal_cut_off_max_gap_df = tissue_optimal_cut_off_max_gap_df, optimised_cutoff_manual = optimised_cutoff_manual, 
                                   feature_type_order_colors = feature_type_order_colors)
  
  return(ERs_validated_plot)
  

}

check_SNCA_region_validation <- function(ERs_w_annotation_all_tissues){
  
  ERs_quant_all_chrs <- 
    quantify_ERs_all_chrs(ERs_w_annotation_all_tissues %>% 
                            filter(ensembl_grch38_v92_region_annot != "exon, intergenic, intron"), 
                          seed_tissue = "frontalcortexba9", 
                          mean_cov_tissue_dir = "/data/recount/SRP051844/mean_cov/controls_by_chr_normalised_40mill/", 
                          gtex_SRP051844_SRP058181 = "SRP051844")
  
  ERs_SNCA_original_details <- 
    ERs_w_annotation_all_tissues %>% 
    filter(ensembl_grch38_v92_region_annot != "exon, intergenic, intron", 
           tissue == "frontalcortexba9",
           seqnames == "chr4", start >= (89724099 - 50000), end <= (89838315 + 10000))
  
  ERs_SNCA_independent_dataset_details <- 
    ERs_quant_all_chrs %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    filter(seqnames == "chr4", start >= (89724099 - 50000), end <= (89838315 + 10000))
  
  stopifnot(identical(ERs_SNCA_original_details$width, ERs_SNCA_independent_dataset_details$width))
  
  ERs_SNCA_original_w_independent_dataset <- 
    ERs_SNCA_original_details %>% 
    mutate(value_independent_dataset = ERs_SNCA_independent_dataset_details$value)
  
  ERs_SNCA_original_w_independent_dataset %>% dplyr::select(ensembl_grch38_v92_region_annot, value, value_independent_dataset, 
                                                            mean_CDTS_percentile, mean_phast_cons_7, mean_phast_cons_100)
  
}

check_validation_jaffe_scz_ERs <- function(ERs_w_annotation_all_tissues, scz_ERs_jaffe_gr){
  
  ERs_w_annotation_all_tissues_tissue_filtered <- 
    ERs_w_annotation_all_tissues %>% 
    filter(tissue == "frontalcortexba9", width >= 12,
           seqnames %in% str_c("chr", c(1:22, "X", "Y")))
  
  ERs_w_annotation_all_tissues_tissue_filtered_gr <- 
    GRanges(str_c(ERs_w_annotation_all_tissues_tissue_filtered$seqnames, ":", 
                  ERs_w_annotation_all_tissues_tissue_filtered$start, "-",
                  ERs_w_annotation_all_tissues_tissue_filtered$end))
  
  hg19_to_Hg38_chain <- import.chain("/data/liftover/hg19/hg19ToHg38.over.chain")
  scz_ERs_jaffe_gr_hg38 <- 
    liftOver(scz_ERs_jaffe_gr, hg19_to_Hg38_chain) %>% 
    unlist()
  
  ERs_w_annotation_all_tissues_tissue_filtered_gr_validated <- 
    mark_overlapping_genes_gr(gr_1 = ERs_w_annotation_all_tissues_tissue_filtered_gr, gr_2 = scz_ERs_jaffe_gr_hg38, identical_gr = F, minoverlap = 1, maxgap = -1)
  
  
  ERs_w_annotation_all_tissues_tissue_filtered_w_valid <- 
    ERs_w_annotation_all_tissues_tissue_filtered %>% 
    mutate(validated_seb_hipp = ERs_w_annotation_all_tissues_tissue_filtered_gr_validated$overlap_gr2,
           annotationType_split_read_annot_fct = 
             ifelse(is.na(annotationType_split_read_annot), "none", "supported") %>% 
             factor() %>% 
             fct_relevel(c("supported", "none"))) 
  
  ERs_w_annotation_all_tissues_tissue_filtered_w_valid_to_plot <- 
    ERs_w_annotation_all_tissues_tissue_filtered_w_valid %>% 
    filter(width > 3) %>% 
    group_by(ensembl_grch38_v92_region_annot, annotationType_split_read_annot_fct) %>% 
    summarise(n_validated_seb_hipp = sum(validated_seb_hipp), 
              n_annotationType_split_read = n()) %>% 
    mutate(percent_validated_seb_hipp = (n_validated_seb_hipp/n_annotationType_split_read) * 100) 
  
  ERs_w_annotation_all_tissues_tissue_filtered_w_valid_plot <- 
    ggplot(ERs_w_annotation_all_tissues_tissue_filtered_w_valid_to_plot, 
           aes(x = annotationType_split_read_annot_fct, y = percent_validated_seb_hipp)) + 
    geom_col(aes(fill = annotationType_split_read_annot_fct), colour = "black") +
    geom_text(aes(x = annotationType_split_read_annot_fct, y = percent_validated_seb_hipp, 
                  label = round(percent_validated_seb_hipp, digits = 2) %>% str_c(., "%")), 
              vjust = -1, size = 3) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Percentage validated - quanitfied", limits = c(0, 100)) + 
    facet_wrap(~ensembl_grch38_v92_region_annot, ncol = 6) +
    scale_fill_manual(name = "Split read support",
                      values = get_palette("jco", 3)[c(1, 3)]) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(ERs_w_annotation_all_tissues_tissue_filtered_w_valid_plot)
  
}

plot_num_unannot_vs_annot_v87_v92 <- function(ERs_w_annotation_all_tissues, feature_type_order_colors, gtex_tissue_name_formatting){ 
  
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
    summarise(total_Kb = sum(width)/1000) %>% 
    ungroup() %>%
    mutate(total_Kb = ifelse(ensembl_v87_v92_annot %in% c("exon in v87", "exon, intergenic in v87", "exon, intron in v87"), -total_Kb, total_Kb)) %>% 
    left_join(gtex_tissue_name_formatting_w_face %>% dplyr::select(OMIM_gtex_name, gtex_tissues_name_to_plot, gtex_tissue_group, tissue_color_hex, face), by = c("tissue" = "OMIM_gtex_name")) %>%
    mutate(gtex_tissues_name_to_plot = gtex_tissues_name_to_plot %>% str_replace(" $", "") %>% str_c("   "))
  
  ERs_v87_vs_v92_to_plot_order <- 
    ERs_v87_vs_v92_to_plot %>% 
    filter(!ensembl_v87_v92_annot %in% c("exon in v87", "exon, intergenic in v87", "exon, intron in v87")) %>% 
    group_by(gtex_tissues_name_to_plot, face, tissue_color_hex, tissue) %>% 
    summarise(intron_inter_total_Kb = sum(total_Kb)) %>% 
    arrange(desc(intron_inter_total_Kb)) %>% 
    ungroup()
  
  ERs_v87_vs_v92_total_Kb <- 
    ggplot(ERs_v87_vs_v92_to_plot %>% 
           mutate(gtex_tissues_name_to_plot = gtex_tissues_name_to_plot %>% 
                    factor() %>% 
                    fct_relevel(ERs_v87_vs_v92_to_plot_order$gtex_tissues_name_to_plot)), 
         aes(x = gtex_tissues_name_to_plot, y = total_Kb)) +
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
  
  return(ERs_v87_vs_v92_total_Kb)
  
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

##### Second level ##### 

quantify_ERs_all_chrs <- function(ERs_w_annotation_all_tissues, seed_tissue, mean_cov_tissue_dir, gtex_SRP051844_SRP058181){
  
  ERs_w_annotation_all_tissues_tissue_filtered <- 
    ERs_w_annotation_all_tissues %>% 
    filter(tissue == seed_tissue, 
           seqnames %in% str_c("chr", c(1:22, "X", "Y")))
  
  ERs_w_annotation_all_tissues_tissue_filtered_gr <- 
    GRanges(str_c(ERs_w_annotation_all_tissues_tissue_filtered$seqnames, ":", 
                  ERs_w_annotation_all_tissues_tissue_filtered$start, "-",
                  ERs_w_annotation_all_tissues_tissue_filtered$end))
  
  mean_cov_tissue_by_chr_paths <- 
    list.files(mean_cov_tissue_dir, full.names = T)
  
  ERs_quant_all_chrs <- GRanges()
  
  for(i in seq_along(mean_cov_tissue_by_chr_paths)){
    
    mean_cov_hipp_seb_path <- mean_cov_tissue_by_chr_paths[i]
    
    if(gtex_SRP051844_SRP058181 == "gtex"){
      
      tissue_to_quantify <- 
        mean_cov_hipp_seb_path %>% 
        str_replace("/.*/", "") %>% 
        str_replace("_chr.*", "") %>% 
        str_replace("gtex_", "")
      
      chr_to_filter <- 
        mean_cov_hipp_seb_path %>% 
        str_replace("/.*/", "") %>% 
        str_replace("_mean_cov\\.rda", "") %>% 
        str_replace(str_c("gtex_", tissue_to_quantify, "_"), "")
      
    }else if(gtex_SRP051844_SRP058181 == "SRP051844"){
      
      chr_to_filter <- 
        mean_cov_hipp_seb_path %>% 
        str_replace("/.*/", "") %>% 
        str_replace("SRP051844_control_mean_cov_40mill_", "") %>% 
        str_replace("\\.rda", "")
      
    }else if(gtex_SRP051844_SRP058181 == "SRP058181"){
      
      chr_to_filter <- 
        mean_cov_hipp_seb_path %>% 
        str_replace("/.*/", "") %>% 
        str_replace("SRP058181_control_mean_cov_40mill_", "") %>% 
        str_replace("\\.rda", "")
      
    }
    
    if(chr_to_filter == "chrM"){
      
      next
      
    }
    
    print(str_c(Sys.time(), " - getting coverage for ", seed_tissue, " regions in ", chr_to_filter))
    
    load(mean_cov_hipp_seb_path)
    
    ERs_w_annotation_all_tissues_tissue_filtered_gr_chr_filtered <- 
      ERs_w_annotation_all_tissues_tissue_filtered_gr %>% keepSeqlevels(value = chr_to_filter, pruning.mode = "coarse")
    
    if(gtex_SRP051844_SRP058181 == "gtex"){
      
      ERs_quant_one_chr <- 
        get_region_coverage(regions = ERs_w_annotation_all_tissues_tissue_filtered_gr_chr_filtered, 
                            fstats = tissue_coverage_w_mean_normalised[["meanCoverage"]], 
                            position = Rle(TRUE, length(tissue_coverage_w_mean_normalised[["meanCoverage"]])), 
                            maxRegionGap = 0L, chr = chr_to_filter)
      
    }else if(gtex_SRP051844_SRP058181 == "SRP051844"){
      
      ERs_quant_one_chr <- 
        get_region_coverage(regions = ERs_w_annotation_all_tissues_tissue_filtered_gr_chr_filtered, 
                            fstats = SRP051844_control_mean_cov_normalised[["meanCoverage"]], 
                            position = Rle(TRUE, length(SRP051844_control_mean_cov_normalised[["meanCoverage"]])), 
                            maxRegionGap = 0L, chr = chr_to_filter)
      
    }else if(gtex_SRP051844_SRP058181 == "SRP058181"){
      
      ERs_quant_one_chr <- 
        get_region_coverage(regions = ERs_w_annotation_all_tissues_tissue_filtered_gr_chr_filtered, 
                            fstats = SRP058181_control_mean_cov_normalised[["meanCoverage"]], 
                            position = Rle(TRUE, length(SRP058181_control_mean_cov_normalised[["meanCoverage"]])), 
                            maxRegionGap = 0L, chr = chr_to_filter)
      
    }
    
    ERs_quant_all_chrs <- c(ERs_quant_all_chrs, ERs_quant_one_chr)
    
  }
  
  return(ERs_quant_all_chrs)
  
}

plot_ERs_percentage_validation <- 
  function(ERs_w_annotation_all_tissues, seed_tissue, gtex_SRP051844_SRP058181, tissue_to_quantify, optimised_cutoff_manual,
           ERs_quant_all_chrs, tissue_optimal_cut_off_max_gap_df, feature_type_order_colors){
    
    # check all ER definitions are still identical to previous (i.e. quantification has not changed these)
    ERs_w_annotation_all_tissues_tissue_filtered <- 
      ERs_w_annotation_all_tissues %>% 
      filter(tissue == seed_tissue, 
             seqnames %in% str_c("chr", c(1:22, "X", "Y")))
    
    ERs_w_annotation_all_tissues_tissue_filtered_gr <- 
      GRanges(str_c(ERs_w_annotation_all_tissues_tissue_filtered$seqnames, ":", 
                    ERs_w_annotation_all_tissues_tissue_filtered$start, "-",
                    ERs_w_annotation_all_tissues_tissue_filtered$end))
    
    ERs_quant_all_chrs_sorted <- sort(sortSeqlevels(ERs_quant_all_chrs))

    stopifnot(identical(ranges(ERs_w_annotation_all_tissues_tissue_filtered_gr), ranges(ERs_quant_all_chrs_sorted)))
    
    if(gtex_SRP051844_SRP058181 == "gtex"){
      
      optimised_cutoff <- 
        tissue_optimal_cut_off_max_gap_df %>% 
        filter(tissue == tissue_to_quantify) %>% 
        .[["cut_off"]]
      
    }else{
      
      optimised_cutoff <- optimised_cutoff_manual
      
    }
    
    
    ERs_w_annotation_all_tissues_tissue_filtered_w_valid <- 
      ERs_w_annotation_all_tissues_tissue_filtered %>% 
      mutate(quant_seb_hipp = ERs_quant_all_chrs_sorted$value, 
             validated_seb_hipp = quant_seb_hipp >= optimised_cutoff, 
             annotationType_split_read_annot_fct = 
               ifelse(annotationType_split_read_annot == "none", "none", "split read support") %>% 
               factor() %>% 
               fct_relevel(c("split read support", "none"))) 
    
    ERs_w_annotation_all_tissues_tissue_filtered_w_valid_to_plot <- 
      ERs_w_annotation_all_tissues_tissue_filtered_w_valid %>% 
      group_by(ensembl_grch38_v92_region_annot, annotationType_split_read_annot_fct) %>% 
      summarise(n_validated_seb_hipp = sum(validated_seb_hipp), 
                n_annotationType_split_read = n()) %>% 
      ungroup() %>% 
      mutate(percent_validated_seb_hipp = (n_validated_seb_hipp/n_annotationType_split_read) * 100, 
             ensembl_grch38_v92_region_annot = 
               ensembl_grch38_v92_region_annot %>% 
               factor() %>% 
               fct_relevel(feature_type_order_colors$feature_type_order)) 
    
    ERs_w_annotation_all_tissues_tissue_filtered_w_valid_plot <- 
      ggplot(ERs_w_annotation_all_tissues_tissue_filtered_w_valid_to_plot, 
             aes(x = ensembl_grch38_v92_region_annot, y = percent_validated_seb_hipp)) + 
      geom_col(aes(fill = ensembl_grch38_v92_region_annot, alpha = annotationType_split_read_annot_fct), position = "dodge", 
               colour = "black", width = 0.75) +
      geom_text(aes(x = ensembl_grch38_v92_region_annot, y = percent_validated_seb_hipp, group = annotationType_split_read_annot_fct, 
                    label = round(percent_validated_seb_hipp, digits = 0)), 
                vjust = -0.5, size = 3, position = position_dodge(width = 0.75)) +
      scale_x_discrete(name = "") +
      scale_y_continuous(name = "Validation of ERs (%)", limits = c(0, 100)) + 
      scale_fill_manual(name = "Annotation feature",
                        values = feature_type_order_colors$feature_type_colour) +
      scale_alpha_manual(name = "Split read support",
                         values = c(1, 0.4)) +
      theme_pubr() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    
    return(ERs_w_annotation_all_tissues_tissue_filtered_w_valid_plot)
    
  }

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/mark_overlapping_genes_gr.R")

##### Third level ##### 

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/get_region_coverage.R")

# Main ------------------------------------------------------------------------------------------------

feature_type_order_colors <- 
  data_frame(feature_type_order = c("exon", "exon, intron", "exon, intergenic", 
                                    "exon, intergenic, intron", 
                                    "intergenic", "intron"), 
             feature_type_colour = get_palette("jco", 9)[c(6, 1, 5, 3, 9, 4)])


ERs_w_annotation_all_tissues <- 
  ERs_w_annotation_all_tissues %>%
  filter(width > 3, !str_detect(tissue, "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast"),
         ensembl_grch38_v92_region_annot != "exon, intergenic, intron")

##### DZ HIPP ERs vs SG HIPP ERs #####

mean_cov_hipp_seb_paths <- list.files("/home/sguelfi/projects/R/hipp/data/derfinder/fullCoverage/meanCoverage/", full.names = T, recursive = T)

DZ_hipp_ERs_seb_quant_all_chrs <- get_GTEx_hipp_ERs_quantification_Seb(ERs_w_annotation_all_tissues, mean_cov_hipp_seb_paths)

# load("results/analyse_ER_annotation/validate_ERs/GTEx_hipp_ERs_quant_w_seb_hipp_all_chrs.rda")

ERs_w_annotation_all_tissues_hipp_chr_filtered_w_valid_seb_hipp_plot <- 
  check_GTEx_hipp_ER_validation(ERs_w_annotation_all_tissues_hipp, DZ_hipp_ERs_seb_quant_all_chrs, seb_hipp_optimised_cutoff = 1.2, feature_type_order_colors)

save(DZ_hipp_ERs_seb_quant_all_chrs, file = "results/analyse_ER_annotation/validate_ERs/GTEx_hipp_ERs_quant_w_seb_hipp_all_chrs.rda")

##### CRBL CTX technical duplicates #####

# cerebellum quantified using cerebellar hemisphere
# cerebellum_ERs_validated_cerebellar_hemi  <- 
#   get_percent_validation_ERs(ERs_w_annotation_all_tissues, seed_tissue = "brain_cerebellum", tissue_to_quantify = "brain_cerebellar_hemisphere", 
#                              mean_cov_tissue_dir = "/data/recount/GTEx_SRP012682/gtex_mean_coverage/by_tissue_smfrze_use_me/brain_cerebellar_hemisphere/",
#                              tissue_optimal_cut_off_max_gap_df = tissue_optimal_cut_off_max_gap_df, gtex_SRP051844_SRP058181 = "gtex", 
#                              feature_type_order_colors = feature_type_order_colors)


# cerebellar hemisphere quantified using cerebellum
cerebellar_hemi_ERs_validated_cerebellum <- 
  get_percent_validation_ERs(ERs_w_annotation_all_tissues %>% 
                               filter(ensembl_grch38_v92_region_annot != "exon, intergenic, intron"), 
                             seed_tissue = "brain_cerebellar_hemisphere", tissue_to_quantify = "brain_cerebellum",
                             mean_cov_tissue_dir = "/data/recount/GTEx_SRP012682/gtex_mean_coverage/by_tissue_smfrze_use_me/brain_cerebellum/",
                             tissue_optimal_cut_off_max_gap_df = tissue_optimal_cut_off_max_gap_df, gtex_SRP051844_SRP058181 = "gtex", 
                             feature_type_order_colors = feature_type_order_colors %>% filter(feature_type_order != c("exon, intergenic, intron")))

# frontal cortex quantified using cortex
# frontal_cortex_ERs_quantified <- 
#   get_percent_validation_ERs(ERs_w_annotation_all_tissues %>% filter(ensembl_grch38_v92_region_annot != "exon, intergenic, intron"), 
#                              seed_tissue = "frontalcortexba9", tissue_to_quantify = "cortex",
#                              mean_cov_tissue_dir = "/data/recount/GTEx_SRP012682/gtex_mean_coverage/by_tissue_smfrze_use_me/brain_cortex/",
#                              tissue_optimal_cut_off_max_gap_df = tissue_optimal_cut_off_max_gap_df, gtex_SRP051844_SRP058181 = "gtex", 
#                              feature_type_order_colors = feature_type_order_colors %>% filter(feature_type_order != "exon, intergenic, intron"))

# cortex quantified using frontal cortex
cortex_ERs_quantified_frontal_crtx <-
  get_percent_validation_ERs(ERs_w_annotation_all_tissues %>% 
                               filter(ensembl_grch38_v92_region_annot != "exon, intergenic, intron"),
                             seed_tissue = "cortex", tissue_to_quantify = "frontalcortexba9",
                             mean_cov_tissue_dir = "/data/recount/GTEx_SRP012682/gtex_mean_coverage/by_tissue_smfrze_use_me/brain_frontal_cortex_ba9/",
                             tissue_optimal_cut_off_max_gap_df = tissue_optimal_cut_off_max_gap_df, gtex_SRP051844_SRP058181 = "gtex",
                             feature_type_order_colors = feature_type_order_colors %>% filter(feature_type_order != c("exon, intergenic, intron")))



##### Using SRP051844 #####

frontal_cortex_ERs_quantified_SRP051844 <- 
  get_percent_validation_ERs(ERs_w_annotation_all_tissues %>% 
                               filter(ensembl_grch38_v92_region_annot != "exon, intergenic, intron"), 
                             seed_tissue = "frontalcortexba9", 
                             mean_cov_tissue_dir = "/data/recount/SRP051844/mean_cov/controls_by_chr_normalised_40mill/", 
                             gtex_SRP051844_SRP058181 = "SRP051844", optimised_cutoff_manual = 1.4, 
                             feature_type_order_colors = feature_type_order_colors %>% filter(feature_type_order != c("exon, intergenic, intron")))

# check if SNCA is one of the validated fctx ERs

##### Using SRP058181 #####

# frontal_cortex_ERs_quantified_SRP058181 <- 
#   get_percent_validation_ERs(ERs_w_annotation_all_tissues %>% filter(ensembl_grch38_v92_region_annot != "exon, intergenic, intron"), 
#                              seed_tissue = "frontalcortexba9", 
#                              mean_cov_tissue_dir = "/data/recount/SRP058181/mean_cov/controls_by_chr_normalised_40mill/", 
#                              gtex_SRP051844_SRP058181 = "SRP058181", optimised_cutoff_manual = 1.4, 
#                              feature_type_order_colors %>% filter(feature_type_order != "exon, intergenic, intron"))

##### Andrew jaffe SCZ ERs ######

validation_jaffe_scz_ERs_plot <- 
  check_validation_jaffe_scz_ERs(ERs_w_annotation_all_tissues, scz_ERs_jaffe_gr)

# Save data -------------------------------------------------------------------------------------------

ggsave(plot = ER_validation_independent_data_set_tech_rep, filename = "ER_validation_independent_data_set_tech_rep.png", 
       path = "OMIM_paper/figures/validation_of_ERs/",
       width =  8.27, height = (11.69/2.6), dpi = 600, units = "in")

ggsave(plot = frontal_cortex_ERs_quantified_SRP051844, filename = "frontal_cortex_ERs_quantified_SRP051844.png", 
       path = "OMIM_paper/figures/validation_of_ERs/",
       width =  8.27/1.3, height = (11.69/2.2), dpi = 600, units = "in")

