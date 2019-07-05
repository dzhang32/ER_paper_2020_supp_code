library(tidyverse)
library(stringr)
library(annotatER)
library(GenomicRanges)
library(SummarizedExperiment)
library(derfinder)
library(rtracklayer)
library(ggpubr)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

exon_details_all_tissues <- 
  read_delim("results/optimise_derfinder_cut_off/exon_delta_details_all_tissues_cut_off_1_10_0.2_maxgap_0_100_10.csv", delim = ",")

gtex_tissue_name_formatting <- read_delim("raw_data/gtex_tissue_name_formatting/OMIM_gtex_tissue_name_formatting.csv", delim = ",")

ER_min_exon_delta_max_exon_delta_eq_0 <- 
  read_delim("results/optimise_derfinder_cut_off/exon_delta_details_optimised_maxgap_cutoff.csv", delim = ",")


##### First level #####

check_exon_biotypes <- function(){
  
  ensembl_grch38_v92_gtf_gr <- 
    import("/data/references/ensembl/gtf_gff3/v92/Homo_sapiens.GRCh38.92.gtf") 
  
  ensembl_grch38_v92_gtf_gr_exons <- ensembl_grch38_v92_gtf_gr[ensembl_grch38_v92_gtf_gr$type == "exon"]
  
  ensembl_grch38_v92_TxDb <- 
    generate_txDb_from_gtf(gtf_gff3_path = "/data/references/ensembl/gtf_gff3/v92/Homo_sapiens.GRCh38.92.gtf",
                           output_path = "/data/references/ensembl/txdb_sqlite/v92//ensembl_grch38_v92_txdb.sqlite",
                           seq_levels_to_keep = c(1:22, "X", "Y", "MT"), 
                           genome_build = "hg38")
  
  ensembl_grch38_v92_all_exons_gr <- 
    exons(ensembl_grch38_v92_TxDb, columns=c("EXONNAME", "GENEID"))
  
  ensembl_grch38_v92_exons_gr_non_overlapping <- 
    get_non_overlapping_exons(ensembl_grch38_v92_TxDb)
  
  ensembl_grch38_v92_all_genes_w_biotype <- 
    query_biomart(mart = 38, attributes = c("ensembl_gene_id", "gene_biotype"), filters = "ensembl_gene_id", 
                  values = ensembl_grch38_v92_all_exons_gr$GENEID %>% unlist() %>% unique())
  
  ensembl_grch38_v92_all_exons_w_biotype <- 
    get_exon_w_gene_biotype_data(ensembl_grch38_v92_all_exons_gr, ensembl_grch38_v92_all_genes_w_biotype) 
  
  ensembl_grch38_v92_all_exons_w_biotype_propor <- 
    ensembl_grch38_v92_all_exons_w_biotype %>% 
    mutate(exons = "all") %>% 
    group_by(exons, gene_biotype) %>% 
    summarise(propor_gene_biotype = n()/nrow(ensembl_grch38_v92_all_exons_w_biotype))
  
  ensembl_grch38_v92_exons_gr_non_overlapping_w_biotype <- 
    get_exon_w_gene_biotype_data(ensembl_grch38_v92_exons_gr_non_overlapping, ensembl_grch38_v92_all_genes_w_biotype) 
  
  ensembl_grch38_v92_exons_gr_non_overlapping_w_biotype_propor <- 
    ensembl_grch38_v92_exons_gr_non_overlapping_w_biotype %>% 
    mutate(exons = "nonoverlapping") %>% 
    group_by(exons, gene_biotype) %>% 
    summarise(propor_gene_biotype = n()/nrow(ensembl_grch38_v92_exons_gr_non_overlapping_w_biotype))
  
  ensembl_grch38_v92_all_exons_w_biotype_propor_w_non_overlapping <- 
    ensembl_grch38_v92_all_exons_w_biotype_propor %>% 
    bind_rows(ensembl_grch38_v92_exons_gr_non_overlapping_w_biotype_propor)
  
  gene_biotype_ab_0.01 <- 
    ensembl_grch38_v92_all_exons_w_biotype_propor_w_non_overlapping %>% 
    filter(propor_gene_biotype >= 0.01) %>% 
    .[["gene_biotype"]] %>% 
    unique()
  
  ensembl_grch38_v92_all_exons_w_biotype_propor_w_non_overlapping_summarised <- 
    ensembl_grch38_v92_all_exons_w_biotype_propor_w_non_overlapping %>% 
    mutate(gene_biotype_collapsed = ifelse(gene_biotype %in% gene_biotype_ab_0.01, gene_biotype, "all other biotypes")) %>% 
    group_by(exons, gene_biotype_collapsed) %>% 
    summarise(propor_gene_biotype = sum(propor_gene_biotype))
  
  gene_biotype_collapsed_desc_propor_gene_biotype <- 
    ensembl_grch38_v92_all_exons_w_biotype_propor_w_non_overlapping_summarised %>% 
    filter(exons == "all") %>% arrange(desc(propor_gene_biotype))
  
  ensembl_grch38_v92_exons_biotype_plot <- 
    ggplot(ensembl_grch38_v92_all_exons_w_biotype_propor_w_non_overlapping_summarised %>% 
             mutate(gene_biotype_collapsed = gene_biotype_collapsed %>% 
                      factor() %>% fct_relevel(gene_biotype_collapsed_desc_propor_gene_biotype$gene_biotype_collapsed)), 
           aes(x = gene_biotype_collapsed, y = propor_gene_biotype, fill = exons)) +
    geom_col(width = 1, colour = "black", position = "dodge") +
    scale_fill_manual(name = "Biotype", labels = c("All exons", "Non-overlapping exons"), 
                      values = get_palette("jco", 2)) +
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Proportion of exons (ensembl v92)") +
    theme_pubr(x.text.angle = 90, legend = "right", base_size = 10) +
    theme(axis.text.x = element_text(vjust = 0.5))
  
  return(ensembl_grch38_v92_exons_biotype_plot)
  
}

get_optimised_exon_graphs <- function(exon_details_all_tissues, ER_min_exon_delta_max_exon_delta_eq_0, tissue_to_plot, gtex_tissue_name_formatting){
  
  exon_details_all_tissues_tissue_filtered <- 
    exon_details_all_tissues %>% 
    filter(tissue == tissue_to_plot) %>% 
    left_join(gtex_tissue_name_formatting, by = c("tissue" = "OMIM_gtex_name"))
  
  ER_min_exon_delta_max_exon_delta_eq_0_filtered <- 
    ER_min_exon_delta_max_exon_delta_eq_0 %>% 
    filter(tissue == tissue_to_plot)
  
  exon_details_all_tissues_tissue_filtered_maxgap <- 
    exon_details_all_tissues_tissue_filtered %>% 
    filter(maxgap == ER_min_exon_delta_max_exon_delta_eq_0_filtered$maxgap)
  
  maxgaps_colours <- 
    data_frame(maxgap = seq(0, 100, 10),
               colours = get_palette("Blues", 14)[4:14]) %>% 
    mutate(colours = ifelse(maxgap == ER_min_exon_delta_max_exon_delta_eq_0_filtered$maxgap, "red", colours))
  
  exon_delta_median_plot <- 
    ggplot(exon_details_all_tissues_tissue_filtered, 
           aes(x = cut_off, y = exon_delta_median)) +
    geom_line(aes(colour = as.character(maxgap) %>% factor() %>% fct_relevel(as.character(seq(0, 100, 10))))) +
    geom_line(data = exon_details_all_tissues_tissue_filtered_maxgap, aes(x = cut_off, y = exon_delta_median), colour = "red") +
    geom_vline(xintercept = ER_min_exon_delta_max_exon_delta_eq_0_filtered$cut_off, colour = "#177D87", linetype = 2) +
    scale_x_continuous(name = "MCC") +
    scale_y_continuous(name = expression("Median"~Delta)) +
    scale_colour_manual("MRG", values = maxgaps_colours$colours) +
    theme_pubr(legend = "right") +
    theme(legend.title = element_text(colour = "red"), 
          axis.title.x = element_text(colour = "#177D87"))
  
  num_exon_delta_eq_0_plot <- 
    ggplot(exon_details_all_tissues_tissue_filtered, 
           aes(x = cut_off, y = num_exon_delta_eq_0)) +
    geom_line(aes(colour = as.character(maxgap) %>% factor() %>% fct_relevel(as.character(seq(0, 100, 10))))) +
    geom_line(data = exon_details_all_tissues_tissue_filtered_maxgap, aes(x = cut_off, y = num_exon_delta_eq_0), colour = "red") +
    geom_vline(xintercept = ER_min_exon_delta_max_exon_delta_eq_0_filtered$cut_off, colour = "#177D87", linetype = 2) +
    scale_x_continuous(name = "MCC") +
    scale_y_continuous(name = expression("Number of ERs with "~Delta~"= 0")) +
    scale_colour_manual("MRG", values = maxgaps_colours$colours) +
    theme_pubr(legend = "right") +
    theme(legend.title = element_text(colour = "red"), 
          axis.title.x = element_text(colour = "#177D87"))
  
  optimised_exon_delta_plots <- ggarrange(exon_delta_median_plot, num_exon_delta_eq_0_plot, nrow = 2, ncol = 1, align = "v", 
                                          common.legend = T, legend = "right")
  
  return(optimised_exon_delta_plots)
  
}

generate_cut_off_max_gap_example_plot <- function(){
  
  eg_position_vs_coverage_to_plot <- 
  data_frame(`Mean coverage` = 
               c(sample(x = 0:5, 250, replace = T), 
                 sample(x = 0:10, 50, replace = T), 
                 sample(x = 10:45, 50, replace = T), 
                 sample(x = 80:120, 250, replace = T), 
                 sample(x = 90:120, 100, replace = T),
                 sample(x = 30:45, 50, replace = T),
                 sample(x = 10:40, 50, replace = T),
                 sample(x = 10:15, 250, replace = T), 
                 sample(x = 10:40, 50, replace = T),
                 sample(x = 30:45, 50, replace = T),
                 sample(x = 120:140, 400, replace = T), 
                 sample(x = 0:10, 50, replace = T), 
                 sample(x = 0:5, 250, replace = T)), 
             `Genomic position` = 1:length(`Mean coverage`)) %>% 
    filter(`Genomic position` %in% seq(0, length(`Genomic position`), 10))

  
  eg_position_vs_coverage_plot_line <- 
    ggplot(eg_position_vs_coverage_to_plot, aes(x = `Genomic position`, y = `Mean coverage`)) +
    geom_line(colour = "black") + 
    geom_hline(yintercept = 55, colour = "#177D87", linetype = 2) +
    theme_pubr() +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  eg_position_vs_coverage_plot_polygon <- 
    ggplot(eg_position_vs_coverage_to_plot, aes(x = `Genomic position`, y = `Mean coverage`)) +
    geom_polygon(fill = "#1D71A8", colour = "black") + 
    geom_hline(yintercept = 55, colour = "#177D87", linetype = 2) +
    theme_pubr() +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank())
    
  eg_position_vs_coverage_plot <- 
    ggarrange(eg_position_vs_coverage_plot_polygon, eg_position_vs_coverage_plot_line + rremove("ylab"), nrow = 2, ncol = 1, align = "v")
  
  return(eg_position_vs_coverage_plot)
  
}

##### Second level #####

get_non_overlapping_exons <- function(ensembl_grch38_v92_TxDb){
  
  ensembl_grch38_v92_exons_gr <- 
    ensembl_grch38_v92_TxDb %>% exons(columns = c("EXONNAME", "GENEID"))
  
  ensembl_grch38_v92_exons_gr_marked_overlapping <- 
    mark_overlapping_genes_gr(gr_1 = ensembl_grch38_v92_exons_gr, gr_2 = ensembl_grch38_v92_exons_gr, identical_gr = T, 
                              maxgap = -1L, minoverlap = 1L)
  
  ensembl_grch38_v92_exons_gr_non_overlapping <- 
    ensembl_grch38_v92_exons_gr_marked_overlapping[ensembl_grch38_v92_exons_gr_marked_overlapping$overlap_gr2 == F]
  
  return(ensembl_grch38_v92_exons_gr_non_overlapping)
  
}

get_exon_w_gene_biotype_data <- function(ensembl_grch38_v92_all_exons_gr, ensembl_grch38_v92_all_genes_w_biotype){
  
  ensembl_grch38_v92_all_exons_w_biotype <- 
    elementMetadata(ensembl_grch38_v92_all_exons_gr) %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(GENEID = unlist(GENEID)) %>% 
    left_join(ensembl_grch38_v92_all_genes_w_biotype, by = c("GENEID" = "ensembl_gene_id"))
  
  return(ensembl_grch38_v92_all_exons_w_biotype)
  
}

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_txDb_from_gtf.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/query_biomart.R")

##### Third level #####

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/mark_overlapping_genes_gr.R")

# Main ------------------------------------------------------------------------------------------------

ensembl_grch38_v92_exons_biotype_plot <- check_exon_biotypes()

optimised_exon_delta_plots_crbl <- 
  get_optimised_exon_graphs(exon_details_all_tissues, ER_min_exon_delta_max_exon_delta_eq_0, 
                            tissue_to_plot = "brain_cerebellum", gtex_tissue_name_formatting)

eg_position_vs_coverage_plot <- 
  generate_cut_off_max_gap_example_plot()

# Load data -------------------------------------------------------------------------------------------

ggsave(plot = ensembl_grch38_v92_exons_biotype_plot, path = "OMIM_paper/supp_figures/exon_biotypes/",
       filename = "optimised_exon_delta_plots_crbl.png", width = (8.27), height = (11.69/2), dpi = 600, units = "in")

ggsave(plot = optimised_exon_delta_plots_crbl, path = "OMIM_paper/figures/optimisation_of_transcription_discovery/",
       filename = "optimised_exon_delta_plots_crbl.png", width = (8.27/2), height = (11.69/2), dpi = 600, units = "in")

ggsave(plot = eg_position_vs_coverage_plot, path = "OMIM_paper/figures/optimisation_of_transcription_discovery/",
       filename = "eg_position_vs_coverage_plot.png", width = (8.27/3.5), height = (11.69/3.5), dpi = 600, units = "in")

