library(tidyverse)
library(stringr)
library(GenomicRanges)
library(SummarizedExperiment)
library(GenomicFeatures)
library(recount)
library(derfinder)
library(annotatER)
library(data.table)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

load("results/optimise_derfinder_cut_off/ERs_optimal_cut_off_maxgap_all_tissues.rda")

load("results/download_tidy_OMIM_data/OMIM_gene_start_stop_gr_marked_overlapping_genes_ens_v92_2015_05_29.Rdata")

# Functions -------------------------------------------------------------------------------------------

source("/home/dzhang/projects/constraint_conservation_wd/constraint_conservation/constraint/constraint_general/constraint_general_functions.R")
source("/home/dzhang/projects/constraint_conservation_wd/constraint_conservation/conservation/conservation_general/conservation_general_functions.R")

get_gtex_split_read_table_mean_cov_df <- function(){
  
  gtex_split_read_table_annotated_paths <- 
    list.files("/data/recount/GTEx_SRP012682/gtex_split_read_table_annotated/", full.names = T)
  
  gtex_split_read_table_df <- 
    data_frame(gtex_split_read_table_annotated_paths = gtex_split_read_table_annotated_paths, 
               tissue = 
                 gtex_split_read_table_annotated_paths %>% 
                 str_replace("/.*/", "") %>% 
                 str_replace("_split_read_table_annotated.csv", "") %>% 
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
    inner_join(gtex_mean_cov_df)
  
  return(gtex_split_read_table_mean_cov_df)
  
}

generate_genomic_states_list <- function(txdb_sqlite_paths){
  
  genomic_states_list <- list()
  
  for(i in seq_along(txdb_sqlite_paths)){
    
    txdb_sqlite_path <- txdb_sqlite_paths[i]
    
    
    print(str_c(Sys.time(), " - ", i, " - ", txdb_sqlite_path))
    
    txdb_ref_annot <- loadDb(file = txdb_sqlite_path)
    
    genomic_state_ref_annot <- 
      generate_genomic_state(txdb_ref_annot, 
                             output_path = 
                               txdb_sqlite_path %>% 
                               str_replace("txdb_sqlite", "genomic_state") %>% 
                               str_replace("_txdb\\.sqlite", "_genomic_state.rda"))
    
    genomic_states_list <- 
      c(genomic_states_list, genomic_state_ref_annot)
    
  }
  
  names(genomic_states_list) <- 
    txdb_sqlite_paths %>% 
    str_replace("/.*/", "") %>% 
    str_replace("_txdb\\.sqlite", "")
  
  return(genomic_states_list)
  
}

get_coverage_for_GTEx_ERs <- function(ERs_gr, tissue_to_filter, gtex_split_read_table_mean_cov_df){
  
  # format tissue to get the path of the split read table from the df 
  tissue_to_filter <- 
    tissue_to_filter %>% 
    str_replace_all("_", "") %>% 
    str_replace_all("-", "") %>% 
    str_replace("brain", "")
  
  gtex_split_read_table_mean_cov_df_tissue_filtered <- 
    gtex_split_read_table_mean_cov_df %>% 
    filter(tissue == tissue_to_filter)
  
  if(nrow(gtex_split_read_table_mean_cov_df_tissue_filtered) != 1){
    
    stop(str_c("Either no or more than 1 split read table path found for: ", tissue_to_filter))
    
  }
  
  mean_cov_tissue_by_chr_paths <- 
    list.files(gtex_split_read_table_mean_cov_df_tissue_filtered$gtex_mean_cov_paths, full.names = T)
  
  ERs_quant_all_chrs <- GRanges()
  
  for(j in seq_along(mean_cov_tissue_by_chr_paths)){
    
    mean_cov_tissue_by_chr_path <- mean_cov_tissue_by_chr_paths[j]
    
    tissue_to_quantify <- 
      mean_cov_tissue_by_chr_path %>% 
      str_replace("/.*/", "") %>% 
      str_replace("_chr.*", "") %>% 
      str_replace("gtex_", "")
    
    chr_to_filter <- 
      mean_cov_tissue_by_chr_path %>% 
      str_replace("/.*/", "") %>% 
      str_replace("_mean_cov\\.rda", "") %>% 
      str_replace(str_c("gtex_", tissue_to_quantify, "_"), "")
    
    if(chr_to_filter == "chrM"){
      
      next
      
    }
    
    print(str_c(Sys.time(), " - getting coverage for ", tissue_to_quantify, " ERs in ", chr_to_filter))
    
    load(mean_cov_tissue_by_chr_path)
    
    ERs_gr_chr_filtered <- 
      ERs_gr %>% keepSeqlevels(value = chr_to_filter, pruning.mode = "coarse")
    
    ERs_quant_one_chr <- 
      get_region_coverage(regions = ERs_gr_chr_filtered, 
                          fstats = tissue_coverage_w_mean_normalised[["meanCoverage"]], 
                          position = Rle(TRUE, length(tissue_coverage_w_mean_normalised[["meanCoverage"]])), 
                          maxRegionGap = 0L, chr = chr_to_filter)
    
    ERs_quant_all_chrs <- c(ERs_quant_all_chrs, ERs_quant_one_chr)
    
  }
  
  ERs_quant_all_chrs_sorted <- ERs_quant_all_chrs %>% sortSeqlevels() %>% sort()
  
  stopifnot(identical(ranges(ERs_gr), ranges(ERs_quant_all_chrs_sorted)))
  
  elementMetadata(ERs_quant_all_chrs_sorted)[["area"]] <- NULL
  elementMetadata(ERs_quant_all_chrs_sorted)[["indexStart"]] <- NULL
  elementMetadata(ERs_quant_all_chrs_sorted)[["indexEnd"]] <- NULL
  
  return(ERs_quant_all_chrs_sorted)
  
  
}

annotate_ERs_w_genomic_state <- function(ERs_gr, genomic_states_list){
  
  for(j in 1:length(genomic_states_list)){
    
    reference_annot_ver <- genomic_states_list[j] %>% names()
    genomic_state <- genomic_states_list[[j]]
    
    print(str_c(Sys.time(), " - ", j, " - ", reference_annot_ver))
  
    ER_annotation_count_table <- 
      annotateRegions(regions = ERs_gr, 
                        genomicState = genomic_state$fullGenome, 
                        maxgap = -1L, minoverlap = 1L, 
                        annotate = F)
      
    ER_annotation_count_table_w_region_annot <- 
      convert_annot_count_table_to_region_annot(count_table = ER_annotation_count_table$countTable)
      
    stopifnot(nrow(ER_annotation_count_table_w_region_annot) == length(ERs_gr))
      
    elementMetadata(ERs_gr)[[str_c(reference_annot_ver, "_region_annot")]] <- 
      ER_annotation_count_table_w_region_annot[["region_annot"]]
    

    
  }
  
  return(ERs_gr)
  
}

annotate_ERs_with_overlap_and_dis_to_genes <- function(ERs_gr, genes_gr, column_name_prefix){
  
  # marking those that overlap OMIM genes 
  ERs_gr_w_overlapping_genes_marked <- 
    mark_overlapping_genes_gr(gr_1 = ERs_gr, gr_2 = genes_gr, 
                              identical_gr = F, maxgap = -1L, minoverlap = 1L, ignore.strand = T) 
  
  elementMetadata(ERs_gr_w_overlapping_genes_marked)[[str_c("overlap_", column_name_prefix, "lgl")]] <- 
    ERs_gr_w_overlapping_genes_marked$overlap_gr2
  
  ERs_gr_w_overlapping_genes_marked$overlap_gr2 <- NULL
  
  # get overlapping gene names, here will concat all the gene names that that ER overlaps seperated by a , 
  ERs_overlapping_genes_hits <- findOverlaps(query = ERs_gr, subject = genes_gr, maxgap = -1, minoverlap = 1)
  
  ERs_overlapping_genes_hits_w_gene_name <- 
    ERs_overlapping_genes_hits %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(name = genes_gr %>% names() %>% .[subjectHits(ERs_overlapping_genes_hits)]) %>% 
    group_by(queryHits) %>% 
    summarise(name = str_c(name, collapse = ", ")) %>% 
    left_join(data_frame(ER_index = 1:length(ERs_gr)), ., by = c("ER_index" = "queryHits"))
  
  elementMetadata(ERs_gr_w_overlapping_genes_marked)[[str_c("overlap_", column_name_prefix, "name")]] <- 
    ERs_overlapping_genes_hits_w_gene_name$name
  
  # getting distance to genes, should return 0 for overlapping but also 0 for those exactly adjacent
  # picks a (to my knowledge random) gene name to annotate 
  ERs_gr_dis_to_nearest_hits <- 
    distanceToNearest(ERs_gr, genes_gr)
  
  ERs_gr_dis_to_nearest_details_df <- 
    data_frame(ER_index = 1:length(ERs_gr)) %>% 
    left_join(ERs_gr_dis_to_nearest_hits %>% 
                as.data.frame() %>% 
                as_tibble() %>% 
                mutate(name = genes_gr %>% names() %>% .[subjectHits(ERs_gr_dis_to_nearest_hits)]), 
              by = c("ER_index" = "queryHits")) %>% 
    dplyr::select(-ER_index, -subjectHits)
  
  stopifnot(length(ERs_gr_w_overlapping_genes_marked) == nrow(ERs_gr_dis_to_nearest_details_df))
  
  # correcting colnames
  colnames(ERs_gr_dis_to_nearest_details_df) <- 
    str_c("nearest_", column_name_prefix, colnames(ERs_gr_dis_to_nearest_details_df))
  
  # add distance and nearest OMIM gene to ER_gr
  
  for(j in seq_along(colnames(ERs_gr_dis_to_nearest_details_df))){
    
    col_to_add <- colnames(ERs_gr_dis_to_nearest_details_df)[j]
    
    elementMetadata(ERs_gr_w_overlapping_genes_marked)[[col_to_add]] <- 
      ERs_gr_dis_to_nearest_details_df[[col_to_add]]
    
  }
  
  return(ERs_gr_w_overlapping_genes_marked)
  
}

annotate_ERs_w_split_reads <- function(ERs_gr, tissue_to_filter, gtex_split_read_table_mean_cov_df){
  
  # format tissue to get the path of the split read table from the df 
  tissue_to_filter <- 
    tissue_to_filter %>% 
    str_replace_all("_", "") %>% 
    str_replace_all("-", "") %>% 
    str_replace("brain", "")
  
  gtex_split_read_table_df_tissue_filtered <- 
    gtex_split_read_table_mean_cov_df %>% 
    filter(tissue %>% str_detect(str_c("^", tissue_to_filter, "$"))) 
  
  if(nrow(gtex_split_read_table_df_tissue_filtered) != 1){
    
    print(str_c("Either no or more than 1 split read table path found for: ", tissue_to_filter))
    
  }

  gtex_split_read_table_annotated <- 
    read_format_gtex_split_read_table(gtex_split_read_table_df_tissue_filtered$gtex_split_read_table_annotated_paths)
  
  # get only the junction annotation columns (removing the coverage), here this should probably be incoporated into the annotatER pacakage at some point
  gtex_split_read_table_annotated_only_junc_coverage <-
  gtex_split_read_table_annotated %>% 
  dplyr::select(junID, chr, start, stop, strand, countsSamples, acceptor, donor, junction, precBoundDonor, precBoundAcceptor) %>%
      filter(!is.na(junID))
  
  # annotate ERs with split reads 
  ERs_gr_split_read_annotated <- 
    annotatERJunction(regions = ERs_gr, 
                      annotatedJunction = data.frame(gtex_split_read_table_annotated_only_junc_coverage), return_df_or_gr = "df")
  
  ERs_gr_split_read_annotated_tidy <- 
    ERs_gr_split_read_annotated %>% 
    as_tibble() %>% 
    mutate(annotationType = ifelse(is.na(annotationType), "none",
                                  ifelse(annotationType == 1, "partially annotated split read", "uannotated split read"))) %>% 
    dplyr::select(-one_of("junID","countsSamples","acceptor","donor","junction",
                         "precBoundDonor","precBoundAcceptor","pred_strand",
                         "multi_junc","junER","distance","transcript"))
  
  # add split_read_annot to colnames to distinguish between other kinds of annot
  colnames_raw <- colnames(ERs_gr_split_read_annotated_tidy)
  
  colnames_split_read_annot <- 
    ifelse(colnames_raw %in% 
             c("annotationType", "p_annot_junc_ids", "p_annot_junc_count_samples", "uniq_transcripts", 
               "p_annot_all_transcripts", "unannot_junc_ids", "unannot_junc_count_samples", "unannot_junc_connect_ER_index"), 
         str_c(colnames_raw, "_split_read_annot"), colnames_raw)
  
  colnames(ERs_gr_split_read_annotated_tidy) <- 
    colnames_split_read_annot
  
  # convert df back to gr 
  ERs_gr_split_read_annotated_gr <-
    makeGRangesFromDataFrame(ERs_gr_split_read_annotated_tidy,
                             keep.extra.columns=T,
                             ignore.strand=FALSE,
                             seqinfo=NULL,
                             seqnames.field=c("seqnames", "seqname",
                                              "chromosome", "chrom",
                                              "chr", "chromosome_name",
                                              "seqid"),
                             start.field="start",
                             end.field= c("end", "stop"),
                             strand.field="strand",
                             starts.in.df.are.0based=FALSE)
  
  return(ERs_gr_split_read_annotated_gr)
  

}

check_split_read_connects_to_OMIM_gene <- function(ERs_gr, ensembl_grch38_v92_genes_txdb, OMIM_gene_start_stop_gr_marked_overlapping_genes){
  
  ensembl_grch38_v92_genes_transcripts <- transcripts(x = ensembl_grch38_v92_genes_txdb, columns = c("tx_name", "gene_id"))
  elementMetadata(ensembl_grch38_v92_genes_transcripts)["gene_id"] <- ensembl_grch38_v92_genes_transcripts$gene_id %>% unlist()
  
  # check that all OMIM genes are found within the ensembl txdb
  stopifnot(all(OMIM_gene_start_stop_gr_marked_overlapping_genes$gene_id %in% ensembl_grch38_v92_genes_transcripts$gene_id))
  
  elementMetadata(ERs_gr)["uniq_genes_split_read_annot"] <- 
    map_chr(.x = ERs_gr$uniq_transcripts_split_read_annot, .f = convert_transcript_to_gene, 
            ensembl_grch38_v92_genes_transcripts = ensembl_grch38_v92_genes_transcripts)
  
  elementMetadata(ERs_gr)["connect_to_OMIM_gene_split_read_annot"] <- 
    map_chr(.x = ERs_gr$uniq_genes_split_read_annot, .f = check_split_read_connect_OMIM_genes, 
            OMIM_gene_start_stop_gr_marked_overlapping_genes = OMIM_gene_start_stop_gr_marked_overlapping_genes)
    
  return(ERs_gr)
  
}

make_results_dir <- function(results_path, folder_name){
  
  results_dir_path <- str_c(results_path, "/", folder_name)
  
  if(!dir.exists(results_dir_path)){
    
    dir.create(results_dir_path)
    
  }else {
    
    print(str_c(results_dir_path, " directory already exists.."))
    
  }
  
  return(results_dir_path)
  
}

get_list_annotated_ERs <- function(path_to_annotated_ERs_dir, list_ERs_each_tissue_optimal_cut_off_maxgap){
  
  annotated_ERs_by_tissue_paths <- list.files(path_to_annotated_ERs_dir, recursive = T, full.names = T)
  
  list_ERs_w_annotation <- list()
  
  for(i in seq_along(annotated_ERs_by_tissue_paths)){
    
    print(str_c(Sys.time(), " - ", i))
    
    annotated_ERs_by_tissue_path <- annotated_ERs_by_tissue_paths[i]
    
    load(annotated_ERs_by_tissue_path)
    
    list_ERs_w_annotation[[i]] <- ERs_w_annotation
    
  }
  
  annotated_ERs_by_tissue_names <- 
    annotated_ERs_by_tissue_paths %>% 
    str_replace("results/.*/", "") %>% 
    str_replace("_annotated.rda", "") %>% 
    str_replace_all("_cutoff", " - cutoff") %>% 
    str_replace_all("_maxgap", " - maxgap")
  
  stopifnot(identical(annotated_ERs_by_tissue_names, names(list_ERs_each_tissue_optimal_cut_off_maxgap)))
      
  names(list_ERs_w_annotation) <- annotated_ERs_by_tissue_names
  
  return(list_ERs_w_annotation)
  
}

get_ERs_annotation_df <- function(list_ERs_w_annotation){
  
  ERs_w_annotation_all_tissues <- 
    data_frame()
  
  for(i in 1:length(list_ERs_w_annotation)){
    
    cut_off_tissue <- names(list_ERs_w_annotation[i])
    cut_off <- cut_off_tissue %>% str_replace(".*cutoff:", "") %>% str_replace(" -.*", "")
    maxgap <- cut_off_tissue %>% str_replace(".*maxgap:", "") 
    tissue <- cut_off_tissue %>% str_replace(" .*", "") 
    
    print(str_c(Sys.time(), " - ", i, " - ", tissue))
    
    ERs_w_annotation_tisue_filtered <- 
      list_ERs_w_annotation[[i]] %>% 
      as.data.frame() %>% 
      as_tibble()
    
    ERs_w_annotation_all_tissues <- 
      bind_rows(ERs_w_annotation_all_tissues, ERs_w_annotation_tisue_filtered)
    
  }
  
  return(ERs_w_annotation_all_tissues)
  
}

##### Second Level ##### 

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/get_region_coverage.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_txDb_from_gtf.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_genomic_state.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/convert_annot_count_table_to_region_annot.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/mark_overlapping_genes_gr.R")

read_format_gtex_split_read_table <- function(gtex_split_read_table_annotated_path){
  
  gtex_split_read_table_annotated <- read_delim(gtex_split_read_table_annotated_path, delim = ",", 
                                                col_types = cols(
                                                  .default = col_integer(),
                                                  junID = col_double(),
                                                  chr = col_character(), 
                                                  start = col_double(), 
                                                  stop = col_double(), 
                                                  strand = col_character(),
                                                  acceptor = col_character(),
                                                  donor = col_character(),
                                                  junction = col_character(),
                                                  precBoundDonor = col_logical(),
                                                  precBoundAcceptor = col_logical()
                                                ))
  
  # added to annotatER code
  # gtex_split_read_table_annotated <- 
  #   gtex_split_read_table_annotated %>% 
  #   mutate(acceptor = ifelse(is.na(acceptor), "", acceptor), 
  #          donor = ifelse(is.na(donor), "", donor), 
  #          junction = ifelse(is.na(junction), "", junction))
  
  return(gtex_split_read_table_annotated)
  
}

convert_transcript_to_gene <- function(uniq_transcripts, ensembl_grch38_v92_genes_transcripts){
  
  if(is.na(uniq_transcripts)){
    
    return(NA)
    
  }else{
    
    uniq_transcripts_vec <- str_split(uniq_transcripts, ", ") %>% unlist()
    ensembl_grch38_v92_genes_transcripts_filtered <- 
      ensembl_grch38_v92_genes_transcripts[ensembl_grch38_v92_genes_transcripts$tx_name %in% uniq_transcripts_vec]
    uniq_genes_vec <- 
      ensembl_grch38_v92_genes_transcripts_filtered$gene_id %>% 
      unique() %>% 
      str_c(collapse = ", ")
    
    return(uniq_genes_vec)
    
  }
  
}

check_split_read_connect_OMIM_genes <- function(uniq_genes, OMIM_gene_start_stop_gr_marked_overlapping_genes){
  
  if(is.na(uniq_genes)){
    
    return(NA)
    
  }else{
    
    uniq_genes_vec <- str_split(uniq_genes, ", ") %>% unlist()
    
    uniq_genes_vec_OMIM <- uniq_genes_vec[uniq_genes_vec %in% OMIM_gene_start_stop_gr_marked_overlapping_genes$gene_id]
    
    if(length(uniq_genes_vec_OMIM) == 0){
      
      return(NA)
      
    }else{
      
      uniq_genes_vec_OMIM_collapse <- uniq_genes_vec_OMIM %>% str_c(collapse = ", ")
      
      return(uniq_genes_vec_OMIM_collapse)
      
    }
    
  }
  
}

# Main ------------------------------------------------------------------------------------------------

# generate a list of genomic states from different reference annotations ready to annotate the ERs
txdb_sqlite_paths <- list.files(path = "/data/references/ensembl/txdb_sqlite/", recursive = T, full.names = T)
genomic_states_list <- generate_genomic_states_list(txdb_sqlite_paths)

# get df with split read table path and tissue
gtex_split_read_table_mean_cov_df <- 
  get_gtex_split_read_table_mean_cov_df()

list_ERs_w_annotation <- list()

ensembl_grch38_v92_genes_txdb <- 
  generate_txDb_from_gtf(gtf_gff3_path = "/data/references/ensembl/gtf_gff3/v92/Homo_sapiens.GRCh38.92.gtf", 
                         output_path = "/data/references/ensembl/txdb_sqlite/v92/ensembl_grch38_v92_txdb.sqlite",
                         seq_levels_to_keep = c(1:22, "X", "Y", "MT"), genome_build = "hg38")

for(i in 1:length(list_ERs_each_tissue_optimal_cut_off_maxgap)){
  
  cut_off_tissue <- names(list_ERs_each_tissue_optimal_cut_off_maxgap[i])
  cut_off <- cut_off_tissue %>% str_replace(".*cutoff:", "") %>% str_replace(" -.*", "")
  maxgap <- cut_off_tissue %>% str_replace(".*maxgap:", "") 
  tissue <- cut_off_tissue %>% str_replace(" .*", "") 
  
  results_path_tissue <- make_results_dir(results_path = "results/annotate_ERs/by_tissue/", tissue)
  
  ERs_one_tissue_optimal_cut_off_no_scaffold <- 
    list_ERs_each_tissue_optimal_cut_off_maxgap[[i]] %>% 
    keepSeqlevels(str_c("chr", c(1:22, "X", "Y")), pruning.mode = "coarse") %>% 
    sortSeqlevels() %>% 
    sort()
  
  print(str_c(Sys.time(), " - ", i, " - ", names(list_ERs_each_tissue_optimal_cut_off_maxgap[i])))
  
  print(str_c(Sys.time(), " - adding ER index"))
  
  elementMetadata(ERs_one_tissue_optimal_cut_off_no_scaffold)[["ER_index"]] <- 1:length(ERs_one_tissue_optimal_cut_off_no_scaffold)
  
  print(str_c(Sys.time(), " - getting quantification of coverage for ERs"))
  
  ERs_w_annotation <- 
    get_coverage_for_GTEx_ERs(ERs_gr = ERs_one_tissue_optimal_cut_off_no_scaffold, tissue_to_filter = tissue, 
                              gtex_split_read_table_mean_cov_df = gtex_split_read_table_mean_cov_df)
  
  print(str_c(Sys.time(), " - annotating ERs by annotation feature"))
  
  ERs_w_annotation <- 
    annotate_ERs_w_genomic_state(ERs_gr = ERs_w_annotation, genomic_states_list = genomic_states_list)
  
  print(str_c(Sys.time(), " - annotating ERs by nearest or overlapping genes"))
  
  ERs_w_annotation <- 
    annotate_ERs_with_overlap_and_dis_to_genes(ERs_gr = ERs_w_annotation, 
                                              genes_gr = genes(ensembl_grch38_v92_genes_txdb), 
                                              column_name_prefix = "any_gene_v92_")
  
  ERs_w_annotation <- 
    annotate_ERs_with_overlap_and_dis_to_genes(ERs_gr = ERs_w_annotation, 
                                              genes_gr = OMIM_gene_start_stop_gr_marked_overlapping_genes, 
                                              column_name_prefix = "any_OMIM_gene_v92_")
  
  ERs_w_annotation <- 
    annotate_ERs_with_overlap_and_dis_to_genes(ERs_gr = ERs_w_annotation, 
                                               genes_gr = OMIM_gene_start_stop_gr_marked_overlapping_genes
                                               [OMIM_gene_start_stop_gr_marked_overlapping_genes$overlap_any_other_gene_v92 == F], 
                                               column_name_prefix = "non_overlapping_OMIM_gene_v92_")
  
  print(str_c(Sys.time(), " - annotating ERs using split reads"))
  
  ERs_w_annotation <- 
    annotate_ERs_w_split_reads(ERs_gr = ERs_w_annotation, tissue_to_filter = tissue, 
                               gtex_split_read_table_mean_cov_df = gtex_split_read_table_mean_cov_df)
  
  print(str_c(Sys.time(), " - checking ERs that link to OMIM genes through split reads"))
  
  ERs_w_annotation <- 
    check_split_read_connects_to_OMIM_gene(ERs_gr = ERs_w_annotation, ensembl_grch38_v92_genes_txdb = ensembl_grch38_v92_genes_txdb, 
                                           OMIM_gene_start_stop_gr_marked_overlapping_genes = OMIM_gene_start_stop_gr_marked_overlapping_genes)
  
  print(str_c(Sys.time(), " - annotating ERs with constratint and conservation scores"))
  
  ERs_w_annotation <- get_constraint_score_for_regions(region_details_gr = ERs_w_annotation, 
                                                       constraint_gr = CDTS_percentile_N7794_unrelated_all_chrs_gr)
  
  ERs_w_annotation <- get_conservation_score_for_regions(conserv_score_to_load = "phast_cons_7", 
                                                         gr = ERs_w_annotation)
  
  ERs_w_annotation <- get_conservation_score_for_regions(conserv_score_to_load = "phast_cons_100", 
                                                         gr = ERs_w_annotation)
  
  print(str_c(Sys.time(), " - adding tissue"))

  elementMetadata(ERs_w_annotation)[["tissue"]] <- tissue

  list_ERs_w_annotation[[i]] <-
    ERs_w_annotation
  
  save(ERs_w_annotation, file = str_c(results_path_tissue, "/", tissue, "_cutoff:", cut_off, "_maxgap:", maxgap, "_annotated.rda"))
  
}

names(list_ERs_w_annotation) <-
  names(list_ERs_each_tissue_optimal_cut_off_maxgap)

# only use if have parrallelised by tissue to save time - just collects all annotated ERs into one object
list_ERs_w_annotation <-
  get_list_annotated_ERs(path_to_annotated_ERs_dir = "results/annotate_ERs/by_tissue/", list_ERs_each_tissue_optimal_cut_off_maxgap)

# merge all annotations into one df and add a column to show where
ERs_w_annotation_all_tissues <-
  get_ERs_annotation_df(list_ERs_w_annotation)

# Save data -------------------------------------------------------------------------------------------

save(list_ERs_w_annotation,
     file = "results/annotate_ERs/ERs_optimised_cut_off_max_gap_all_tissues_w_annot_list.rda")

save(ERs_w_annotation_all_tissues, file = "results/annotate_ERs/ERs_optimised_cut_off_max_gap_all_tissues_w_annot_df.rda")

