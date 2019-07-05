library(tidyverse)
library(stringr)
library(RColorBrewer)
library(ggpubr)
library(SummarizedExperiment)
library(forcats)
library(data.table)
library(annotatER)
# library(devtools)
# install("~/packages/annotatER/")

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

GTEx_split_read_counts_paths <- list.files(path = "/data/recount/GTEx_SRP012682/gtex_junction_coverage_by_tissue/datatable/", 
                                           full.names = T, recursive = T)

gtf_path <- "/data/references/ensembl/gtf_gff3/v92/Homo_sapiens.GRCh38.92.gtf"

# Main ------------------------------------------------------------------------------------------------

##### annotate gtex split reads #####

# format split reads info
splitReadsInfo_2 <- fread("/data/recount/GTEx_SRP012682/SRP012682.junction_id_with_transcripts.bed")
splitReadsInfo_2[,"junID" := tstrsplit(V4, "|", fixed=TRUE,keep=1)]
# setDF(splitReadsInfo)
# save(splitReadsInfo, file = "/data/recount/GTEx_SRP012682/splitReadsInfo_formatted.rda")

load("/data/recount/GTEx_SRP012682/splitReadsInfo_formatted.rda")

for(i in seq_along(GTEx_split_read_counts_paths)){
  
  GTEx_split_read_counts_path <- GTEx_split_read_counts_paths[i]
  tissue <- 
    GTEx_split_read_counts_path %>% 
    str_replace("/.*/", "") %>% 
    str_replace("\\.csv", "") %>% 
    str_replace_all("-", "_") %>% 
    str_to_lower()
  
  print(str_c(i, " - ", tissue))
  
  GTEx_split_read_count <- read_delim(GTEx_split_read_counts_path, delim = ",", 
                                      col_types = cols(.default = "i", juncID = "n")) %>% 
    dplyr::rename(junID = juncID)
  
  GTEx_split_read_count <- GTEx_split_read_count %>% data.frame()
  
  GTEx_split_read_table <- filterSplitReads(splitReadCounts = GTEx_split_read_count, 
                                            splitReadsInfo = splitReadsInfo, minSamples = 5, minCounts = 1) 
  
  # convert UCSC to ensembl 
  GTEx_split_read_table$start <- GTEx_split_read_table$start +1
  GTEx_split_read_table$stop <- GTEx_split_read_table$stop +1
  
  GTEx_split_read_table_annotated <- annotateSplitReads(GTFPath = gtf_path, 
                                                        splitReadTable = GTEx_split_read_table)
  
  gtex_split_read_table_annotated_only_junc_coverage <-
    GTEx_split_read_table_annotated %>% 
    dplyr::select(junID, chr, start, stop, strand, countsSamples, acceptor, donor, junction, precBoundDonor, precBoundAcceptor) %>%
    filter(!is.na(junID))
  
  save(gtex_split_read_table_annotated_only_junc_coverage, 
       file = str_c("/data/recount/GTEx_SRP012682/gtex_split_read_table_annotated_rda/", tissue, "_split_read_table_annotated.rda"))
  
  rm(GTEx_split_read_count)
  rm(GTEx_split_read_table)
  rm(GTEx_split_read_table_annotated)
  
}



