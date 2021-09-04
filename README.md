## Background

This repository contains analysis code for the publication found [here](https://www.science.org/doi/10.1126/sciadv.aay8299), which aimed at improving the annotation of disease-relevant genes using RNA-sequencing data. The accompanying web-based tool [vizER](https://snca.atica.um.es/browser/app/vizER) can be used to visualise individual genes of interest for evidence of incomplete annotation. 

## Citation

If you use any code part of this repository please cite the Science Advances publication: DOI [10.1126/sciadv.aay8299](https://www.science.org/doi/10.1126/sciadv.aay8299).

## Code contents

| Directory | Description |
| --------- | ----------- |
| [analyse_ER_annotation](analyse_ER_annotation) | ER related analyisis including number of quantifying OMIM gene re-annotation and total ER Mb across annotation features. Validation of ERs across Ensembl versions and within an independent dataset |
| [annotate_ERs](annotate_ERs) | Annotating ERs with metrics such as association to genes through junctions, annotation features, conservation and constraint |
| [check_protein_coding_potential](check_protein_coding_potential) | Checking protein potential of ERs |
| [complex_disorders](complex_disorders) | Re-annotation of GWAS hits from [STOPGAP](https://pubmed.ncbi.nlm.nih.gov/28472345-stopgap-a-database-for-systematic-target-opportunity-assessment-by-genetic-association-predictions/) |
| [download_tidy_OMIM_data](download_tidy_OMIM_data) | Download details of Mendelian disease genes via [OMIM](https://omim.org/) API |
| [export_ER_details](export_ER_details) | Formatting ER details for publication |
| [generate_ERs_varying_cut_offs_maxgaps_GTEx_tissues](generate_ERs_varying_cut_offs_maxgaps_GTEx_tissues) | Using [derfinder](https://bioconductor.org/packages/release/bioc/html/derfinder.html) to define tissue-specific expressed regions (ERs) for each GTEx tissue* |
| [generate_randomised_intron_inter_regions](generate_randomised_intron_inter_regions) | Generating tissue-specific randomised length-matched regions |
| [GTEx_split_read_reformatting](GTEx_split_read_reformatting) | Re-format the raw GTEx junction data dowloaded from [recount2](https://jhubiostatistics.shinyapps.io/recount/) for input into [annotatER](https://github.com/SebGuelfi/annotatER) |
| [optimising_derfinder_cutoff](optimising_derfinder_cutoff) | Optimising the definitions of ERs using a gold-standard set of non-overlapping exons* |

*These elements of the pipeline have been wrapped into an R package that can be found [here](https://github.com/dzhang32/ODER). 
