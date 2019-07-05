library(tidyverse)
library(stringr)
library(httr) # package for R to use API
library(jsonlite) # convert downloaded json formatted txt to list in R

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- Sys.getenv("OMIM_wd")
setwd(OMIM_wd)

# Functions -------------------------------------------------------------------------------------------
  
get_pheno_with_molecular_basis <- function(omim_api_url){
  
  # search path for phenotypes with a "#" prefix
  # limit to 100000 to return all (expecting ~5000)
  search_prefix_number_sign_path <-
    "api/entry/search?search=prefix%3A%23&start=0&limit=100000&format=json&apiKey=XXXXXXXX"

  # download list of phenotype MIM numbers using API
  pheno_mol_known_raw <- GET(url = omim_api_url, path = search_prefix_number_sign_path)
  pheno_mol_known_list <- convert_raw_API_data_to_list(pheno_mol_known_raw$content)
  
  MIM_number_pheno <- pheno_mol_known_list[["omim"]][["searchResponse"]][["entryList"]][["entry"]][["mimNumber"]]
  omim_pheno_preferred_title <- pheno_mol_known_list[["omim"]][["searchResponse"]][["entryList"]][["entry"]][["titles"]][["preferredTitle"]]
  
  pheno_mol_known_df <- 
    data_frame(MIM_number_pheno = MIM_number_pheno, 
               omim_pheno_preferred_title = omim_pheno_preferred_title) %>% 
    mutate(index = seq_along(MIM_number_pheno))
  
  return(pheno_mol_known_df)
  
}

download_details_by_MIM_num_pheno <- function(omim_api_url, MIM_number_pheno, omim_pheno_preferred_title, data_to_download){
  
  # Get all details available for phenotype entry
  MIM_num_pheno_path <- str_c("api/entry?mimNumber=", MIM_number_pheno, 
                              "&include=all&format=json&apiKey=CyS6AQ6JTjOcWE1W27WHuQ")
  
  MIM_num_pheno_details_raw <- GET(url = omim_api_url, path = MIM_num_pheno_path)
  MIM_num_pheno_details_list <- convert_raw_API_data_to_list(MIM_num_pheno_details_raw$content)
  
  omim_pheno_preferred_title <- 
    MIM_num_pheno_details_list[["omim"]][["entryList"]][["entry"]][["titles"]][["preferredTitle"]]
  
  PM_date_created <- 
    MIM_num_pheno_details_list[["omim"]][["entryList"]][["entry"]][["dateCreated"]]
  
  PM_date_updated <- 
    MIM_num_pheno_details_list[["omim"]][["entryList"]][["entry"]][["dateUpdated"]]
  
  # Extract the clinical synopsis or phenotype map
  if(data_to_download == "clinical_synopsis"){
    
    OMIM_pheno_details <- 
      MIM_num_pheno_details_list[["omim"]][["entryList"]][["entry"]][["clinicalSynopsis"]] 
    
  }else if(data_to_download == "phenotype_map"){
    
    OMIM_pheno_details <- 
      MIM_num_pheno_details_list[["omim"]][["entryList"]][["entry"]][["phenotypeMapList"]][[1]][["phenotypeMap"]] 
    
  }
  
  # check if there details exist (if none, then previous step will return NULL)
  # add variable to describe if details exist
  if(!is.null(OMIM_pheno_details)){
    
    # sometimes clinical synopsis will be stored as a complete data.frame
    # other times it may have further data.frames within it
    # this function unravels all data into one complete data frame
    OMIM_pheno_details_tibble <- 
      OMIM_pheno_details %>% convert_data_to_tibble()
    
    OMIM_pheno_details_with_MIM_num <- 
      OMIM_pheno_details_tibble %>% 
      mutate(MIM_number_pheno = MIM_number_pheno,
             omim_pheno_preferred_title = omim_pheno_preferred_title, 
             PM_date_created = PM_date_created, 
             PM_date_updated = PM_date_updated, 
             details_exist = TRUE)
    
  }else{
    
    OMIM_pheno_details_with_MIM_num <-
      data_frame(MIM_number_pheno = MIM_number_pheno,
                 omim_pheno_preferred_title = omim_pheno_preferred_title, 
                 PM_date_created = PM_date_created, 
                 PM_date_updated = PM_date_updated,
                 details_exist = FALSE)
    
  }
  
  return(OMIM_pheno_details_with_MIM_num)
  
}

convert_raw_API_data_to_list <- function(x){
  
  # convert raw binary data to char then to R list 
  x_char <- rawToChar(x)
  x_list <- fromJSON(x_char)
  
  return(x_list)
  
}

convert_data_to_tibble <- function(x){
    
  for(i in seq_along(x)){
    
    # if the element within position i is data.frame
    if(is.data.frame(x[[i]])){
      
      if(i == 1){
        
        # if first element, store as overall df 
        # this bypasses problems of cbinding dataframes of unequal rows, since empty data frame will have 0 rows
        # note the [[i]] subsetting has two brackets to extract the element at i 
        x_df <- x[[i]]
        
      }else{
        
        # otherwise cbind to existing overall df
        x_df <- cbind(x_df, x[[i]])
        
      }
      
    # if the element within position i is vector
    }else if(is.vector(x[[i]])){
      
      if(i == 1){
        
        # same as above but [i] subsets the dataframe containing element i, instead of element i itself
        x_df <- x[i]
        
      }else{
        
        x_df <- bind_cols(x_df, x[i])
        
      }
      
    }

  }
  
  # convert to tibble for easier tidyverse manipulation
  x_tibble <- x_df %>% as_tibble()
  
  return(x_tibble)
  
}

# Main ------------------------------------------------------------------------------------------------

omim_api_url <- "http://api.omim.org"

current_date <- Sys.Date() %>% str_replace_all("-", "_")

# get all phenotype MIM numbers 
pheno_mol_known_df <- get_pheno_with_molecular_basis(omim_api_url)

clinical_synopsis_all <- data_frame()
phenotype_map_all <- data_frame()

# loop across all MIM numbers and extract data from OMIM using API, appending to overall dataframe
for(i in seq_along(pheno_mol_known_df$MIM_number_pheno)){
  
  MIM_number_pheno <- pheno_mol_known_df$MIM_number_pheno[i]
  
  print(str_c(i, " - ", MIM_number_pheno))
  
  clinical_synopsis <- download_details_by_MIM_num_pheno(omim_api_url, MIM_number_pheno, omim_pheno_preferred_title, 
                                                  data_to_download = "clinical_synopsis")
  
  phenotype_map <- download_details_by_MIM_num_pheno(omim_api_url, MIM_number_pheno, omim_pheno_preferred_title, 
                                              data_to_download = "phenotype_map")
  
  clinical_synopsis_all <- bind_rows(clinical_synopsis_all, clinical_synopsis)
  phenotype_map_all <- bind_rows(phenotype_map_all, phenotype_map)
  
}

# Save data -------------------------------------------------------------------------------------------

write_delim(clinical_synopsis_all, 
            str_c("results/download_tidy_OMIM_data//OMIM_clinical_synopsis_prefix_num_sign_raw_", current_date, ".csv"), 
            delim = ",")

write_delim(phenotype_map_all, 
            str_c("results/download_tidy_OMIM_data/OMIM_phenotype_map_prefix_num_sign_raw_", current_date, ".csv"),
            delim = ",")
