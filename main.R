#!/usr/bin/Rscript
## Author: Swathy Selvakumar


#### Bioconductor ####
# it is standard among R packages to define libraries and packages at the 
# beginning of a script. Also note that a package should NOT be installed every 
# time a script runs.
# The bioconductor repository has installation instructions for biomaRt: 
# https://bioconductor.org/install/

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!require("biomaRt", quietly = TRUE)){
  BiocManager::install("biomaRt")
}
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(tidyverse))

#### Loading and processing data ####
#' Load Expression Data
#'
#' @param filepath A text string of the full filepath to the file to load.
#'
#' @return A tibble containing the data loaded from the CSV in `filepath`. 
#' 
#' @details Note that not all CSVs are created equal, and there are often cases where 
#' the data will not load in correctly on the first try. You may want to write this functon to 
#' adjust the CSV file being loaded into this assignment so that it can be formed into a 
#' tibble correctly.
#'
#' @examples 
#' `data <- load_expression('/project/bf528/project_1/data/example_intensity_data.csv')`
load_expression <- function(filepath) {
  data <- read_csv(filepath, show_col_types = FALSE)
  
  if (startsWith(colnames(data)[1], "...")) {
    colnames(data)[1] <- "probe_id"
  }
  
  return(data)
}

#' Filter 15% of the gene expression values.
#'
#' @param tibble A tibble of expression values, rows by probe and columns by sample.
#'
#' @return A tibble of affymetrix probe names from the input expression data tibble. 
#' These names match the rows with 15% or more of the expression values about log2(15).
#' 
#'
#' @examples `samples <- filter_15(data_tib)`
#' `> str(samples)`
#' `tibble [40,158 × 1] (S3: tbl_df/tbl/data.frame)`
#' `$ probe: chr [1:40158] "1007_s_at" "1053_at" "117_at" "121_at" ...`
filter_15 <- function(tibble){
  # 1. Separate the probe names from the numeric data
  # Assumes the first column contains the probe IDs
  probe_names <- tibble[[1]]
  expr_matrix <- as.matrix(tibble[, -1])
  
  # 2. Define the threshold
  threshold <- log2(15)
  
  # 3. Calculate the proportion of samples per row exceeding the threshold
  # rowSums returns how many samples > log2(15)
  # Dividing by ncol(expr_matrix) gives the percentage/proportion
  passing_indices <- rowSums(expr_matrix > threshold) / ncol(expr_matrix) >= 0.15
  
  # 4. Return as a tibble with a single column named 'probe'
  filtered_probes <- tibble(probe = probe_names[passing_indices])
  
  return(filtered_probes)
}

#### Gene name conversion ####

#' Convert affymetrix array names into hgnc_symbol IDs using biomaRt. Inputs and 
#' outputs will likely not be the same size.
#'
#' @param affy_tib A single column tibble of strings containing array names.
#'
#' @return A 2 column tibble that contains affy IDs in the first column,
#' and their corresponding HGNC gene ID in the second column. Note that not all affy IDs 
#' will necessarily correspond with a gene ID, and one gene may have multiple affy IDs.
#' 
#' @details Connecting to ensembl via biomaRt can be...hit or miss...so you may 
#' want to check if data was correctly returned (or if it was just empty). The 
#' `getBM()` function may not accept a tibble, so you might need to convert your 
#' input into a flat vector.
#'
#' @examples 
#' `> affy_to_hgnc(tibble(c('202860_at', '1553551_s_at')))`
#' `affy_hg_u133_plus_2 hgnc_symbol`
#' `1        1553551_s_at      MT-ND1`
#' `2        1553551_s_at       MT-TI`
#' `3        1553551_s_at       MT-TM`
#' `4        1553551_s_at      MT-ND2`
#' `5           202860_at     DENND4B`
affy_to_hgnc <- function(affy_vector) {
  # 1. Convert tibble column to a flat character vector
  # Using [[1]] ensures we get the vector regardless of column name
  affy_ids <- affy_vector[[1]]
  
  # 2. Connect to Ensembl (using the Human dataset)
  # Clear cache
  biomaRt::biomartCacheClear()
  mart <- useMart("ensembl", 
                          dataset = "hsapiens_gene_ensembl", 
                          host = "https://www.ensembl.org")
  
  # 3. Perform the query
  # filters: the type of ID we are providing
  # attributes: the data we want to retrieve
  res <- getBM(
    attributes = c('affy_hg_u133_plus_2', 'hgnc_symbol'),
    filters = 'affy_hg_u133_plus_2',
    values = affy_ids,
    mart = mart
  )
  
  # 4. Convert result to tibble and handle potential empty returns
  res_tibble <- as_tibble(res)
  
  # Optional: Filter out empty HGNC symbols if desired
  # res_tibble <- res_tibble %>% filter(hgnc_symbol != "")
  
  return(res_tibble)
}

#' Reduce a tibble of expression data to only the rows in good_genes or bad_genes.
#'
#' @param expr_tibble A tibble holding the expression data, each row corresponding to
#' one affymetrix probe ID and each column to a sample.
#' @param names_ids A two column tibble that associates affy IDs with HGNC gene IDs. 
#' Generated `with affy_to_hgnc()`.
#' @param good_genes A list of gene names stored as a vector of strings.
#' @param bad_genes A list of gene names stored as a vector of strings.
#'
#' @return A tibble with two additional columns added:
#' 1. HGNC gene IDs 
#' 2. Does the gene is this row fall into "good" or "bad" genes?
#' This tibble should be reduced to only rows present in good or bad genes. All
#' other rows can be discarded.
#' 
#' @details In order to plot only our genes of interest, we need to rearrange our 
#' data to include only the elements we want to see. We also want to add to columns, 
#' one that associates the probeids with the HGNC gene name, and one that says if 
#' that gene is in the good or bad sets of genes.
#'
#' @examples 
#' `plot_tibble <- reduce_data(expr_tibble = expr, names_ids = sample_names,`
#' `                           goodGenes, badGenes)`
#' `> head(plot_tibble)`
#' `A tibble: 6 × 38`
#' `  probeids    hgnc    gene_set    GSM972389 ...`
#' `  <chr>       <chr>   <chr>       <dbl>     ...`
#' `1 202860_at   DENND4B good        7.16      ...`
#' `2 204340_at   TMEM187 good        6.40      ...`
reduce_data <- function(expr_tibble, names_ids, good_genes, bad_genes){
  # 1. Identify existing column names dynamically
  first_col_expr <- colnames(expr_tibble)[1]
  map_affy_col <- colnames(names_ids)[1]
  map_hgnc_col <- colnames(names_ids)[2]
  
  # 2. Join and transform
  reduced <- expr_tibble %>%
    # Join the expression data with the mapping table
    inner_join(names_ids, by = setNames(map_affy_col, first_col_expr)) %>%
    # Rename the HGNC column to the required 'hgnc'
    rename(hgnc = !!sym(map_hgnc_col)) %>%
    # Categorize the genes
    mutate(
      gene_set = case_when(
        hgnc %in% good_genes ~ "good",
        hgnc %in% bad_genes ~ "bad",
        TRUE ~ NA_character_
      )
    ) %>%
    # Drop genes not in our lists
    filter(!is.na(gene_set)) %>%
    # IMPORTANT: Ensure the first column is named 'probeids' for the next function
    select(probeids = !!sym(first_col_expr), hgnc, gene_set, everything())
  
  return(reduced)
}

#' Convert a wide format tibble to long for easy plotting
#'
#' @param tibble A tibble of data in wide format. Specifically, it should be 
#' operating on the reduced tibble created by the previous function
#'
#' @return A tibble properly converted from wide to long format, the old sample 
#' columns should now be contained within a single column called "sample"
#' @export
#'
#' @examples
convert_to_long <- function(tibble) {
  # This will pivot all columns EXCEPT probeids, hgnc, and gene_set
  long_tibble <- tibble %>%
    pivot_longer(
      cols = -c(probeids, hgnc, gene_set), 
      names_to = "sample", 
      values_to = "expression"
    )
  
  return(long_tibble)
}

