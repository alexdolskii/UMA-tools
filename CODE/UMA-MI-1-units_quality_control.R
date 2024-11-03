# Data Quality Control Script for Biological Data Types: GS, pFAK, PALLD, pSMAD
# Author: [Your Name]
# Date: [Date]
#
# This script performs quality control on a specified data file by filtering out rows
# where the Matrix/WIM (Whole Image) area ratio is below a user-defined threshold.
# The script prompts the user for the file path, data type, and threshold value.
# It processes the data accordingly and saves the results in a uniquely named folder
# within a 'QC' directory located in the same directory as the script.
# The 'QC' directory contains:
# - Subfolders for each run named with date, time, filename, marker, and threshold
#   - Filtered good data
#   - A separate file with the discarded data
#   - A log file including the threshold value
#
# Requirements:
# - R version >= 3.6.0
# - Packages: dplyr, stringr
#
# Usage:
# - Run the script in R or RStudio.
# - When prompted, enter the path to the data file, the data type (GS, pFAK, PALLD, pSMAD),
#   and the threshold value.
#
# Example:
# source("quality_control_script.R")
#
# License: MIT License

# Load necessary libraries
library(dplyr)
library(stringr)
library(readxl)

pattern_fib <- "^Cell\\.\\.[^WIM]\\S+[Mm]atrix\\S+Total\\.Area_Sum"
pattern_WIM <- "^Cell\\.\\.WIM\\S+[Mm]atrix\\S+Total\\.Area_Sum"
pattern_Nuclei_counts <- "^Cell\\.\\.\\S+nuc\\.\\S+msk\\S+Features\\.Count_Sum"
pattern_fib_intensity <- "^Cell\\.\\.[^WIM]\\S+[Mm]atrix\\S+Integrated\\.Intensity_Sum"
pattern_intensity <- "^Cell\\.\\.\\S+_obj\\.?msk_Integrated\\.Intensity_Sum"

# Function to process the file
process_file <- function(file_path, data_type_key, data_type, threshold) {
  # Read the file
  if (endsWith(file_path, '.txt')) {
    prel_data <- read.table(file_path, 
                          header = TRUE, 
                          sep = "\t", 
                          stringsAsFactors = FALSE)
  } else if (endsWith(file_path, '.xlsx')) {
    prel_data <- read_excel(file_path)
    # Rename columns in excel file to match with patterns
    colnames(prel_data) <- gsub("\\s", ".", 
                                gsub(":", ".", 
                                     gsub("\\(", ".",
                                          gsub("\\)", ".", colnames(prel_data)))))
  } else {
    stop('Please use file with .txt or .xlsx extensions')
  }
  
  # Check if data_type argument corresponds to file contents
  coulumn_substring <- c('GS'='GS', 'pFAK'='pFAK', 'PALLD'='iso3', 'pSMAD'='psmad')
  if (!any(grepl(coulumn_substring[data_type], colnames(prel_data)))){
    stop('The marker name does not match the file content!')
  }
  
  # Find the names of columns that correspond to our patterns
  fib <- str_subset(colnames(prel_data), regex(pattern_fib))
  wim <- str_subset(colnames(prel_data), regex(pattern_WIM))
  nc <- str_subset(colnames(prel_data), regex(pattern_Nuclei_counts))
  fib_int <- str_subset(colnames(prel_data), regex(pattern_fib_intensity))
  mar_int <- str_subset(colnames(prel_data), regex(pattern_intensity))
  
  # Filter our df and choose only features that are necessary for the following analysis
  data <- prel_data[c('Plate.ID', 
                    'Well.Name', 
                    'Site.ID', 
                    'MEASUREMENT.SET.ID', 
                    fib, 
                    wim, 
                    nc, 
                    fib_int, 
                    mar_int
  )]
  
  # Make names of our df readable
  names(data) <- c('Plate.ID', 
                   'Well.Name', 
                   'Site.ID', 
                   'MEASUREMENT.SET.ID', 
                   paste0(data_type, '_fibronectin_mask'), 
                   paste0(data_type, '_wim_mask'), 
                   paste0(data_type, '_nuclei_counts'), 
                   paste0(data_type, '_fibronectin_integrated_intensity_sum'),
                   paste0(data_type, '_marker_integrated_intensity_sum'))
  
  #Check if multiple or no columns are matched
  if (any(length(fib) == 0,
         length(wim) == 0,
         length(nc) == 0,
         length(fib_int) == 0,
         length(mar_int) == 0)) {
    message <- sprintf("Required columns not found in the dataset for file: %s", basename(file_path))
    return(message)
  } else if (any(length(fib) > 1,
                length(wim) > 1,
                length(nc) > 1,
                length(fib_int) > 1,
                length(mar_int) > 1)) {
    message <- sprintf("Multiple columns matched the pattern in file: %s. Please check columns!", basename(file_path))
    return(message)
  }
  
  # Add NEW_NAMES column if it doesn't exist
  if (!"NEW_NAMES" %in% colnames(data)) {
    data$NEW_NAMES <- paste(data$Plate.ID, data$Well.Name, data$Site.ID, data$MEASUREMENT.SET.ID)
  }
  
  # Calculate the Matrix/WIM Area ratio
  ratio_values <- ((data[[paste0(data_type, '_fibronectin_mask')]] 
                    / data[[paste0(data_type, '_wim_mask')]]) * 100)
  
  # Add the ratio to the data
  data <- data %>%
    mutate(Matrix_WIM_Area_ratio = ratio_values)
  
  # Filter data based on the threshold
  data_filtered <- data %>%
    filter(!is.na(Matrix_WIM_Area_ratio) & Matrix_WIM_Area_ratio > threshold)
  
  # Data that did not meet the threshold
  data_discarded <- data %>%
    filter(is.na(Matrix_WIM_Area_ratio) | Matrix_WIM_Area_ratio <= threshold)
  
  # Calculate statistics
  total_rows <- nrow(data)
  removed_rows <- nrow(data_discarded)
  percent_removed <- (removed_rows / total_rows) * 100
  names_removed <- data_discarded$NEW_NAMES
  
  
  # Create 'QC' directory in the same directory as the script
  script_dir <- tryCatch({
    dirname(normalizePath(sys.frames()[[1]]$ofile))
  }, error = function(e) {
    getwd()
  })
  qc_dir <- file.path(script_dir, "QC")
  dir.create(qc_dir, showWarnings = FALSE)
  
  # Create output directory within 'QC' with date, time, filename, marker, and threshold
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  base_filename <- tools::file_path_sans_ext(basename(file_path))
  output_dir_name <- paste0("QC_", base_filename, "_", data_type, "_", "thresh", threshold, "_", timestamp)
  output_dir <- file.path(qc_dir, output_dir_name)
  dir.create(output_dir, showWarnings = FALSE)
  
  # Save the filtered data
  output_file_name <- strsplit(basename(file_path), split = "\\.")[[1]][1]
  output_file_filtered <- file.path(output_dir, paste0('filtered_', output_file_name, '.csv'))
  write.csv(data_filtered, file = output_file_filtered, row.names = FALSE, quote = FALSE)
  
  # Save the discarded data
  output_file_discarded <- file.path(output_dir, paste0('discarded_', output_file_name, '.csv'))
  write.csv(data_discarded, file = output_file_discarded, row.names = FALSE, quote = FALSE)
  
  # Prepare log information
  log_info <- sprintf('File: %s\nData Type: %s\nThreshold: %.2f%%\n%d rows removed (%.2f%%) because Matrix/WIM ratio < %.2f%%.\nRemoved rows: %s\n', 
                      basename(file_path), data_type, threshold, removed_rows, percent_removed, threshold, toString(names_removed))
  
  # Save log information to a file in the output directory
  log_file <- file.path(output_dir, "Data_Quality_Log.txt")
  writeLines(log_info, con = log_file)
  
  message(sprintf("Processing complete. Results saved in folder: %s", output_dir))
}


# #Ask for user input
# file_path <- readline(prompt = "Enter the path to the data file: ")
# data_type <- readline(prompt = "Enter the data type (GS, pFAK, PALLD, pSMAD): ")
# threshold_input <- readline(prompt = "Enter the threshold value: ")

#Kate_example

file_path <- '/Users/ekaterinashitik/UMA-tools/DATA/for-umi-ma/data_pt_cu_far_pl_09242024/excel-raw-data/Cukierman_TL_GS_ptCufar_pl09242024_pl2298_SITE.xlsx'
data_type <- 'GS'
threshold_input <- '50'


# Convert threshold to numeric
threshold <- as.numeric(threshold_input)

# Validate data type
valid_data_types <- c("GS", "pFAK", "PALLD", "pSMAD")
valid_data_types_lower <- tolower(valid_data_types)
data_type_lower <- tolower(data_type)
if (!data_type_lower %in% valid_data_types_lower) {
  stop("Invalid data type specified. Please choose from 'GS', 'pFAK', 'PALLD', or 'pSMAD'.")
}

data_type_key <- data_type_lower

# Process file path to handle backslashes and quotes
file_path <- gsub('\"', '', file_path)
file_path <- gsub('\\\\', '/', file_path)

# Check if file exists
if (!file.exists(file_path)) {
  stop("The specified data file does not exist.")
}

# Process the file
process_file(file_path, data_type_key, data_type, threshold)
