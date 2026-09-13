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

# Define patterns for each data type
patterns <- list(
  gs = list(
    matrix_pattern = "^Cell\\.\\.GS_\\.Matrix\\.objmsk_Total\\.Area_Sum.*",
    wim_pattern = "^Cell\\.\\.WIM_GS_matrix_Total\\.Area_Sum.*"
  ),
  pfak = list(
    matrix_pattern = "^Cell\\.\\.pFAK\\.Matrix\\.objmsk_Total\\.Area_Sum.*",
    wim_pattern = "^Cell\\.\\.WIM_pFAK_Matrix\\.Area_Sum.*"
  ),
  palld = list(
    # Updated patterns for PALLD to match your data file
    matrix_pattern = "^Cell\\.\\.iso3_\\.Matrix\\.objmsk_Area_Sum.*",
    wim_pattern = "^Cell\\.\\.WIM_iso3_matrix_Area_Sum.*"
  ),
  psmad = list(
    matrix_pattern = "^Cell\\.\\.psmad_matrix\\.objmsk_Total\\.Area_Sum.*",
    wim_pattern = "^Cell\\.\\.WIM_psmad_matrix_Total\\.Area_Sum.*"
  )
)

# Function to process the file
process_file <- function(file_path, data_type_key, data_type, threshold) {
  # Read the file
  data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Add NEW_NAMES column if it doesn't exist
  if (!"NEW_NAMES" %in% colnames(data)) {
    data$NEW_NAMES <- paste(data$Plate.ID, data$Well.Name, data$Site.ID, data$MEASUREMENT.SET.ID)
  }
  
  # Get patterns for the specified data type
  pattern_info <- patterns[[data_type_key]]
  matrix_pattern <- pattern_info$matrix_pattern
  wim_pattern <- pattern_info$wim_pattern
  
  # Identify columns
  matrix_col <- str_subset(colnames(data), regex(matrix_pattern))
  wim_col <- str_subset(colnames(data), regex(wim_pattern))
  
  # Check if multiple or no columns are matched
  if (length(matrix_col) == 0 || length(wim_col) == 0) {
    message(sprintf("Required columns not found in the dataset for file: %s", basename(file_path)))
    return(sprintf('File: %s\nRequired columns not found in the dataset.\n', basename(file_path)))
  } else if (length(matrix_col) > 1 || length(wim_col) > 1) {
    message(sprintf("Multiple columns matched the pattern in file: %s. Please refine the patterns.", basename(file_path)))
    return(sprintf('File: %s\nMultiple columns matched the pattern.\n', basename(file_path)))
  }
  
  # Calculate the Matrix/WIM Area ratio
  ratio_values <- (as.numeric(data[[matrix_col]]) / as.numeric(data[[wim_col]])) * 100
  
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
  output_file_filtered <- file.path(output_dir, paste0('filtered_', basename(file_path)))
  write.table(data_filtered, file = output_file_filtered, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save the discarded data
  output_file_discarded <- file.path(output_dir, paste0('discarded_', basename(file_path)))
  write.table(data_discarded, file = output_file_discarded, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Prepare log information
  log_info <- sprintf('File: %s\nData Type: %s\nThreshold: %.2f%%\n%d rows removed (%.2f%%) because Matrix/WIM ratio < %.2f%%.\nRemoved rows: %s\n', 
                      basename(file_path), data_type, threshold, removed_rows, percent_removed, threshold, toString(names_removed))
  
  # Save log information to a file in the output directory
  log_file <- file.path(output_dir, "Data_Quality_Log.txt")
  writeLines(log_info, con = log_file)
  
  message(sprintf("Processing complete. Results saved in folder: %s", output_dir))
  
  return(invisible(NULL))
}

# Main code
# Ask for user input
file_path <- readline(prompt = "Enter the path to the data file: ")
data_type <- readline(prompt = "Enter the data type (GS, pFAK, PALLD, pSMAD): ")
threshold_input <- readline(prompt = "Enter the threshold value: ")

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
