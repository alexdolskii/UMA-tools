# Data Analysis Script for Markers: pSMAD, PALLD, pFAK, GS
# Author: [Your Name]
# Date: [Date]
#
# This script processes and analyzes biological data for four markers: pSMAD, PALLD, pFAK, and GS.
# It prompts the user for the location of the raw data file and the Excel file containing conditions and groups.
# Depending on the specified marker, it performs data processing, plotting, and statistical analysis.
# All outputs are saved in a separate folder named 'Marker_intensity_level_analysis' in the script's directory.
# Each run creates a new subfolder with the date, time, input file name, and marker.
# A log file is created in the output folder containing the names and paths of the input files.
#
# Requirements:
# - R version >= 3.6.0
# - Packages: ggplot2, ggpubr, readxl, dplyr, stringr
#
# Usage:
# - Run the script in R or RStudio.
# - Before calling the 'analyze_marker' function, please note that you need to provide a marker as a parameter.
#   For example: analyze_marker("pSMAD")
# - If no marker is specified, the script will prompt you to enter one.
# - The script will prompt you to enter the file paths for the raw data and the Excel file.
# - You will be prompted to select the desired image quality.
#
# Example:
# analyze_marker("pSMAD")
#
# License: MIT License

# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(stringr)  # Added for string manipulation

# Function to process data based on the marker
process_marker_data <- function(file_path, marker) {
  data <- read.table(file_path, header = TRUE, sep = "\t", dec = ".")
  
  if (marker == "pSMAD") {
    part_data <- data %>%
      select(
        Well.Name,
        psmad_nuc_count = matches("Cell\\.\\.psmad_nuc\\.obj\\.msk_Features\\.Count_Sum"),
        psmad_posnuc_intensity = matches("Cell\\.\\.psmad_posnuc_objmsk_Integrated\\.Intensity_Sum"),
        psmad_matrix_intensity = matches("Cell\\.\\.psmad_matrix\\.objmsk_Integrated\\.Intensity_Sum")
      ) %>%
      mutate(
        intensity_to_nuclei = psmad_posnuc_intensity / psmad_nuc_count,
        intensity_to_matrix = psmad_posnuc_intensity / psmad_matrix_intensity
      )
  } else if (marker == "PALLD") {
    part_data <- data %>%
      select(
        Well.Name,
        palld_nuc_count = matches("^Cell\\.\\.iso3_nuc\\.obj\\.msk__Features\\.Count_Sum.*"),
        palld_intensity = matches("^Cell\\.\\.iso3_obj\\.msk_Integrated\\.Intensity_Sum.*"),
        palld_matrix_intensity = matches("^Cell\\.\\.iso3_\\.Matrix\\.objmsk_Integrated\\.Intensity_Sum.*")
      )
    
    # Check for missing columns
    required_columns <- c("palld_nuc_count", "palld_intensity", "palld_matrix_intensity")
    missing_columns <- setdiff(required_columns, names(part_data))
    if (length(missing_columns) > 0) {
      stop(paste("The following required columns are missing for PALLD:", paste(missing_columns, collapse = ", ")))
    }
    
    part_data <- part_data %>%
      mutate(
        intensity_to_nuclei = palld_intensity / palld_nuc_count,
        intensity_to_matrix = palld_intensity / palld_matrix_intensity
      )
  } else if (marker == "pFAK") {
    part_data <- data %>%
      select(
        Well.Name,
        pfak_intensity = matches("Cell\\.\\.Tritc_pFAK_obj\\.msk_Integrated\\.Intensity_Sum"),
        pfak_matrix_intensity = matches("Cell\\.\\.pFAK\\.Matrix\\.objmsk_Integrated\\.Intensity_Sum"),
        pfak_nuc_count = matches("Cell\\.\\.pFAK_nuc\\.obj\\.msk__Features\\.Count_Sum"),
        pfak_nuc_area = matches("Cell\\.\\.pFAK_nuc\\.obj\\.msk__Total\\.Area_Sum")
      ) %>%
      mutate(
        intensity_to_nuclei = pfak_intensity / pfak_nuc_area,
        intensity_to_matrix = pfak_intensity / pfak_matrix_intensity
      )
  } else if (marker == "GS") {
    part_data <- data %>%
      select(
        Well.Name,
        gs_nuc_count = matches("Cell\\.\\.GS_nuc\\.objmsk__Features\\.Count_Sum"),
        gs_intensity = matches("Cell\\.\\.GS_objmsk_Integrated\\.Intensity_Sum"),
        gs_matrix_intensity = matches("Cell\\.\\.GS_\\.Matrix\\.objmsk_Integrated\\.Intensity_Sum")
      ) %>%
      mutate(
        intensity_to_nuclei = gs_intensity / gs_nuc_count,
        intensity_to_matrix = gs_intensity / gs_matrix_intensity
      )
  } else {
    stop("Invalid marker specified. Please choose from 'pSMAD', 'PALLD', 'pFAK', or 'GS'.")
  }
  
  return(part_data)
}

# Function to define conditions and groups from an Excel file
define_conditions_groups_from_excel <- function(data, excel_file) {
  conditions_groups <- read_excel(excel_file)
  
  data <- data %>%
    left_join(conditions_groups, by = "Well.Name") %>%
    mutate(
      Conditions = as.factor(Conditions),
      Group = as.factor(Group)
    )
  
  if (any(is.na(data$Conditions)) || any(is.na(data$Group))) {
    stop("Conditions or Group columns contain NA values after merging. Please check the input Excel file.")
  }
  
  return(data)
}

# General plotting function
plot_data <- function(data, y_var, title, y_label) {
  ggplot(data, aes(x = Conditions, y = .data[[y_var]], fill = Group)) +
    geom_violin(trim = FALSE) +
    geom_jitter(position = position_jitter(0.05), size = 0.8) +
    theme_minimal() +
    guides(fill = "none") +  # Updated line to fix warning
    ggtitle(title) +
    xlab("Conditions") +
    ylab(y_label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(~ Group, scales = "free_x", space = "free_x") +
    scale_fill_brewer(palette = "Accent")
}

# Function to plot data based on the marker
plot_marker_data <- function(data, marker, image_resolution, output_dir) {
  # Define plotting variables based on the marker
  if (marker == "pSMAD") {
    intensity_col <- "psmad_posnuc_intensity"
    nuc_count_col <- "psmad_nuc_count"
    matrix_intensity_col <- "psmad_matrix_intensity"
  } else if (marker == "PALLD") {
    intensity_col <- "palld_intensity"
    nuc_count_col <- "palld_nuc_count"
    matrix_intensity_col <- "palld_matrix_intensity"
  } else if (marker == "pFAK") {
    intensity_col <- "pfak_intensity"
    nuc_count_col <- "pfak_nuc_count"
    matrix_intensity_col <- "pfak_matrix_intensity"
  } else if (marker == "GS") {
    intensity_col <- "gs_intensity"
    nuc_count_col <- "gs_nuc_count"
    matrix_intensity_col <- "gs_matrix_intensity"
  }
  
  # Remove non-finite values from data
  data <- data %>%
    filter(
      is.finite(.data[[intensity_col]]) &
        is.finite(.data[[nuc_count_col]]) &
        is.finite(.data[[matrix_intensity_col]]) &
        is.finite(intensity_to_nuclei) &
        is.finite(intensity_to_matrix)
    )
  
  # Create plots
  plots <- list(
    intensity_to_nuclei = plot_data(
      data,
      y_var = "intensity_to_nuclei",
      title = paste(marker, "Intensity to Nuclei Counts"),
      y_label = "Intensity"
    ),
    intensity_to_matrix = plot_data(
      data,
      y_var = "intensity_to_matrix",
      title = paste(marker, "Intensity to Matrix"),
      y_label = "Intensity Ratio"
    ),
    nuclei_counts = plot_data(
      data,
      y_var = nuc_count_col,
      title = "Nuclei Counts in Wells",
      y_label = "Count"
    ),
    matrix_intensity = plot_data(
      data,
      y_var = matrix_intensity_col,
      title = "Matrix Intensity in Wells",
      y_label = "Intensity"
    )
  )
  
  # Arrange plots
  plot_grid <- ggarrange(
    plots$intensity_to_nuclei,
    plots$intensity_to_matrix,
    plots$nuclei_counts,
    plots$matrix_intensity,
    nrow = 2, ncol = 2
  )
  
  # Save the combined plot as a JPEG file in the output directory
  jpeg_filename <- file.path(output_dir, paste0(marker, "_combined_plots.jpeg"))
  jpeg(jpeg_filename, width = image_resolution$width, height = image_resolution$height, units = "px", res = image_resolution$res)
  print(plot_grid)
  dev.off()
  
  cat("Combined plot saved as", jpeg_filename, "\n")
  
  # Save individual plots
  plot_names <- c("intensity_to_nuclei", "intensity_to_matrix", "nuclei_counts", "matrix_intensity")
  for (i in seq_along(plot_names)) {
    plot_name <- plot_names[i]
    plot_obj <- plots[[plot_name]]
    individual_jpeg_filename <- file.path(output_dir, paste0(marker, "_", plot_name, ".jpeg"))
    jpeg(individual_jpeg_filename, width = image_resolution$width / 2, height = image_resolution$height / 2, units = "px", res = image_resolution$res)
    print(plot_obj)
    dev.off()
    cat("Individual plot saved as", individual_jpeg_filename, "\n")
  }
}

# Function to perform pairwise comparisons and save results to files
perform_pairwise_tests <- function(data, value_columns, marker, output_dir) {
  results <- list()
  for (value_column in value_columns) {
    test_results <- list()
    for (group in unique(data$Group)) {
      group_data <- subset(data, Group == group)
      if (nlevels(group_data$Conditions) > 1) {
        test_result <- pairwise.t.test(
          group_data[[value_column]],
          group_data$Conditions,
          p.adjust.method = "bonferroni"
        )
        test_results[[as.character(group)]] <- test_result
      }
    }
    results[[value_column]] <- test_results
  }
  
  # Save results to files in the output directory
  for (test_name in names(results)) {
    output_file <- file.path(output_dir, paste0(marker, "_", test_name, "_pairwise_results.txt"))
    sink(output_file)
    cat("Results of pairwise comparisons for", test_name, ":\n")
    test_results <- results[[test_name]]
    
    for (group in names(test_results)) {
      cat("\nGroup:", group, "\n")
      test_result <- test_results[[group]]
      p_values <- test_result$p.value
      significant_pairs <- which(p_values < 0.05, arr.ind = TRUE)
      
      if (length(significant_pairs) > 0) {
        for (i in seq_len(nrow(significant_pairs))) {
          row <- significant_pairs[i, "row"]
          col <- significant_pairs[i, "col"]
          condition1 <- rownames(p_values)[row]
          condition2 <- colnames(p_values)[col]
          p_value <- p_values[row, col]
          cat("Comparison:", condition1, "vs", condition2, "- p-value:", p_value, "\n")
        }
      } else {
        cat("No significant pairwise comparisons found.\n")
      }
    }
    sink()
    cat("Statistical results saved to", output_file, "\n")
  }
  
  return(results)
}

# Function to create a log file with input file names and paths
create_log_file <- function(data_file, excel_file, output_dir) {
  log_file <- file.path(output_dir, "analysis_log.txt")
  sink(log_file)
  cat("Analysis Log\n")
  cat("Date and Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Data File:", normalizePath(data_file), "\n")
  cat("Excel File:", normalizePath(excel_file), "\n")
  sink()
}

# Main function to analyze a marker
analyze_marker <- function(marker = NULL) {
  # Informing the user before using the analyze_marker function
  cat("To use this script, please call the 'analyze_marker' function with one of the following markers as a parameter: 'pSMAD', 'PALLD', 'pFAK', or 'GS'.")
  cat("\nExample: analyze_marker('pSMAD')\n")
  
  # If marker is missing or empty, prompt the user
  if (is.null(marker) || marker == "") {
    marker <- readline(prompt = "Please enter a marker (pSMAD, PALLD, pFAK, GS): ")
  }
  # Ensure marker is valid
  while (!marker %in% c("pSMAD", "PALLD", "pFAK", "GS")) {
    marker <- readline(prompt = "Invalid marker. Please enter one of 'pSMAD', 'PALLD', 'pFAK', or 'GS': ")
  }
  
  # Prompt user for file paths
  data_file <- readline(prompt = "Enter the path to the raw data file: ")
  excel_file <- readline(prompt = "Enter the path to the Excel file with conditions and groups: ")
  
  # Process file paths
  data_file <- gsub('\"', '', data_file)
  excel_file <- gsub('\"', '', excel_file)
  data_file <- gsub('\\\\', '/', data_file)
  excel_file <- gsub('\\\\', '/', excel_file)
  
  # Check if files exist
  if (!file.exists(data_file)) {
    stop("The specified data file does not exist.")
  }
  if (!file.exists(excel_file)) {
    stop("The specified Excel file does not exist.")
  }
  
  # Create main output directory 'Marker_intensity_level_analysis' in the script directory
  # Determine the script directory
  script_dir <- tryCatch({
    this_file <- dirname(normalizePath(sys.frames()[[1]]$ofile))
  }, error = function(e) {
    getwd()
  })
  main_output_dir <- file.path(script_dir, "Marker_intensity_level_analysis")
  dir.create(main_output_dir, showWarnings = FALSE)
  
  # Create output directory with date, time, data file name, and marker
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  data_file_name <- tools::file_path_sans_ext(basename(data_file))
  output_dir_name <- paste0("Analysis_", data_file_name, "_", marker, "_", timestamp)
  output_dir <- file.path(main_output_dir, output_dir_name)
  dir.create(output_dir, showWarnings = FALSE)
  
  # Create log file
  create_log_file(data_file, excel_file, output_dir)
  
  # Prompt for image quality
  cat("Select image quality for output plots:")
  cat("\n1: High (width=4000, height=3000)")
  cat("\n2: Medium (width=2000, height=1500)")
  cat("\n3: Low (width=1000, height=750)")
  quality_choice <- readline(prompt = "\nEnter the number corresponding to the desired quality (default is High): ")
  if (quality_choice == "") {
    quality_choice <- "1"
  }
  if (quality_choice == "1") {
    width <- 4000
    height <- 3000
    res <- 300
  } else if (quality_choice == "2") {
    width <- 2000
    height <- 1500
    res <- 150
  } else if (quality_choice == "3") {
    width <- 1000
    height <- 750
    res <- 72
  } else {
    # If invalid input, default to High
    width <- 4000
    height <- 3000
    res <- 300
    cat("Invalid choice. Defaulting to High quality.\n")
  }
  image_resolution <- list(width = width, height = height, res = res)
  
  # Process data
  part_data <- process_marker_data(data_file, marker)
  part_data <- define_conditions_groups_from_excel(part_data, excel_file)
  
  # Plot data and save as JPEG
  plot_marker_data(part_data, marker, image_resolution, output_dir)
  
  # Define variables to test
  variables_to_test <- switch(marker,
                              "pSMAD" = c("intensity_to_nuclei", "intensity_to_matrix", "psmad_nuc_count", "psmad_matrix_intensity"),
                              "PALLD" = c("intensity_to_nuclei", "intensity_to_matrix", "palld_nuc_count", "palld_matrix_intensity"),
                              "pFAK" = c("intensity_to_nuclei", "intensity_to_matrix", "pfak_nuc_count", "pfak_matrix_intensity"),
                              "GS" = c("intensity_to_nuclei", "intensity_to_matrix", "gs_nuc_count", "gs_matrix_intensity")
  )
  
  # Perform pairwise tests and save results
  pairwise_tests <- perform_pairwise_tests(part_data, variables_to_test, marker, output_dir)
  
  cat("All outputs have been saved to the folder:", output_dir, "\n")
}

# Informing the user before using the analyze_marker function
cat("To use this script, please call the 'analyze_marker' function with one of the following markers as a parameter: 'pSMAD', 'PALLD', 'pFAK', or 'GS'.")
cat("\nExample: analyze_marker('pSMAD')\n")

# Example usage:
# analyze_marker()  # If you don't specify a marker, you'll be prompted to enter one.
# analyze_marker("pSMAD")  # Replace "pSMAD" with "PALLD", "pFAK", or "GS" as needed
