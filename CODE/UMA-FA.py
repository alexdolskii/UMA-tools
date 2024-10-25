import imagej
import os
from pathlib import Path
from scyjava import jimport
import pandas as pd
from datetime import datetime

# Initialize ImageJ
print("Initializing ImageJ...")
ij = imagej.init(r"C:\Users\dolsk\Desktop\Orientation assay script\Fiji.app", mode='interactive')
print("ImageJ initialization completed.")

# Import Java classes
IJ = jimport('ij.IJ')
ImagePlus = jimport('ij.ImagePlus')
ZProjector = jimport('ij.plugin.ZProjector')
WindowManager = jimport('ij.WindowManager')
ResultsTable = jimport('ij.measure.ResultsTable')

# Desired image size
desired_width = 500
desired_height = 500

# Step 1: Prompt the user for the input folder
input_folder = input("Enter the full path to the folder containing the primary image files: ")

# Verify if the folder exists
if not os.path.exists(input_folder):
    print(f"The folder '{input_folder}' does not exist. Please try again.")
    exit()

# Create a new processed_images folder with date and time
now = datetime.now()
timestamp = now.strftime("%Y%m%d_%H%M%S")
processed_folder_name = f"processed_images_{timestamp}"
processed_folder = os.path.join(input_folder, processed_folder_name)
Path(processed_folder).mkdir(parents=True, exist_ok=True)
print(f"\nProcessed images will be saved in: {processed_folder}")

# Create output folders for Excel files and images
excel_folder = os.path.join(processed_folder, "Excel")
images_folder = os.path.join(processed_folder, "Images")
Path(excel_folder).mkdir(parents=True, exist_ok=True)
Path(images_folder).mkdir(parents=True, exist_ok=True)
print(f"Excel and Images folders created in {processed_folder}")

# Part 1: Creating 2D Projections
print("\nStarting Part 1: Creating 2D Projections...")
z_stacks_info = {}  # Dictionary to store the number of Z-stacks for each processed image

for filename in os.listdir(input_folder):
    if filename.lower().endswith(('.tif', '.tiff', '.nd2')):
        file_path = os.path.join(input_folder, filename)
        print(f"\nProcessing file: {file_path}")

        # Try to open the image in ImageJ
        try:
            # Use Bio-Formats to open the image
            img = ij.io().open(file_path)
            if img is None:
                print(f"Failed to open image: {file_path}")
                continue
            else:
                print(f"Image '{filename}' opened successfully.")
        except Exception as e:
            print(f"Error opening image '{filename}': {e}")
            continue

        # Convert ImgPlus to ImagePlus
        imp = ij.convert().convert(img, ImagePlus)

        # Get image dimensions
        dimensions = imp.getDimensions()
        width, height, channels, slices, frames = dimensions
        stack_size = imp.getStackSize()
        print(f"Image '{filename}' dimensions: width={width}, height={height}, channels={channels}, slices={slices}, frames={frames}, stack size={stack_size}")

        # Determine the initial number of Z-stacks (slices or channels)
        initial_z_stacks = slices if slices > 1 else channels
        z_stack_type = 'slices' if slices > 1 else 'channels'
        print(f"Initial number of Z-stacks ({z_stack_type}) in '{filename}': {initial_z_stacks}")

        # Check the file extension to determine processing steps
        file_ext = os.path.splitext(filename)[1].lower()

        if file_ext == '.nd2':
            # Process .nd2 files: extract the red (Cy5) channel
            if channels > 1:
                print(f"Image '{filename}' has {channels} channels. Extracting the red (Cy5) channel.")
                # Adjust the channel index based on your data
                red_channel_index = 1  # Index starts at 0; adjust if necessary
                imp.setC(red_channel_index + 1)  # Channels in ImageJ are 1-based
                IJ.run(imp, "Make Substack...", f"channels={red_channel_index+1} slices=1-{slices} frames=1-{frames}")
                imp.close()
                imp = IJ.getImage()
                channels = imp.getNChannels()  # Update channel count after extraction
                print(f"Red channel extracted from image '{filename}'.")
            else:
                print(f"Image '{filename}' has a single channel.")
        else:
            # For .tiff files, check if Z-stack is stored in channels
            if channels > 1:
                print(f"Image '{filename}' has {channels} channels. Performing Max Intensity Projection across channels.")
                imp.show()
                IJ.run("Z Project...", "projection=[Max Intensity] all")
                imp.close()
                imp = IJ.getImage()
                channels = imp.getNChannels()
                print(f"Max Intensity Projection across channels completed for image '{filename}'.")
            else:
                print(f"No channel extraction or projection needed for '{filename}'.")

        # Update slices after potential channel extraction or projection
        slices = imp.getNSlices()
        print(f"Number of Z-stacks (slices) in '{filename}' after processing: {slices}")

        # Check if it's a Z-stack in slices dimension
        if slices > 1:
            print(f"Image '{filename}' is a Z-stack with {slices} slices.")
            # Perform Max Intensity Projection along slices
            zp = ZProjector(imp)
            zp.setMethod(ZProjector.MAX_METHOD)
            zp.doProjection()
            projected_imp = zp.getProjection()
            imp.close()  # Close the original image to free memory
            imp = projected_imp
            print(f"Max Intensity Projection along slices completed for image '{filename}'.")
        else:
            print(f"Image '{filename}' is already a 2D image.")

        # Resize the image
        print(f"Resizing image '{filename}' to {desired_width}x{desired_height} pixels.")
        IJ.run(imp, "Size...", f"width={desired_width} height={desired_height} constrain average interpolation=Bilinear")

        # Save the processed image
        output_filename = os.path.splitext(filename)[0] + '_processed.tif'
        output_path = os.path.join(processed_folder, output_filename)
        IJ.saveAs(imp, "Tiff", output_path)
        print(f"Processed image saved to '{output_path}'.")
        imp.close()

        # Store the number of Z-stacks for the processed image
        processed_base_name = os.path.splitext(output_filename)[0]
        z_stacks_info[processed_base_name] = {
            'original_filename': filename,
            'number_of_z_stacks': initial_z_stacks,
            'z_stack_type': z_stack_type
        }

print("\nPart 1 completed successfully.")

# Part 2: Applying OrientationJ Plugin
print("\nStarting Part 2: Applying OrientationJ Plugin...")
processed_files = [f for f in os.listdir(processed_folder) if f.lower().endswith('_processed.tif')]

if not processed_files:
    print("No processed images found. Ensure that Part 1 completed successfully.")
    exit()

for filename in processed_files:
    file_path = os.path.join(processed_folder, filename)
    print(f"\nProcessing file: {file_path}")

    # Try to open the image in ImageJ
    try:
        imp = IJ.openImage(file_path)
        if imp is None:
            print(f"Failed to open image: {file_path}")
            continue
        else:
            print(f"Image '{filename}' opened successfully.")
    except Exception as e:
        print(f"Error opening image '{filename}': {e}")
        continue

    try:
        # Ensure imp is the active image
        imp.show()

        # Perform OrientationJ Analysis
        print(f"Performing OrientationJ Analysis on '{filename}'.")
        IJ.run("OrientationJ Analysis", "tensor=3.0 gradient=4 color-survey=on hsb=on hue=Orientation sat=Coherency bri=Original-Image radian=on")
        IJ.wait(500)  # Wait for the plugin to finish

        # Get the analysis image
        analysis_title = "OJ-Color-survey-1"
        analysis_imp = WindowManager.getImage(analysis_title)
        if analysis_imp is not None:
            analysis_filename = os.path.splitext(filename)[0] + '_oj_analysis.tif'
            analysis_path = os.path.join(images_folder, analysis_filename)
            IJ.saveAs(analysis_imp, "Tiff", analysis_path)
            print(f"OrientationJ Analysis image saved to '{analysis_path}'.")
            analysis_imp.close()
        else:
            print(f"No image found for OrientationJ Analysis for '{filename}'.")

        # Close all images
        IJ.run("Close All")

        # Re-open the image
        imp = IJ.openImage(file_path)
        if imp is None:
            print(f"Failed to re-open image: {file_path}")
            continue
        imp.show()

        # Perform OrientationJ Distribution
        print(f"Performing OrientationJ Distribution on '{filename}'.")
        IJ.run("OrientationJ Distribution", "tensor=3.0 gradient=4 radian=on histogram=on table=on min-coherency=0.0 min-energy=0.0")
        IJ.wait(500)  # Wait for the plugin to finish

        # Save the Results table
        try:
            excel_filename = os.path.splitext(filename)[0] + '_oj_distribution.csv'
            excel_path = os.path.join(excel_folder, excel_filename)
            # Save the Results table
            IJ.saveAs("Results", excel_path)
            print(f"OrientationJ Distribution results saved to '{excel_path}'.")
            # Clear the Results table
            IJ.run("Clear Results")
        except Exception as e:
            print(f"Error saving OrientationJ Distribution results for '{filename}': {e}")

        # Close all images and windows
        IJ.run("Close All")

    except Exception as e:
        print(f"Error during processing of '{filename}': {e}")
        # Close any open images or windows
        IJ.run("Close All")

print("\nPart 2 completed successfully.")

# Part 3: Processing the Excel Files
print("\nStarting Part 3: Processing Excel Files...")

# Verify if the Excel folder exists and has CSV files
if not os.path.exists(excel_folder):
    print(f"The folder '{excel_folder}' does not exist. Ensure that Part 2 completed successfully.")
    exit()

# Get list of CSV files in the Excel folder
file_list = [f for f in os.listdir(excel_folder) if f.endswith('.csv')]

if not file_list:
    print(f"No CSV files found in '{excel_folder}'. Ensure that Part 2 completed successfully.")
    exit()

# Create alignment_analysis_results folder inside the processed folder
analysis_folder = os.path.join(processed_folder, 'alignment_analysis_results')
if not os.path.exists(analysis_folder):
    os.makedirs(analysis_folder)

# Create a new Results folder with the current date and time
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
results_folder = os.path.join(analysis_folder, f'Results_{timestamp}')
os.makedirs(results_folder)

# Initialize a list to store summary data
summary_data = []

# Function to process individual file
def process_file(file_path):
    # Read the CSV file into a DataFrame
    read_file = pd.read_csv(file_path)

    # Rename the first and second columns
    read_file.rename(columns={read_file.columns[0]: "orientation_angle", read_file.columns[1]: "occurrence_value"}, inplace=True)

    # Find the maximum value in the "occurrence_value" column
    max_occurrence_value = read_file['occurrence_value'].max()

    # Find the corresponding value in the first column "orientation_angle"
    angle_of_max_occurrence_value = read_file['orientation_angle'][read_file['occurrence_value'].idxmax()]

    # Normalize angles to angle of maximum occurrence value (MOV)
    read_file['angles_normalized_to_angle_of_MOV'] = read_file['orientation_angle'] - angle_of_max_occurrence_value

    # Apply the transformation based on conditions for 'angles_normalized_to_angle_of_MOV'
    read_file['corrected_angles'] = read_file['angles_normalized_to_angle_of_MOV'].apply(
        lambda x: x + 180 if x < -90 else (x - 180 if x > 90 else x))

    # Rank the 'corrected_angles' and store the result
    read_file['rank_of_angle_occurrence_value'] = read_file['corrected_angles'].rank(method='min')

    # Calculate the sum of occurrence values
    sum_of_occurrence_values = read_file['occurrence_value'].sum()

    # Calculate percentage occurrence values
    read_file['percent_occurrence_value_to_sum_of_occurrence_value'] = (read_file['occurrence_value'] / sum_of_occurrence_values) * 100

    # Filter rows where 'corrected_angles' are between -15 and 15
    filtered_data = read_file[(read_file['corrected_angles'] >= -15) & (read_file['corrected_angles'] <= 15)]

    # Sum the 'percent_occurrence_value_to_sum_of_occurrence_value' for filtered rows
    percentage_of_fibers_aligned_within_15_degree = filtered_data['percent_occurrence_value_to_sum_of_occurrence_value'].sum()

    # Determine orientation mode based on percentage
    orientation_mode = 'disorganized' if percentage_of_fibers_aligned_within_15_degree < 55 else 'aligned'

    return read_file, percentage_of_fibers_aligned_within_15_degree, orientation_mode

# Process each CSV file
for file_name in file_list:
    file_path = os.path.join(excel_folder, file_name)
    print(f"\nProcessing Excel file: {file_name}")
    try:
        read_file, percentage_of_fibers_aligned_within_15_degree, orientation_mode = process_file(file_path)

        # Sort the DataFrame by 'rank_of_angle_occurrence_value'
        read_file_sorted = read_file.sort_values(by='rank_of_angle_occurrence_value')

        # Generate output file name and save the sorted DataFrame as an Excel file
        output_file_name = f'{os.path.splitext(file_name)[0]}_processed_{timestamp}.xlsx'
        output_file_path = os.path.join(results_folder, output_file_name)
        read_file_sorted.to_excel(output_file_path, index=False)
        print(f"Processed data saved to: {output_file_path}")

        # Get the processed base name corresponding to the CSV file
        processed_base_name = os.path.splitext(file_name.replace('_oj_distribution.csv', ''))[0]

        # Retrieve the number of Z-stacks for the processed image
        z_stacks_info_entry = z_stacks_info.get(processed_base_name, None)
        if z_stacks_info_entry is not None:
            number_of_z_stacks = z_stacks_info_entry['number_of_z_stacks']
            z_stack_type = z_stacks_info_entry['z_stack_type']
        else:
            number_of_z_stacks = 'N/A'
            z_stack_type = 'N/A'

        # Append summary data
        summary_data.append({
            'File_Name': file_name,
            'Number_of_Z_Stacks': number_of_z_stacks,
            'Z_Stack_Type': z_stack_type,
            'Percentage_Fibers_Aligned_Within_15_Degree': percentage_of_fibers_aligned_within_15_degree,
            'Orientation_Mode': orientation_mode
        })
    except Exception as e:
        print(f"Error processing file '{file_name}': {e}")

# Save the summary data as an Excel file
summary_df = pd.DataFrame(summary_data)
summary_file_path = os.path.join(results_folder, f'Summary_{timestamp}.xlsx')
summary_df.to_excel(summary_file_path, index=False)
print(f"\nSummary data saved to: {summary_file_path}")

print(f"\nProcessing completed. All results are saved in folder: {processed_folder}")
