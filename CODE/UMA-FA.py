import imagej
import os
from pathlib import Path
from scyjava import jimport
import pandas as pd
from datetime import datetime
import sys
import time

def main():
    # Initialize ImageJ in interactive mode
    print("Initializing ImageJ...")
    ij = imagej.init(r"C:\Users\dolsk\Desktop\Orientation_assay_script\Fiji.app", mode='interactive')
    print("ImageJ initialization completed.")

    # Import Java classes
    IJ = jimport('ij.IJ')
    ImagePlus = jimport('ij.ImagePlus')
    ZProjector = jimport('ij.plugin.ZProjector')
    WindowManager = jimport('ij.WindowManager')
    ResultsTable = jimport('ij.measure.ResultsTable')
    System = jimport('java.lang.System')  # For System.exit(0)

    # Desired image size
    desired_width = 500
    desired_height = 500

    # Step 1: Request txt file with folder paths
    txt_file_path = input("Enter the full path to the .txt file containing folder paths for processing: ").strip()

    # Check if the file exists
    if not os.path.isfile(txt_file_path):
        print(f"The file '{txt_file_path}' does not exist. Please try again.")
        sys.exit(1)

    # Read folder paths from the txt file
    with open(txt_file_path, 'r', encoding='utf-8') as f:
        folder_paths = [line.strip() for line in f if line.strip()]

    # Check if there are any folder paths
    if not folder_paths:
        print("The file does not contain any folder paths. Please check the file content.")
        sys.exit(1)

    # Print how many folder paths were found
    print(f"\nFound {len(folder_paths)} folder paths.")

    # For each folder, print the number of files and their formats
    folder_info = []
    for idx, folder in enumerate(folder_paths, start=1):
        print(f"\nFolder {idx}: {folder}")
        if not os.path.exists(folder):
            print(f"Folder '{folder}' does not exist.")
            continue
        files_in_folder = os.listdir(folder)
        file_types = {}
        for filename in files_in_folder:
            ext = os.path.splitext(filename)[1].lower()
            if ext in file_types:
                file_types[ext] += 1
            else:
                file_types[ext] = 1
        total_files = sum(file_types.values())
        print(f"  Total files: {total_files}")
        print("  File formats:")
        for ext, count in file_types.items():
            print(f"    {ext}: {count}")
        folder_info.append({
            'folder': folder,
            'total_files': total_files,
            'file_types': file_types,
        })

    # Ask the user for the angle to use in analysis (default is 15)
    angle_input = input("\nEnter the angle to use for analysis (default is 15): ").strip()
    if angle_input == '':
        angle_value = 15
    else:
        try:
            angle_value = float(angle_input)
        except ValueError:
            print("Invalid input. Using default angle value of 15.")
            angle_value = 15

    # Convert angle_value to a string suitable for folder naming
    angle_str = f"{angle_value}".replace('.', '_')

    # Ask whether all folders have the same fibronectin channel
    same_channel = input("\nDo all folders have the same fibronectin channel? (yes/no): ").strip().lower()
    fibronectin_channel_indices = {}  # Dictionary to store channel index per folder

    if same_channel in ('yes', 'y'):
        # Ask for the fibronectin channel index once
        channel_input = input("Enter the fibronectin channel index (starting from 0): ").strip()
        try:
            fibronectin_channel_index = int(channel_input)
        except ValueError:
            print("Invalid input. Using default channel index of 0.")
            fibronectin_channel_index = 0
        # Set the same channel index for all folders
        for folder in folder_paths:
            fibronectin_channel_indices[folder] = fibronectin_channel_index
    else:
        # Ask for the fibronectin channel index for each folder
        for folder in folder_paths:
            print(f"\nFolder: {folder}")
            channel_input = input("Enter the fibronectin channel index (starting from 0): ").strip()
            try:
                fibronectin_channel_index = int(channel_input)
            except ValueError:
                print("Invalid input. Using default channel index of 0.")
                fibronectin_channel_index = 0
            fibronectin_channel_indices[folder] = fibronectin_channel_index

    # Ask whether to start processing
    start_processing = input("\nDo you want to start processing the files? (yes/no): ").strip().lower()
    if start_processing not in ('yes', 'y'):
        print("File processing canceled by the user.")
        sys.exit(0)

    # Get current date and time for folder naming
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S")

    # Part 1: Creating 2D Projections for each folder
    print("\nStarting Part 1: Creating 2D Projections...")

    z_stacks_info = {}  # Dictionary to store the number of Z-stacks for each processed image

    for folder_idx, input_folder in enumerate(folder_paths, start=1):
        print(f"\nProcessing folder {folder_idx}: {input_folder}")
        if not os.path.exists(input_folder):
            print(f"The folder '{input_folder}' does not exist. Skipping this folder.")
            continue

        # Retrieve the fibronectin channel index for this folder
        fibronectin_channel_index = fibronectin_channel_indices.get(input_folder, 0)

        # Create a new Alignment_assay_results folder with date, time, and angle inside each input folder
        results_folder_name = f"Alignment_assay_results_angle_{angle_str}_{timestamp}"
        results_folder = os.path.join(input_folder, results_folder_name)
        Path(results_folder).mkdir(parents=True, exist_ok=True)
        print(f"Results will be saved in: {results_folder}")

        # Create output folders for Excel files and images
        excel_folder = os.path.join(results_folder, "Excel")
        images_folder = os.path.join(results_folder, "Images")
        Path(excel_folder).mkdir(parents=True, exist_ok=True)
        Path(images_folder).mkdir(parents=True, exist_ok=True)
        print(f"Excel and Images folders created in {results_folder}")

        # Process files in the current folder
        z_stacks_info_folder = {}
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
                    # Process .nd2 files: extract the fibronectin channel
                    if channels > 1:
                        print(f"Image '{filename}' has {channels} channels. Extracting the fibronectin channel (index {fibronectin_channel_index}).")
                        # Adjust the channel index based on user input
                        imp.setC(fibronectin_channel_index + 1)  # Channels in ImageJ are 1-based
                        IJ.run(imp, "Make Substack...", f"channels={fibronectin_channel_index+1} slices=1-{slices} frames=1-{frames}")
                        imp.close()
                        imp = IJ.getImage()
                        channels = imp.getNChannels()  # Update channel count after extraction
                        print(f"Fibronectin channel extracted from image '{filename}'.")
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
                output_path = os.path.join(results_folder, output_filename)
                IJ.saveAs(imp, "Tiff", output_path)
                print(f"Processed image saved to '{output_path}'.")
                imp.close()

                # Store the number of Z-stacks for the processed image
                processed_base_name = os.path.splitext(output_filename)[0]
                z_stacks_info_folder[processed_base_name] = {
                    'original_filename': filename,
                    'number_of_z_stacks': initial_z_stacks,
                    'z_stack_type': z_stack_type
                }
        # Update the main z_stacks_info dictionary
        z_stacks_info.update(z_stacks_info_folder)

        print(f"\nPart 1 completed successfully for folder {input_folder}.")

        # Part 2: Applying OrientationJ Plugin
        print("\nStarting Part 2: Applying OrientationJ Plugin...")
        processed_files = [f for f in os.listdir(results_folder) if f.lower().endswith('_processed.tif')]

        if not processed_files:
            print("No processed images found. Ensure that Part 1 completed successfully.")
            continue

        for filename in processed_files:
            file_path = os.path.join(results_folder, filename)
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

        print(f"\nPart 2 completed successfully for folder {input_folder}.")

        # Part 3: Processing the Excel Files
        print("\nStarting Part 3: Processing Excel Files...")

        # Verify if the Excel folder exists and has CSV files
        if not os.path.exists(excel_folder):
            print(f"The folder '{excel_folder}' does not exist. Ensure that Part 2 completed successfully.")
            continue

        # Get list of CSV files in the Excel folder
        file_list = [f for f in os.listdir(excel_folder) if f.endswith('.csv')]

        if not file_list:
            print(f"No CSV files found in '{excel_folder}'. Ensure that Part 2 completed successfully.")
            continue

        # Create a new Analysis folder inside the results folder
        analysis_folder = os.path.join(results_folder, 'Analysis')
        if not os.path.exists(analysis_folder):
            os.makedirs(analysis_folder)

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

            # Filter rows where 'corrected_angles' are between -angle_value and angle_value
            filtered_data = read_file[(read_file['corrected_angles'] >= -angle_value) & (read_file['corrected_angles'] <= angle_value)]

            # Sum the 'percent_occurrence_value_to_sum_of_occurrence_value' for filtered rows
            percentage_of_fibers_aligned_within_angle = filtered_data['percent_occurrence_value_to_sum_of_occurrence_value'].sum()

            # Determine orientation mode based on percentage
            orientation_mode = 'disorganized' if percentage_of_fibers_aligned_within_angle < 55 else 'aligned'

            return read_file, percentage_of_fibers_aligned_within_angle, orientation_mode

        # Process each CSV file
        for file_name in file_list:
            file_path = os.path.join(excel_folder, file_name)
            print(f"\nProcessing Excel file: {file_name}")
            try:
                read_file, percentage_of_fibers_aligned_within_angle, orientation_mode = process_file(file_path)

                # Sort the DataFrame by 'rank_of_angle_occurrence_value'
                read_file_sorted = read_file.sort_values(by='rank_of_angle_occurrence_value')

                # Generate output file name and save the sorted DataFrame as an Excel file
                output_file_name = f'{os.path.splitext(file_name)[0]}_processed_{timestamp}.xlsx'
                output_file_path = os.path.join(analysis_folder, output_file_name)
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
                    f'Percentage_Fibers_Aligned_Within_{angle_value}_Degree': percentage_of_fibers_aligned_within_angle,
                    'Orientation_Mode': orientation_mode
                })
            except Exception as e:
                print(f"Error processing file '{file_name}': {e}")

        # Save the summary data as an Excel file
        summary_df = pd.DataFrame(summary_data)
        summary_file_path = os.path.join(analysis_folder, f'Alignment_Summary_{timestamp}.xlsx')
        summary_df.to_excel(summary_file_path, index=False)
        print(f"\nSummary data saved to: {summary_file_path}")

        print(f"\nProcessing completed for folder {input_folder}. All results are saved in folder: {results_folder}")

    print("\nAll folders have been processed.")

    # Dispose of the ImageJ context
    print("Disposing of ImageJ context...")
    ij.dispose()
    print("Script execution completed.")

if __name__ == "__main__":
    main()
