import imagej
import os
from pathlib import Path
import pandas as pd
from datetime import datetime
from scyjava import jimport
import sys
import time

def main():
    # Initialize ImageJ
    print("Initializing ImageJ...")
    try:
        ij = imagej.init(r"C:\Users\dolsk\Desktop\Orientation assay script\Fiji.app", mode='interactive')
        print("ImageJ initialization completed.")
    except Exception as e:
        print(f"Error initializing ImageJ: {e}")
        sys.exit(1)

    # Import necessary Java classes
    try:
        IJ = jimport('ij.IJ')
        ImagePlus = jimport('ij.ImagePlus')
        WindowManager = jimport('ij.WindowManager')
        ResultsTable = jimport('ij.measure.ResultsTable')
        Duplicator = jimport('ij.plugin.Duplicator')
    except Exception as e:
        print(f"Error importing Java classes: {e}")
        sys.exit(1)

    # Prompt user for file type
    print("\nSelect the file type to process:")
    print("1. .nd2")
    print("2. .tiff")
    file_type_choice = input("Enter 1 for .nd2 or 2 for .tiff: ").strip()

    if file_type_choice == '1':
        file_extension = '.nd2'
    elif file_type_choice == '2':
        file_extension = '.tiff'
    else:
        print("Invalid choice. Please run the script again and select 1 or 2.")
        sys.exit(1)

    # Prompt user for number of channels
    while True:
        try:
            num_channels = int(input("Enter the number of channels in the files: ").strip())
            if num_channels < 1:
                print("Number of channels must be at least 1.")
                continue
            break
        except ValueError:
            print("Please enter a valid integer for the number of channels.")

    # If multiple channels, prompt for fibronectin channel
    if num_channels > 1:
        while True:
            try:
                fibronectin_channel = int(input(f"Enter the channel number (1-{num_channels}) representing fibronectin: ").strip())
                if 1 <= fibronectin_channel <= num_channels:
                    break
                else:
                    print(f"Channel number must be between 1 and {num_channels}.")
            except ValueError:
                print("Please enter a valid integer for the channel number.")
    else:
        fibronectin_channel = 1  # Default to channel 1 if only one channel exists

    # Request path to folder with image files
    input_folder = input("\nEnter the full path to the folder with image files for analysis: ").strip()

    if not os.path.exists(input_folder):
        print(f"Folder '{input_folder}' does not exist.")
        sys.exit(1)

    # Find all files with the specified extension in the folder
    image_files = [f for f in os.listdir(input_folder) if f.lower().endswith(file_extension)]

    if not image_files:
        print(f"No {file_extension} files found in the specified folder.")
        sys.exit(1)

    # Create results folder with timestamp
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S")
    results_folder = os.path.join(input_folder, f"Results_{timestamp}")
    Path(results_folder).mkdir(parents=True, exist_ok=True)
    print(f"\nResults will be saved in: {results_folder}\n")

    # Initialize list for summary data
    summary_data = []

    # Process each file
    for filename in image_files:
        file_path = os.path.join(input_folder, filename)
        print(f"Processing file: {filename}")

        try:
            # Open image
            print("  Opening image...")
            img = ij.io().open(file_path)
            if img is None:
                print(f"  Failed to open image: {filename}")
                continue
            else:
                print(f"  Image '{filename}' successfully opened.")

            # Convert to ImagePlus
            imp = ij.convert().convert(img, ImagePlus)
            if imp is None:
                print(f"  Failed to convert image '{filename}' to ImagePlus.")
                continue

            # Extract fibronectin channel using Duplicator
            print("  Extracting fibronectin channel...")
            print(f"  Number of channels in image: {imp.getNChannels()}")
            if fibronectin_channel > imp.getNChannels() or fibronectin_channel < 1:
                print(f"  Specified fibronectin channel ({fibronectin_channel}) is invalid for image with {imp.getNChannels()} channels.")
                imp.close()
                continue

            imp_fibronectin = Duplicator().run(
                imp, 
                fibronectin_channel, 
                fibronectin_channel, 
                1, 
                imp.getNSlices(), 
                1, 
                imp.getNFrames()
            )
            imp.close()
            if imp_fibronectin is None:
                print(f"  Failed to extract fibronectin channel from '{filename}'.")
                continue
            imp_fibronectin.setTitle(f"{filename}_C{fibronectin_channel}")

            # Reslice the image
            print("  Performing Reslice...")
            IJ.run(imp_fibronectin, "Reslice [/]...", "output=0.500 start=Top flip rotate avoid")
            imp_fibronectin.close()
            time.sleep(2)  # Wait for ImageJ to complete the operation

            resliced_imp = IJ.getImage()
            if resliced_imp is None:
                print(f"  Failed to perform Reslice for '{filename}'.")
                continue
            resliced_imp.setTitle(f"Reslice_of_{filename}")

            # Perform Z projection
            print("  Performing Z projection...")
            IJ.run(resliced_imp, "Z Project...", "projection=[Max Intensity]")
            resliced_imp.close()
            time.sleep(2)

            projected_imp = IJ.getImage()
            if projected_imp is None:
                print(f"  Failed to perform Z projection for '{filename}'.")
                continue
            projected_imp.setTitle(f"MAX_Reslice_of_{filename}")

            # Apply Maximum filter
            print("  Applying Maximum filter...")
            IJ.run(projected_imp, "Maximum...", "radius=2")

            # Apply Gaussian Blur
            print("  Applying Gaussian Blur...")
            IJ.run(projected_imp, "Gaussian Blur...", "sigma=2 scaled")

            # Subtract background
            print("  Subtracting background...")
            IJ.run(projected_imp, "Subtract Background...", "rolling=50 sliding")

            # Apply threshold
            print("  Applying threshold...")
            IJ.setAutoThreshold(projected_imp, "Otsu dark no-reset")
            IJ.run("Options...", "black")  # Set BlackBackground=true
            IJ.run(projected_imp, "Convert to Mask", "")

            # Save the mask image
            mask_image_path = os.path.join(results_folder, f"Mask_{filename}.tif")
            IJ.saveAs(projected_imp, "Tiff", mask_image_path)
            print(f"  Mask image saved to '{mask_image_path}'.")

            # Launch Local Thickness plugin
            print("  Running Local Thickness...")
            # Capture image titles before running the plugin
            images_before = set(WindowManager.getImageTitles())
            
            # Run the Local Thickness plugin via macro
            macro_code = 'run("Local Thickness (masked, calibrated, silent)");'
            IJ.runMacro(macro_code)
            time.sleep(5)  # Wait for the plugin to finish

            # Capture image titles after running the plugin
            images_after = set(WindowManager.getImageTitles())
            new_images = images_after - images_before

            if not new_images:
                print(f"  Local Thickness plugin did not create a new image for '{filename}'.")
                continue

            # Get the new image title
            new_image_title = new_images.pop()
            local_thickness_imp = WindowManager.getImage(new_image_title)
            if local_thickness_imp is None:
                print(f"  Failed to get Local Thickness image for '{filename}'.")
                continue
            local_thickness_imp.setTitle(f"Local_Thickness_of_{filename}")

            # Save the Local Thickness image
            local_thickness_image_path = os.path.join(results_folder, f"Local_Thickness_{filename}.tif")
            IJ.saveAs(local_thickness_imp, "Tiff", local_thickness_image_path)
            print(f"  Local Thickness image saved to '{local_thickness_image_path}'.")

            # Reset ResultsTable before measurement to avoid data overlap
            rt = ResultsTable.getResultsTable()
            if rt is not None:
                rt.reset()

            # Measure thickness
            print("  Measuring thickness...")
            IJ.run(local_thickness_imp, "Measure", "")
            rt = ResultsTable.getResultsTable()
            if rt is None or rt.getCounter() == 0:
                area = None
                mean_thickness = None
                min_thickness = None
                max_thickness = None
                print(f"  No measurements obtained for '{filename}'.")
            else:
                # Extract all necessary measurements from the last row
                area = rt.getValue("Area", rt.getCounter() - 1)
                mean_thickness = rt.getValue("Mean", rt.getCounter() - 1)
                min_thickness = rt.getValue("Min", rt.getCounter() - 1)
                max_thickness = rt.getValue("Max", rt.getCounter() - 1)
                print(f"  Measurements for '{filename}': Area={area}, Mean={mean_thickness}, Min={min_thickness}, Max={max_thickness}")

            # Add measurements to summary_data
            summary_data.append({
                'File_Name': filename,
                'Area': area,
                'Mean': mean_thickness,
                'Min': min_thickness,
                'Max': max_thickness
            })

            # Close all images
            IJ.run("Close All")
            print("  Closed all images.\n")

        except Exception as e:
            print(f"  Error processing '{filename}': {e}")
            IJ.run("Close All")
            continue

    # After processing all files, save summary_data to CSV
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_file_path = os.path.join(results_folder, 'Thickness_Summary.csv')
        summary_df.to_csv(summary_file_path, index=False)
        print(f"Analysis complete. Summary data saved to '{summary_file_path}'.")
    else:
        print("No data to save in summary.")

if __name__ == "__main__":
    main()
