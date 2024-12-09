# UMA-tools Setup Guide

## For Windows Users

### Step 1: Install WSL (Windows Subsystem for Linux)

1. Open Powerell with administrator privileges:
    Press the Windows key, type "Powerell," right-click on it, and select "Run as administrator".
2. Run the following command to install WSL:
    wsl --install
3. After the installation, restart your computer.
4. Once WSL is installed, open your WSL terminal and run:
    sudo apt update
    sudo apt upgrade -y

### Step 2: Clone a GitHub Repository on WSL
1. Check if Git is installed:
   git --version
   sudo apt update
   sudo apt install git
2. Choose the Directory for Cloning
   cd ~
   Or create another folder if preferred:
   mkdir projects
   cd projects
3. Clone the Specific Branch of the Repository:
   git clone -b UMA-tools-linux https://github.com/alexdolskii/UMA-tools.git
4. Verify the Result:
   cd UMA-tools
5. git branch (You will see a list of branches, and an asterisk (`*`) will appear next to `UMA-tools-linux`, confirming you are on the correct branch.)
6. git pull (This will fetch all new changes from the remote repository in the `UMA-tools-linux` branch.)
   

### Step 3: Install Miniconda
1. Download the Miniconda installer for Linux:
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.
2. Make the installer executable:
    chmod +x Miniconda3-latest-Linux-x86_64.
3. Run the installer script:
    ./Miniconda3-latest-Linux-x86_64.
4. Restart the ell to apply changes:
    source ~/.barc
5. Verify the installation:
    conda --version
6. Update Conda to the latest version:
    conda update conda

### Step 4: Create Environment for UMA Tools
1. Create the environment:
    conda env create -f uma_environment.yml -n uma_environment
    
    To check conda environments:
    conda info --envs
    To delete conda envioment:
    conda remove --name <environment_name> --all

2. Activate the environment:
    conda activate uma_environment
3. Make the main script executable:
    chmod +x code/alignment_analysis.py
    chmod +x code/thickness_analysis.py

### Step 5: Installing OrientationPy (v3+)
1. Clone the Repository. Download the  code for OrientationPy:
    git clone https://gitlab.com/epfl-center-for-imaging/orientationpy.git
2. Enter the cloned directory:
    cd orientationpy
3. Install OrientationPy and its dependencies:
    pip install .
4. Return to the UMA-tools directory:
    cd ..

### Step 6:  Running the UMA Tools Script
1. Modify `input_paths.json` to include your paths to the `.nd2` files.
2. Run the main analysis script:
   python ./code/alignment_analysis.py -i input_paths.json



# Dependencies and Tools Used

This program utilizes the following tools:

1. **OrientationPy**  
   [OrientationPy](https://epfl-center-for-imaging.gitlab.io/orientationpy/introduction.html) is a Python-based plugin used in this project for calculating the **alignment of fibronectin fibers**.

   - Repository: [OrientationPy](https://gitlab.com/epfl-center-for-imaging/orientationpy/)  
   - License: [The GNU General Public License, Version 3, 29 June 2007 (GPLv3)](https://gitlab.com/epfl-center-for-imaging/orientationpy/-/blob/main/LICENSE.md?ref_type=heads)

2. **Fiji**  
   [Fiji](https://fiji.sc/) (Fiji Is Just ImageJ) is an open-source distribution of ImageJ with a focus on image analysis. In this project, Fiji was used for **preprocessing image stacks**, including tasks such as contrast enhancement, filtering, and segmentation, to ensure the data is optimized for analysis in OrientationPy.

   - Repository: [Fiji](https://github.com/fiji/fiji)  
   - License: [GPL License](https://imagej.net/licensing/)

## Usage

- **OrientationPy**: Used for calculating the alignment of fibronectin fibers in the processed images. This helps quantify structural organization and provides insights into [specific biological context, e.g., "tumor microenvironment"].
- **Fiji**: Serves as the preprocessing tool to optimize the images.

