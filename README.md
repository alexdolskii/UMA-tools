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

### Step 2: Install Miniconda
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

### Step 3: Create Environment for UMA Tools
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

### Step 4: Installing OrientationPy (v3+)
1. Clone the Repository. Download the  code for OrientationPy:
    git clone https://gitlab.com/epfl-center-for-imaging/orientationpy.git
2. Enter the cloned directory:
    cd orientationpy
3. Install OrientationPy and its dependencies:
    pip install .
4. Return to the UMA-tools directory:
    cd ..

### Step 5:  Running the UMA Tools Script
1. Modify `input_paths.json` to include your paths to the `.nd2` files.
2. Run the main analysis script:
   python ./code/alignment_analysis.py -i input_paths.json

