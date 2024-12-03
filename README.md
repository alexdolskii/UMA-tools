# UMA-tools

For windows:

Install wsl
    Open PowerShell with administrator privileges. To do this, press the Windows key, type "PowerShell," right-click on PowerShell, and select "Run as administrator".
    - wsl --install
    After running the command, a computer restart will be required

Update wsl
    - sudo apt update
    - sudo apt upgrade -y

Download the Miniconda installer for Linux
    - wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
After downloading the file, make it executable
    - chmod +x Miniconda3-latest-Linux-x86_64.sh
Run the installer script
    - ./Miniconda3-latest-Linux-x86_64.sh

After the installation is complete, restart the shell to apply changes
    - source ~/.bashrc

Verify that conda is installed by running
    - conda --version

Update Conda to the latest version
    - conda update conda

Create enviroment for uma tools
    - conda env create -f uma_environment.yml -n uma_enviroment1
    - conda activate <name of enviroment>
    - chmod +x code/fibronectin_orientation_analysis.py

Modify input_paths.json with your paths to .nd2 files

Run code
    -./code/Step1.py -i input_paths.json