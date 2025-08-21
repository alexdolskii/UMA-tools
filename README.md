# UMA-tools


The program can be run on Linux or macOS machines. If you are using Windows, you will need to install WSL (Windows Subsystem for Linux).

For the code to function properly, specific package versions are required. Therefore, it is essential to work within a designated environment. As an example, commands for working with Miniconda are provided, but you are free to use any environment management tool that you find convenient. For macOS users, skip the WSL installation and work directly in the terminal.
## For Windows Users

The project is developed in the [Edna (Eti) Cukierman lab](https://www.foxchase.org/edna-cukierman).


## Installation 
To download and install *git* please visit [Git Download page](https://git-scm.com/downloads).

To download and install *conda* please visit [Miniforge github](https://github.com/conda-forge/miniforge)

Installing OrientationPy (v3+) The program requires a version higher than the one available via pip. Therefore, it needs to be installed via git (version 0.3.0 is used here). The official website is: https://epfl-center-for-imaging.gitlab.io/orientationpy/introduction.html 

1. Clone the Repository. Download the  code for OrientationPy:
    - git clone https://gitlab.com/epfl-center-for-imaging/orientationpy.git
2. Enter the cloned directory:
    - cd orientationpy
3. Install OrientationPy and its dependencies:
    - pip install .
4. Return to the UMA-tools directory:
    - cd ..

## Usage
1. Before running the program, you need to modify a `input_paths.json` file. This file should contain a list of folders with .nd2 images, and you can include as many folders as needed.

Additionally, before starting the program, make sure you know how many fluorescence channels you have (e.g., DAPI, Cy5) and their order in the file. You can check this by opening the image using the standard method in the GPU application FiJi (https://imagej.net/software/fiji/downloads).


2. Run the main analysis script:
   - python ./code/alignment_analysis.py -i input_paths.json
   - python ./code/thickness_analysis.py -i input_paths.json



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


## References

This program utilizes the following tools:

1. **Fiji** 
    This project used Fiji for preprocessing into preprocess image stacks as contrast enhancement, filtering, and particle analysis.

    [Fiji](https://fiji.sc/) is an open-source distribution of ImageJ focusing on image analysis. 
    
    - Repository: [Fiji](https://github.com/fiji/fiji)  
    - License: [GPL License](https://imagej.net/licensing/)

2. **StarDist**
    In this project, the standard StarDist model was employed to generate high-quality nuclei masks from image data, significantly improving segmentation accuracy and reducing background noise issues commonly encountered in immunofluorescence (IF) image analysis.
    
    [StarDist](https://stardist.net/)

    - Repository: [StarDist](https://github.com/stardist/stardist)  
    - License: [BSD 3-Clause License](https://github.com/stardist/stardist/blob/main/LICENSE.txt)


## Contributors

- [Aleksandr Dolskii](aleksandr.dolskii@fccc.edu)

- [Ekaterina Shitik](mailto:shitik.ekaterina@gmail.com) 

Enjoy your use 💫

