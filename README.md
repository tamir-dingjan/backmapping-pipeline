# backmapping-pipeline

## Table of Contents
- [Intro](#intro)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)

## Intro
`backmapping-pipeline` is an organized set of tools for generating all-atom parameters for lipid molecules to allow resolution transformation from a MARTINI coarse-grained simulation.

## Requirements
* Python 3.11
* A local installation of the Schrodinger Software Suite (this is used for controlling stereochemistry)
* CHARMM-to-GMX format conversion tool (forked by me [here](https://github.com/tamir-dingjan/cgenff_charmm2gmx/blob/main/cgenff_charmm2gmx_py3_nx3.py) to work with NetworkX version 3.4.2; note the original is available from the Lemkul lab repository [here](https://github.com/Lemkul-Lab/cgenff_charmm2gmx/blob/main/cgenff_charmm2gmx_py3_nx2.py))

## Installation
1. Install the required packages by running `pip install -r requirements.txt` in your terminal.
2. Clone this repository using the URL

## Usage
The backmapping pipeline follows these steps:

1. Preparation
    - Create `config.yaml` with all file paths and binary locations
    - Assemble all-atom structures of the lipids to be simulated in mol2 format, stored in the path defined in the configuration file under filepaths > aa_coords
    - Create a virtual environment using the Schrodinger suite and install MDAnalysis:
    ```
    $SCHRODINGER/run schrodinger_virtualenv.py schro.ve
    source schro.vs/bin/activate
    python3 -m pip install --upgrade pip
    pip install --upgrade MDAnalysis
    ```


2. Parameterise lipids using CGenFF: `python scripts/prepare/0_generate_cgenff_params.py --config config.yaml`

3. Convert CHARMM-formatted parameters to GMX format: `python scripts/prepare/1_convert_cgenff_params_to_gmx.py --config config.yaml`

4. Generate stereochemistry reference: 


## License
This project is licensed under the MIT License.