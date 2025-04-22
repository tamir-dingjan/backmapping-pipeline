# backmapping-pipeline

## Table of Contents
- [Intro](#intro)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Common issues and how to fix them](#common-issues-and-how-to-fix-them)
- [License](#license)

## Intro
`backmapping-pipeline` is an organized set of tools for generating all-atom parameters for lipid molecules to allow resolution transformation from a Martini coarse-grained simulation.

## Requirements
* Python 3.11
* A local installation of the Schrodinger Software Suite (this is used for controlling stereochemistry)
* CHARMM-to-GROMACS format conversion tool (forked by me [here](https://github.com/tamir-dingjan/cgenff_charmm2gmx/blob/main/cgenff_charmm2gmx_py3_nx3.py) to work with NetworkX version 3.4.2; note the original is available from the Lemkul lab repository [here](https://github.com/Lemkul-Lab/cgenff_charmm2gmx/blob/main/cgenff_charmm2gmx_py3_nx2.py))
* *(Recommended)* CGenFF binary. This is proprietary software available from [SilcsBio](https://app.cgenff.com/homepage), used to generate parameter files for lipid molecules. Access to a local binary file allows the parameters to be generated using `backmapping-pipeline`, however if this tool is inaccessible it is also possible to provide your own parameter files and simply skip the step which requires CGenFF.

## Installation
1. Install the required packages by running `pip install -r requirements.txt` in your terminal.
2. Clone this repository using the URL

## Usage
The backmapping pipeline follows these steps:

1. Preparation
    - Create `config.yaml` with all directory paths and binary locations
        - Note that the `filepaths` section of `config.yaml` contains paths for directories and for binaries. The directory locations are defined relative to the location of the file `config.yaml`. These directories should be empty at this point, and the user needs write permissions to them:
            - `aa_coords`: Directory for the all-atom lipid coordinate files. Default value: "data/aa_coords"
            - `aa_cgenff_parameters`: Directory for the CGenFF parameter files. Default value: "data/aa_cgenff"
            - `aa_gmx_parameters`: Directory for the GROMACS format parameter files. Default value: "data/aa_gmx"
            - `stereo`: Directory for the concatenated stereoconformation record file. Default value: "data/stereo"
            - `mapping`: Directory for the mapping files to be written to. Default value: "data/Mapping"
            - `cg_martini_parameters`: Directory for the coarse-grained Martini parameter files. Default value: "data/cg_martini"
        - The locations of binaries or external Python tools:
            - `cgenff`: Location of the CGenFF binary. Default value: "/opt/silcsbio.2024.1/cgenff/cgenff"
            - `charmm2gmx`: Location of the CHARMM-to-GROMACS parameter file conversion tool. Default value: "/opt/charmm2gmx/cgenff_charmm2gmx/cgenff_charmm2gmx_py3_nx3.py"
            - `charmm36_ff`: Location of the CHARMM36 force field in GROMACS format, available from the [MacKerell group website](https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs). Default value: "/usr/local/gromacs/share/gromacs/top/charmm36-jul2022.ff"
            - `gmx`: Location of the GROMACS binary `gmx` or `gmx_mpi`. Default value: "/usr/local/gromacs-2021.3/bin/gmx_mpi"
            - `write_martini_itp`: Location of the Martini parameter file generation tool. Default value: "/opt/martini/lipid-martini-itp-v06.py"
    - Assemble all-atom structures of the lipids to be simulated in mol2 format, stored in the path defined in `config.yaml` under `filepaths:aa_coords`
    - Create a virtual environment using the Schrodinger suite and install MDAnalysis:
    ```
    $SCHRODINGER/run schrodinger_virtualenv.py schro.ve
    source schro.ve/bin/activate
    python3 -m pip install --upgrade pip
    pip install --upgrade MDAnalysis
    ```


2. Parameterise lipids using CGenFF: `python scripts/prepare/0_generate_cgenff_params.py --config config.yaml`

3. Convert CHARMM-formatted parameters to GMX format: `python scripts/prepare/1_convert_cgenff_params_to_gmx.py --config config.yaml`

4. Generate stereochemistry reference: `source schro.ve/bin/activate; python scripts/prepare/2_make_stereo_reference.py --config config.yaml`


## Common issues and how to fix them
- The all-atom structures of lipids in mol2 format should have reasonable geometry, because the way bond valences are inferred by RDKit does involve the molecular conformation. If using the Schrodinger Suite to prepare these input structures, a quick minimisation seems to suffice.
- The Martini force field parameters for sterols are not able to be generated automatically using the tools available from the force field developers. However, `backmapping-pipeline` still needs to work with all-atom coordinate files and CGenFF-generated all-atom force field parameters in order to control the stereochemistry for sterols during backmapping. To make this possible, please add your own Martini force field parameters for any sterol species into the `filepaths:cg_martini` directory, and ensure the filenames match with the all-atom coordinate files.


## License
This project is licensed under the MIT License.