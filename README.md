# backmapping-pipeline

## Table of Contents
- [Intro](#intro)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Common issues and how to fix them](#common-issues-and-how-to-fix-them)
- [Notes on cholesterol](#notes-on-cholesterol)
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

This has only been tested with Schrodinger Suite version 2024-3.

## Usage
The backmapping pipeline follows these steps:

1. Preparation
    - Create `config.yaml` with all directory paths and binary locations
        - Note that the `filepaths` section of `config.yaml` contains paths for directories and for binaries. The directory locations are defined relative to the location of the file `config.yaml`. These directories should be empty at this point, and the user needs write permissions to them:
            - `aa_coords`: Directory for the all-atom lipid coordinate files. Default value: `data/aa_coords`
            - `aa_cgenff_parameters`: Directory for the CGenFF parameter files. Default value: `data/aa_cgenff`
            - `aa_gmx_parameters`: Directory for the GROMACS format parameter files. Default value: `data/aa_gmx`
            - `stereo`: Directory for the concatenated stereoconformation record file. Default value: `data/stereo`
            - `mapping`: Directory for the mapping files to be written to. Default value: `data/Mapping`
            - `cg_martini_parameters`: Directory for the coarse-grained Martini parameter files. Default value: `data/cg_martini`
        - The locations of binaries or external Python tools:
            - `cgenff`: Location of the CGenFF binary. Default value: `/opt/silcsbio.2024.1/cgenff/cgenff`
            - `charmm2gmx`: Location of the CHARMM-to-GROMACS parameter file conversion tool. Default value: `/opt/charmm2gmx/cgenff_charmm2gmx/cgenff_charmm2gmx_py3_nx3.py`
            - `charmm36_ff`: Location of the CHARMM36 force field in GROMACS format, available from the [MacKerell group website](https://mackerell.umaryland.edu/charmm_ff.shtml#gromacs). Default value: `/usr/local/gromacs/share/gromacs/top/charmm36-jul2022.ff`
            - `gmx`: Location of the GROMACS binary `gmx` or `gmx_mpi`. Default value: `/usr/local/gromacs-2021.3/bin/gmx_mpi`
            - `write_martini_itp`: Location of the Martini parameter file generation tool. Default value: `/opt/martini/lipid-martini-itp-v06.py`
    - Assemble all-atom structures of the lipids to be simulated in mol2 format, stored in the path defined in `config.yaml` under `filepaths:aa_coords`
    - Create a virtual environment using the Schrodinger suite and install MDAnalysis:
    ```
    $SCHRODINGER/run schrodinger_virtualenv.py schro.ve
    source schro.ve/bin/activate
    python3 -m pip install --upgrade pip
    pip install --upgrade MDAnalysis
    ```
    - Create a second virtual environment to install other required packages:
    ```
    python -m venv ./venv
    source ./venv/bin/activate
    uv pip install pyyaml numpy networkx polars rdkit mdtraj scikit-learn
    ```

Prior to each step in the protocol, the `PYTHONPATH` must be set prior to running the Python script to allow the Python itnerpreter to find the `backmapping` module:
`export PYTHONPATH=$(pwd)`


2. Parameterise lipids using CGenFF: `python scripts/prepare/0_generate_cgenff_params.py --config config.yaml`

3. Convert CHARMM-formatted parameters to GMX format: `python scripts/prepare/1_convert_cgenff_params_to_gmx.py --config config.yaml`

4. Generate stereochemistry reference: `source schro.ve/bin/activate; python scripts/prepare/2_make_stereo_reference.py --config config.yaml`

5. Exit the Schrodinger virtual environment and reactivate the `./venv` environment, including resetting the `PYTHONPATH`.Generate Martini parameters: `python scripts/prepare/3_generate_martini_params.py --config config.yaml`. Ensure you have a `python2` in your path for this step, since the Martini ITP generator script is run in a subprocess using `python2`.

6. Build mapping files: `python scripts/prepare/4_build_backmaps.py --config config.yaml`

7. Generate backmapped patches: `python scripts/run/5_make_patches.py --config config.yaml`

## Common issues and how to fix them
- The all-atom structures of lipids in mol2 format should have reasonable geometry, because the way bond valences are inferred by RDKit does involve the molecular conformation. If using the Schrodinger Suite to prepare these input structures, a quick minimisation seems to suffice.
- The Martini force field parameters for sterols are not able to be generated automatically using the tools available from the force field developers. However, `backmapping-pipeline` still needs to work with all-atom coordinate files and CGenFF-generated all-atom force field parameters in order to control the stereochemistry for sterols during backmapping. To make this possible, please add your own Martini force field parameters for any sterol species into the `filepaths:cg_martini` directory, and ensure the filenames match with the all-atom coordinate files.

## Notes on cholesterol
As noted above, `backmapping-pipeline` relies on CGenFF-generated all-atom force field parameters to generate topology files which are used to record the stereoconformation of each lipid molecule. This does not mean that each all-atom patch simulation must use the CGenFF parameters, however. For cholesterol (and also for any other lipid molecule), any parameter file can be used as long as it matches the atom naming and ordering in the coordinate file in `filepaths:aa_coords`.
For cholesterol specifically, the current recommendation is to use the following files:
    - Coordinate file available at [CHARMM-GUI](https://www.charmm-gui.org/archive/csml/chl1.pdb)
    - Mapping file available at the [Martini site](https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini2/lipidome/sterols/chol/CHOL.amber.map)
    - Martini parameter file also available at the [Martini site](https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini2/lipidome/sterols/chol/martini_v2.0_CHOL_02.itp)

This mapping file is for the AMBER forcefield, used because the atom namping is a closer match for the Gromacs port of CHARMM36 July 2022 force field. A few atom names are ordered differently in the mapping file, and need to be updated to match the Gromacs port (this just means moving the atom lines in the `[ atoms ]` section of the mapping file to match the Gromacs port order):

```
        AMBER map       GROMACS port of CHARMM36
2       O3              H3
3       H3'             O3
4       H3              H3'
```

To use these files for cholesterol in the `backmapping-pipeline`:
- Download the coordinate, mapping, and Martini parameter files
- Rename the coordinate file residue name to match the CHARMM36 forcefield cholesterol, i.e., to `CHL1`
- Generate simulation parameters using the Gromacs version specified in `config.yaml` under `filepaths:gmx`:
    - `gmx pdb2gmx -f chl1.pdb -o chl1.gro`
    - Select the CHARMM forcefield specificied in `config.yaml`
    - Extract the cholesterol topology and parameter components from the generated `topol.top` to `CHL1.itp`, and rename the molecule name in the `[ moleculetype ]` section to `CHL1`. The relevant sections are:
        - `moleculetype`, `atoms`, `bonds`, `pairs`, `angles`, `dihedrals`
- Convert `chl1.pdb` to `CHL1.mol2`, and change the molecule name from `chl1` to `CHL1`
    - Ensure the correct bond order for the double bond linking atoms `C5` and `C6`
- Edit the mapping file:
    - Rename `CHOL` to `CHL1` in the `[ molecule ]` section
    - Replace the mapping forcefield labels with `charmm36`
    - Adjust atom ordering in the `[ atoms ]` section to match the Gromacs CHARMM36 forcefield atom order in `CHL1.mol2`:
        - Atom 2 change from `O3` to `H3`
        - Atom 3 change from `H3'` to `O3`
        - Atom 4 change from `H3` to `H3'`
    - Save as `CHL1.charmm36.map`
- Edit the Martini parameter file `martini_v2.0_CHOL_02.itp` to change all occurances of `CHOL` to `CHL1` in the `[ moleculetype ]` and `[ atoms ]` sections
- Copy the coordinate file `CHL1.mol2` to `filepaths:aa_coords`
- Complete steps 2 through 5 in [Usage](#usage) to generate CGenFF parameters, stereoconformer references, Martini parameters and mapping files for all lipids. This will fail to generate Martini parameters and a mapping file for cholesterol.
- Copy the Martini parameter file for cholesterol `martini_v2.0_CHOL_02.itp` to `filepaths:cg_martini_parameters/CHL1.itp`
- Copy the Charmm parameter file `CHL1.itp` to `filepaths:aa_gmx_parameters/CHL1.itp`, overwriting the CGenFF-generated file
- Copy the mapping file `CHL1.charmm36.map` to `filepaths:mapping/CHL1.charmm36.map`
- Complete step 6 to build backmaps.
- Proceed to the final step 7 to generate backmapped patches:
    - If the coarse-grained trajectory uses the residue name `CHOL` for cholesterol, rename this to `CHL1` in the coordinate file specified in `patches:top`

## License
This project is licensed under the MIT License.