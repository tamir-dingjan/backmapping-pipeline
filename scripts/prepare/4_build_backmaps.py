import argparse
import os.path
import glob
from rdkit import Chem
from backmapping.config import load_config, check_fields
from backmapping.logger import logger
from backmapping.workflow import BackmappingWorkflow


def build_mapping(filepath, config):
    """
    Build a mapping file for the provided mol2 file.
    This will attempt to load the corresponding all-atom and coarse-grained
    topology files.
    """
    logger.info(f"Generating mapping file for: {filepath}")
    basename = os.path.splitext(os.path.basename(filepath))[0]
    mapping_file = os.path.join(
        config["filepaths"]["mapping"], basename + ".charmm36.map"
    )

    # Check if coord and parameter files exist

    # PDB
    aa_coord = os.path.join(
        config["filepaths"]["aa_gmx_parameters"], basename + "_box.pdb"
    )

    # mol2
    # aa_coord = filepath

    aa_param = os.path.join(config["filepaths"]["aa_gmx_parameters"], basename + ".itp")
    cg_param = os.path.join(
        config["filepaths"]["cg_martini_parameters"], basename + ".itp"
    )
    stereo = os.path.join(config["filepaths"]["aa_gmx_parameters"], basename + ".json")

    if not os.path.isfile(aa_param):
        logger.error(f"Couldn't find all-atom topology file: {aa_param}")
        raise FileNotFoundError
    if not os.path.isfile(cg_param):
        logger.error(f"Couldn't find coarse-grained topology file: {cg_param}")
        raise FileNotFoundError

    backmapper = BackmappingWorkflow(
        aa_coord, basename, aa_param, cg_param, stereo, mapping_file
    )

    backmapper.set_lipid_class()

    # Check for supported lipid class before proceeding
    if backmapper.is_supported_lipid_class:
        backmapper.assign_regions()
        backmapper.allocate_head()
        backmapper.allocate_linker()
        backmapper.allocate_both_tails()
        backmapper.catch_unmapped_atoms()
        backmapper.write_mapping_file()


def main():

    parser = argparse.ArgumentParser(
        description="Load config file with file input and output paths."
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",
        help="Path to the YAML config file.",
    )
    args = parser.parse_args()

    config = load_config(args.config)

    required_paths = [
        "aa_coords",
        "aa_gmx_parameters",
        "cg_martini_parameters",
        "stereo",
        "mapping",
    ]
    check_fields(config["filepaths"], required_paths)

    logger.info("Building mapping files for lipids")
    coord_files = glob.glob(os.path.join(config["filepaths"]["aa_coords"], "*mol2"))

    for file_path in coord_files:
        build_mapping(file_path, config)


if __name__ == "__main__":
    logger.info(f"Starting {__file__}")
    main()
    logger.info(f"Finished {__file__}")
