import argparse
import os.path
import subprocess
import glob
import sys
import json
import MDAnalysis as mda
from io import StringIO
from rdkit import Chem
from schrodinger import adapter
from schrodinger.rdkit import rdkit_adapter
from schrodinger.structutils import analyze
from backmapping.config import load_config, check_fields
from backmapping.logger import logger, check_if_file_exists
from backmapping.stereo.utils import get_PDBBlock, get_stereo


def make_tpr_file(basename: str, wdir, mdp, gmx):
    box_args = [
        f"{gmx}",
        "editconf",
        "-f",
        f"{basename}_ini.pdb",
        "-o",
        f"{basename}_box.gro",
        "-d",
        "1",
        "-c",
        "-bt",
        "cubic",
    ]
    subprocess.run(box_args, cwd=wdir)

    if not check_if_file_exists(os.path.join(wdir, f"{basename}_box.gro")):
        sys.exit()

    if not check_if_file_exists(os.path.join(wdir, f"{basename}.top")):
        sys.exit()

    tpr_args = [
        f"{gmx}",
        "grompp",
        "-f",
        f"{mdp}",
        "-c",
        f"{basename}_box.gro",
        "-p",
        f"{basename}.top",
        "-o",
        f"{basename}.tpr",
        "-maxwarn",
        "1",
    ]
    subprocess.run(tpr_args, cwd=wdir)

    if not check_if_file_exists(os.path.join(wdir, f"{basename}.tpr")):
        sys.exit()

    pdb_args = [
        f"{gmx}",
        "editconf",
        "-f",
        f"{basename}.tpr",
        "-o",
        f"{basename}_box.pdb",
        "-conect",
    ]
    subprocess.run(pdb_args, cwd=wdir)

    if not check_if_file_exists(os.path.join(wdir, f"{basename}_box.pdb")):
        sys.exit()

    return os.path.join(wdir, f"{basename}.tpr")


def get_stereo_labels(filepath, config):
    """
    Get stereoisomer lables for this molecule.
    """
    logger.debug(f"Getting stereo labels for file: {filepath}")

    # Split off the "_ini.pdb" part of the filename to get the usable basename
    basename = os.path.splitext(os.path.basename(filepath))[0]
    if not "_" in basename:
        logger.error(f"Couldn't get basename for file: {filepath}")
        return
    else:
        basename = basename.split("_")[0]

    tpr = make_tpr_file(
        basename,
        config["filepaths"]["aa_gmx_parameters"],
        config["filepaths"]["mdp"],
        config["filepaths"]["gmx"],
    )

    u = mda.Universe(tpr, filepath)

    PDBBlock = get_PDBBlock(u, atom_group=u.select_atoms("all"), frame=0)

    if PDBBlock == None:
        logger.error(f"Couldn't make PDBBlock for: {basename}")
        return

    stereo_reference = {basename: get_stereo(PDBBlock)}

    logger.info("Stereoisomer record: ", stereo_reference)

    with open(
        os.path.join(config["filepaths"]["aa_gmx_parameters"], f"{basename}.json"),
        "w",
    ) as f:
        json.dump(stereo_reference, f, indent=4)


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
        "aa_gmx_parameters",
        "mdp",
        "gmx",
    ]
    check_fields(config["filepaths"], required_paths)

    logger.info("Generating stereoisomer reference from coordinate files")
    coord_files = glob.glob(
        os.path.join(config["filepaths"]["aa_gmx_parameters"], "*ini.pdb")
    )

    for file_path in coord_files:
        get_stereo_labels(file_path, config)

    # Combine the stereoisomer records into a single file
    stereo = {}
    stereo_files = glob.glob(
        os.path.join(config["filepaths"]["aa_gmx_parameters"], "*json")
    )
    for stereo_file in stereo_files:
        with open(stereo_file, "r") as f:
            stereo.update(json.load(f))

    with open(os.path.join(config["filepaths"]["stereo"], "stereo.json"), "w") as f:
        json.dump(stereo, f, indent=4)


if __name__ == "__main__":
    logger.info(f"Starting {__file__}")
    main()
    logger.info(f"Finished {__file__}")
