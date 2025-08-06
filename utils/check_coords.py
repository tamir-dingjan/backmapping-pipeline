import argparse
import os.path
import glob
from backmapping.config import load_config, check_fields
from backmapping.logger import logger
import polars as pl


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
    ]
    check_fields(config["filepaths"], required_paths)

    logger.info(
        "Checking that all lipids defined in lipidomics datafile are present in aa_coords"
    )
    coord_files = glob.glob(os.path.join(config["filepaths"]["aa_coords"], "*mol2"))

    df = pl.read_csv(config["filepaths"]["martini_beads"])
    defined_lipids = df.select("lipid_code").to_series().to_list()

    lipid_structures = [os.path.splitext(os.path.basename(x))[0] for x in coord_files]

    for defined_lipid in defined_lipids:
        if not (defined_lipid in lipid_structures):
            logger.warning(
                f"Lipid defined in lipidomics data file ({defined_lipid}) not found in lipid structures!"
            )

    for lipid_structure in lipid_structures:
        if not (lipid_structure in defined_lipids):
            logger.warning(
                f"Lipid structure file ({lipid_structure}) not defined in lipidomics data file!"
            )


if __name__ == "__main__":
    logger.info(f"Starting {__file__}")
    main()
    logger.info(f"Finished {__file__}")
