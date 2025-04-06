import argparse
import os.path
import subprocess
import glob
from backmapping.config import load_config, check_fields
from backmapping.logger import logger


def generate_params(filepath, config):
    """
    Generate CGenFF params for the provided mol2 file.
    """
    logger.debug(f"Generating parameter for file: {filepath}")
    basename = os.path.splitext(os.path.basename(filepath))[0]
    with open(
        os.path.join(config["filepaths"]["aa_cgenff_parameters"], basename + ".str"),
        "w",
    ) as f:

        subprocess.call(
            args=[
                config["filepaths"]["cgenff"],
                filepath,
            ],
            stdout=f,
        )


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
        "aa_cgenff_parameters",
        "cgenff",
    ]
    check_fields(config["filepaths"], required_paths)

    logger.info("Generating CGenFF parameters for coordinate files")
    coord_files = glob.glob(os.path.join(config["filepaths"]["aa_coords"], "*mol2"))

    for file_path in coord_files:
        generate_params(file_path, config)


if __name__ == "__main__":
    logger.info(f"Starting {__file__}")
    main()
    logger.info(f"Finished {__file__}")
