import argparse
import os.path
import subprocess
import glob
import polars as pl
from backmapping.config import load_config, check_fields
from backmapping.logger import logger


def run_martini_itp_builder(data, config):
    subprocess.run(
        args=[
            "python2",
            config["filepaths"]["write_martini_itp"],
            "-alhead",
            str(data["martini_head"]),
            "-allink",
            str(data["martini_linker"]),
            "-altail",
            str(data["martini_tailA"]) + " " + str(data["martini_tailB"]),
            "-alname",
            str(data["lipid_code"]),
            "-o",
            str(data["lipid_code"]) + ".itp",
        ],
        cwd=config["filepaths"]["cg_martini_parameters"],
    )


def generate_params(filepath, config):
    """
    Generate MARTINI params using the bead labels in the provided csv file.

    This function currently filters MARTINI bead assignments to only those
    species which contain a linker, tailA, and tailB beads. This is to avoid
    attempting to generate ITP files for lipids without those beads (i.e., sterols).
    """
    logger.debug(f"Reading MARTINI bead assignment file: {filepath}")

    df = pl.read_csv(config["filepaths"]["martini_beads"]).filter(
        pl.col("martini_linker").is_not_null(),
        pl.col("martini_tailA").is_not_null(),
        pl.col("martini_tailB").is_not_null(),
    )

    for row in df.iter_rows(named=True):
        run_martini_itp_builder(row, config)

    logger.debug(f"Finished with MARTINI bead file.")


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

    required_paths = ["martini_beads", "cg_martini_parameters", "write_martini_itp"]
    check_fields(config["filepaths"], required_paths)

    logger.info("Generating MARTINI parameters.")

    generate_params(config["filepaths"]["martini_beads"], config)


if __name__ == "__main__":
    logger.info(f"Starting {__file__}")
    main()
    logger.info(f"Finished {__file__}")
