import argparse
import os.path
import subprocess
import glob
import sys
from backmapping.config import load_config, check_fields
from backmapping.logger import logger


def truncate_atom_names(param_file: str, basename, config):
    logger.debug(f"Truncating atom names for param file: {param_file}")
    lines = []
    atomnames = []

    # Collect atom names to shorten
    with open(param_file, "r") as f:
        for line in f.readlines():
            lines.append(line)
            if (
                (len(line.split()) > 2)
                and (line.split()[0] == "ATOM")
                and (len(line.split()[1]) > 4)
            ):
                atomnames.append(line.split()[1])

    # Map original names to shortened versions
    # Shortened names take the original first character and
    # replace the final three with a right-padded integer
    shortened_names = {}
    for i, name in enumerate(atomnames):
        short = name[0] + str("{:>03}".format(i))
        n_spaces = len(name) - len(short)
        shortened_names[name] = short + " " * n_spaces

    # Write out parameter file with shortened atom names
    truncated_param_file = os.path.join(
        config["filepaths"]["aa_gmx_parameters"], basename + ".str"
    )
    with open(truncated_param_file, "w") as f:
        for line in lines:
            for name in atomnames:
                line = line.replace(name, shortened_names[name])
            f.write(line)

    if not os.path.isfile(truncated_param_file):
        logger.error(f"Couldn't find truncated parameter file: {truncated_param_file}")
        sys.exit()
    else:
        return truncated_param_file


def convert_params(param_file: str, coord_file: str, basename: str, config: dict):
    """
    Convert CHARMM-formatted CGenFF params to GMX-formatted.
    """
    logger.debug(f"Converting parameter file: {param_file}")

    # Truncate all atom names in the CGenFF STR file to 4 characters
    truncated_param_file = truncate_atom_names(param_file, basename, config)

    # Run the charmm2gmx conversion
    subprocess.run(
        args=[
            "python",
            config["filepaths"]["charmm2gmx"],
            basename,
            coord_file,
            truncated_param_file,
            config["filepaths"]["charmm36_ff"],
        ],
        shell=False,
        cwd=config["filepaths"]["aa_gmx_parameters"],
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
        "aa_gmx_parameters",
        "charmm2gmx",
        "charmm36_ff",
    ]
    check_fields(config["filepaths"], required_paths)

    logger.info("Converting CHARMM format parameters to GMX format")
    charmm_param_files = glob.glob(
        os.path.join(config["filepaths"]["aa_cgenff_parameters"], "*str")
    )

    for param_file in charmm_param_files:
        # Check if the corresponding coordinate file exists
        basename = os.path.splitext(os.path.basename(param_file))[0]
        coord_file = os.path.join(config["filepaths"]["aa_coords"], basename + ".mol2")
        if not os.path.isfile(coord_file):
            logger.error(
                f"Couldn't find corresponding coordinate file: {coord_file} for parameter file: {param_file}"
            )
        else:
            convert_params(param_file, coord_file, basename, config)


if __name__ == "__main__":
    logger.info(f"Starting {__file__}")
    main()
    logger.info(f"Finished {__file__}")
