# Get the patch stereo from the coordinate file
# Compare with the lookup
# Save to disk whether the patch is correct


import os
import argparse
import numpy as np
from backmapping.logger import logger
from backmapping.stereo.io import load_coordinates, load_stereo_lookup
from backmapping.stereo.utils import get_PDBBlock, get_stereo
from backmapping.stereo.corrections import correct_stereo
from backmapping.stereo.validation import (
    validate_patch_stereoconformation,
    InvalidPatchStereo,
    get_nonmatching_unsat,
    get_nonmatching_chiral_centers,
    get_nonmatching_amide,
)


def run(coord, topol, stereo_lookup):
    # Check files exist
    for f in (coord, topol, stereo_lookup):
        if not os.path.isfile(f):
            logger.error(f"Couldn't find file: {f}")
            raise FileNotFoundError

    # Load the coordinates
    u = load_coordinates(coord, topol)

    stereo_lookup = load_stereo_lookup(stereo_lookup)

    # Check stereoconfiguration for each residue in the patch
    resids = np.unique(u.select_atoms("all").resids)
    for res in resids:
        logger.debug(f"Check stereoconfiguration for resid: {res}")
        resname = u.select_atoms(f"resid {str(res)}").resnames[0]
        pdb = get_PDBBlock(u, u.select_atoms(f"resid {str(res)}"), frame=0)

        stereo = get_stereo(pdb)
        stereo_lookup_res = stereo_lookup[resname]

        # Validate that the chiral centers, unsaturations, and amide bonds
        # defined in the stereo lookup are also detected in the patch

        if not validate_patch_stereoconformation(stereo, stereo_lookup_res):
            logger.warning(
                "Patch stereoconformation is not valid. Cannot correct stereoconformation."
            )
            raise InvalidPatchStereo

        # Check for any mismatches between the reference and measured stereoconformation
        if (
            (get_nonmatching_amide(stereo, stereo_lookup_res) == [])
            and (get_nonmatching_unsat(stereo, stereo_lookup_res) == [])
            and (get_nonmatching_chiral_centers(stereo, stereo_lookup_res) == [])
        ):
            stereomatch = True
        else:
            stereomatch = False

        # Save the result of the check to disk - this allows reading beyond the subprocess enrivonment
        with open("stereo.check", "w") as f:
            f.write(f"{stereomatch}")

    logger.info("Finished checking patch stereoconformation.")
    logger.info(f"Matched stereochemistry: {stereomatch}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check if the stereoconformation in the provided coordinate file and lookup match"
    )
    parser.add_argument("coord", type=str, help="Path to coordinate file")
    parser.add_argument("topol", type=str, help="Path to topology file")
    parser.add_argument(
        "stereo", type=str, help="Path to the stereochemistry lookup reference"
    )

    args = parser.parse_args()

    coord = args.coord
    topol = args.topol
    stereo_lookup = args.stereo

    run(coord, topol, stereo_lookup)
