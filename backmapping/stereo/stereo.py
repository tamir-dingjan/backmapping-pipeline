import os
import argparse
import numpy as np
from backmapping.logger import logger
from backmapping.stereo.io import load_coordinates, load_stereo_lookup
from backmapping.stereo.utils import get_PDBBlock, get_stereo
from backmapping.stereo.corrections import correct_stereo, InvalidPatchStereo


def run(coord, topol, stereo_lookup):
    # Check files exist
    for f in (coord, topol, stereo_lookup):
        if not os.path.isfile(f):
            logger.error(f"Couldn't find file: {f}")
            raise FileNotFoundError

    # Load the coordinates
    u = load_coordinates(coord, topol)

    # Get the original positions only once, and mutate for each residue
    positions = u.atoms.positions.copy()

    stereo_lookup = load_stereo_lookup(stereo_lookup)

    # Check stereoconfiguration for each residue in the patch
    resids = np.unique(u.select_atoms("all").resids)
    for res in resids:
        logger.debug(f"Check stereoconfiguration for resid: {res}")
        resname = u.select_atoms(f"resid {str(res)}").resnames[0]
        pdb = get_PDBBlock(u, u.select_atoms(f"resid {str(res)}"), frame=0)

        stereo = get_stereo(pdb)
        stereo_lookup_res = stereo_lookup[resname]

        # Correct any mismatches between the reference and measured stereoconformation
        try:
            positions = correct_stereo(stereo, stereo_lookup_res, u, res, positions)
        except InvalidPatchStereo as e:
            logger.warning("Invalid patch stereo detected.")
            raise (e)

    u.atoms.positions = positions
    u.select_atoms("all").write("stereo.gro")
    logger.info("Finished stereoconformer correction")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform stereochemistry stamping")
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
