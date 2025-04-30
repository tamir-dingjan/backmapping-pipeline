#!/usr/bin/python

# kicker.py
# A tool to give random kicks to atoms which are too close to eachother.
#
# Writes the kicked output coordinate file to the same directory as the
# input file, titled "kicked.gro".
#
# Usage:
# $> python kicker.py file.gro

import numpy as np
import scipy
from sklearn.neighbors import BallTree
import sys
import os
import argparse
import logging

logging.basicConfig(filename="kicker.log", level=logging.DEBUG)


def coordload(filepath, return_box_vectors=False):
    # Returns a numpy array of size N, 3 where N is the number of atoms
    coords = []
    with open(filepath, "r") as infile:
        # Get the file length here so we can identify the last line when reading
        linecount = sum(1 for line in infile) - 1
        infile.seek(0)

        for line_i, line in enumerate(infile.readlines()):
            if (line_i == 0) or (line_i == 1):
                continue
            elif line_i == linecount:
                box_vectors = np.asarray(line.split()[:3], dtype=np.float32)
            else:
                x = np.float32(line[20:28])
                y = np.float32(line[28:36])
                z = np.float32(line[36:44])
            coords.append([x, y, z])
    if return_box_vectors:
        return np.asarray(coords), box_vectors
    else:
        return np.asarray(coords)


def get_kickable_atoms_sklearn(coords, radius):
    """
    Find overlapping atoms using sklearn's BallTree to search for neighbours in coordinate space.

    See Scikit-Learn documentation at: https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.BallTree.html

    The radius in the BallTree.query_radius function is the distance in coordinate space used to check for nearby atom pairs.

    Parameters:
        coords (numpy.ndarray): The coordinates of the atoms.

    Returns:
        numpy.ndarray: The indices of the overlapping atoms.

    """
    # Find overlapping atoms using sklearn's BallTree
    tree = BallTree(coords)
    inds, dist = tree.query_radius(coords, r=radius, return_distance=True)
    # Are there any indices at distances greater than 0?
    kickable = np.unique(
        (np.concatenate(inds).ravel()[np.concatenate(dist).ravel() > 0])
    )
    return kickable


def get_pbc_overlapping_atoms(coords, box_vectors, radius):
    # Find the points in coords which overlap with
    # the PBC-transform versions of coords.
    # We can do this by making a KDTree of the PBC coords
    # and querying the central, real coords.

    # Make PBC transformed coords
    # Check the +Z, -Z, +X, -X, +Y, and -Y transformed coords
    # Also want to consider the double-axis diagonal neighbors.
    #
    # The triple axis neighbors are probably safe to neglect,
    # because this would involve a coordinate crossing the box vertex,
    # which is unlikely in bilayer systems with a tall water layer

    dims = (0, 1, 2)
    directions = (-1, 1)
    PBC_nbors = []
    seen = set()

    for dim1 in dims:

        for direction1 in directions:
            for dim2 in dims:
                if dim1 == dim2:
                    # Do the single-axis neighbors
                    PBC_transform = [0, 0, 0]
                    PBC_transform[dim1] = box_vectors[dim1] * direction1

                    transform_id = "-".join(str(x) for x in PBC_transform)

                    if not transform_id in seen:
                        seen.add(transform_id)
                        PBC_nbors.append(coords + PBC_transform)
                        continue

                for direction2 in directions:
                    # Do diagonal neighbors
                    PBC_transform = [0, 0, 0]
                    PBC_transform[dim1] = box_vectors[dim1] * direction1
                    PBC_transform[dim2] = box_vectors[dim2] * direction2

                    transform_id = "-".join(str(x) for x in PBC_transform)

                    if not transform_id in seen:
                        seen.add(transform_id)
                        PBC_nbors.append(coords + PBC_transform)

    PBC_nbors = np.concatenate(PBC_nbors)

    # Make KDTree out of PBC coords
    tree = BallTree(PBC_nbors)

    # Query the tree using the real coords
    inds, dist = tree.query_radius(coords, r=radius, return_distance=True)

    # The inds values here are indices for the PBC neighbors.
    # Their position in the list indicates which of the coords is overlapping.
    # Most of the elements of inds are empty arrays, indicating no overlaps
    # We want to select the positions of inds which are non-empty arrays
    kickable = np.where(np.asarray([np.size(x) for x in inds]) != 0)[0]

    return kickable


def kick(coords, kickable, kick_amount):
    # Get the indices for overlapping atoms
    # Generate kick array of same shape as coords
    kick_array = np.random.default_rng().uniform(
        -1 * kick_amount, kick_amount, coords.shape
    )
    # Set the kick value for all non-kickable entries to 0
    mask = np.ones(len(coords), bool)
    mask[kickable] = 0
    kick_array[mask] = 0
    # Sum arrays
    coords += kick_array


def write_outfile(filepath, coords, outpath):
    """
    Writes the modified coordinates to a new file.

    Args:
        filepath (str): The path to the input file.
        coords (numpy.ndarray): The modified coordinates.
        outpath (str): The path to the output file.

    Returns:
        None

    This function reads the input file line by line and writes the modified coordinates to the output file. It skips the first two lines and the final line of the input file and writes them directly to the output file. For coordinate lines, it uses the modified coordinates in the `coords` array. The first 20 characters of the line are the atom and residue labels, so they are transferred unchanged to the new line.

    Note:
        The input file should have the same number of lines as the `coords` array.

    """
    # Save out the modified coords
    # Step through the input file
    with open(filepath, "r") as infile:
        linecount = sum(1 for line in infile) - 1
        infile.seek(0)

        with open(outpath, "w") as outfile:
            # Transmit the first two lines and the final line directly to output
            # For coordinate lines, use the modified coordinates in coords
            # The first 20 characters of the line are the atom and residue labels,
            # so transfer these to the new line unchanged.
            for line_i, line in enumerate(infile.readlines()):
                if (line_i == 0) or (line_i == 1) or (line_i == linecount):
                    outfile.writelines(line)
                else:
                    newline = (
                        line[:20]
                        + "".join(["{:8.3f}".format(i) for i in coords[line_i - 2]])
                        + "\n"
                    )
                    outfile.writelines(newline)


def run(filepath, kick_amount, radius):

    logging.info("Kicking file: %s" % filepath)

    if not os.path.isfile(filepath):
        print("Error! Couldn't find file: %s" % filepath)
        sys.exit()

    # Load the coordinate file
    coords, box_vectors = coordload(filepath, return_box_vectors=True)

    # Apply random kicks until atoms are appropriately spaced
    atoms_overlap = True
    cycles = 0

    while atoms_overlap:

        atoms_overlap = False

        # Are there any atoms to kick?
        kickable = get_kickable_atoms_sklearn(coords, radius)

        if len(kickable) == 0:
            # At this point, no atoms in the real box overlap.
            # Let's now consider PBC overlap

            pbc_kickable = get_pbc_overlapping_atoms(coords, box_vectors, radius)

            if len(pbc_kickable) == 0:
                # There's nothing left to kick
                continue
            else:
                atoms_overlap = True
                logging.info("Kicking %s atoms (PBC)" % len(pbc_kickable))
                cycles += 1
                kick(coords, pbc_kickable, kick_amount)

        else:
            atoms_overlap = True
            cycles += 1
            kick(coords, kickable, kick_amount)

    outpath = os.path.join(os.path.dirname(filepath), "kicked.gro")
    write_outfile(filepath, coords, outpath)
    return cycles


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Move overlapping atoms by small amounts to avoid clashes."
    )
    parser.add_argument(
        "filepath", type=str, help="Path to the input file. Must be GRO format."
    )
    parser.add_argument(
        "--kick",
        type=float,
        default=0.05,
        help="Maximum per-axis kick distance (nm) (default: 0.05)",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=0.1,
        help="Neighbour search radius (nm) (default: 0.1)",
    )

    args = parser.parse_args()
    run(args.filepath, args.kick, args.radius)
