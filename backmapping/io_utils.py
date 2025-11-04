from rdkit import Chem
import json
import networkx as nx
import fnmatch
import mdtraj as md
import MDAnalysis as mda
import MDAnalysis.transformations as trans
import os
from backmapping.logger import logger


def load_structure_from_pdb_file(path: str):
    structure = Chem.rdmolfiles.MolFromPDBFile(path, sanitize=True, removeHs=False)
    if structure is None:
        msg = f"Problem loading structure from file: {path}"
        logger.error(msg)
        raise Exception(msg)
    return structure


def load_structure_from_mol2_file(path: str):
    structure = Chem.rdmolfiles.MolFromMol2File(path, sanitize=True, removeHs=False)
    if structure is None:
        msg = f"Problem loading structure from file: {path}"
        logger.error(msg)
        raise Exception(msg)
    return structure


def load_itp_as_network(filepath):
    with open(filepath, "r") as itp:
        # Construct a new network
        topology = nx.Graph()

        current_section = None
        heteroflag = None

        for line in itp.readlines():
            if (line[0] == ";") or (line[0] == "#"):
                continue

            elif len(line.split()) == 0:
                # Blank lines mark the end of a section
                current_section = None

            elif current_section == "atoms":
                # Add this atom as a network node
                atom = line.split()
                # Define heteroatom flag based on the first character of the atomname
                if (atom[1][0] == "H") or (atom[1][0] == "C"):
                    heteroflag = False
                else:
                    heteroflag = True

                topology.add_node(
                    atom[0],
                    element=atom[1][0],
                    atomtype=atom[1],
                    name=atom[4],
                    hetero=heteroflag,
                )

            elif current_section == "bonds":
                # Add this bond as an edge
                bond = line.split()
                topology.add_edge(bond[0], bond[1])

            elif fnmatch.fnmatch(line, "*[[]*atoms*[]]*"):
                current_section = "atoms"

            elif fnmatch.fnmatch(line, "*[[]*bonds*[]]*"):
                current_section = "bonds"

    return topology


def get_bead_order_from_itp(filepath):
    with open(filepath, "r") as itp:
        for line in itp.readlines():
            if fnmatch.fnmatch(line, ";@BEADS*"):
                return line.split("BEADS ")[-1]


def load_stereo_reference(path: str):
    with open(path, "r") as f:
        stereo_reference = json.load(f)
    return stereo_reference


def load_cg_trajectory(traj, top):
    if not os.path.isfile(traj):
        msg = f"Trajectory file {traj} does not exist."
        logger.error(msg)
        raise FileNotFoundError(msg)
    if not os.path.isfile(top):
        msg = f"Topology file {top} does not exist."
        logger.error(msg)
        raise FileNotFoundError(msg)
    try:
        return md.load(traj, top=top)
    except Exception as e:
        msg = f"Error loading trajectory {traj} and topology {top}: {e}"
        logger.error(msg)
        raise


def load_coordinates(coord, topol):
    u = mda.Universe(topol, coord)
    transforms = [trans.unwrap(u.select_atoms("all"))]
    u.trajectory.add_transformations(*transforms)
    return u
