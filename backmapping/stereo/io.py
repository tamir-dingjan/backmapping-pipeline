import MDAnalysis as mda
import MDAnalysis.transformations as trans
import json


def load_coordinates(coord, topol):
    u = mda.Universe(topol, coord)
    transforms = [trans.unwrap(u.select_atoms("all"))]
    u.trajectory.add_transformations(*transforms)
    return u


def load_stereo_lookup(path):
    with open(path, "r") as f:
        stereo_lookup = json.load(f)
    return stereo_lookup
