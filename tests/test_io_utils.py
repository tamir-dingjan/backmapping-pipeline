import os
import networkx as nx
from backmapping.io_utils import load_itp_as_network, load_structure_from_file


def test_load_itp_file_as_network():
    itp_path = os.path.join(os.path.dirname(__file__), "data", "aa_gmx", "PSM.itp")
    network = load_itp_as_network(itp_path)

    assert network != None
    assert network.number_of_nodes() == 127


def test_load_structure_from_file():
    coord_path = os.path.join(
        os.path.dirname(__file__), "data", "aa_gmx", "PSM_box.pdb"
    )
    structure = load_structure_from_file(coord_path)

    assert structure != None
