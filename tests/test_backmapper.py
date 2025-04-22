import pytest
import os
import networkx as nx
from backmapping.workflow import BackmappingWorkflow
from backmapping.io_utils import load_itp_as_network
from backmapping.lipidclass import LipidClass


@pytest.fixture
def backmapper():
    # Define the paths to the input files
    aa_coords = os.path.join(os.path.dirname(__file__), "data", "aa_gmx", "PSM_box.pdb")
    basename = "PSM"
    aa_param = os.path.join(os.path.dirname(__file__), "data", "aa_gmx", "PSM.itp")
    cg_param = os.path.join(os.path.dirname(__file__), "data", "cg_martini", "PSM.itp")
    stereo = os.path.join(os.path.dirname(__file__), "data", "stereo", "stereo.json")
    mapping_file = os.path.join(
        os.path.dirname(__file__), "data", "Mapping", "PSM.charmm36.map"
    )

    # Create a BackmappingWorkflow instance
    backmapper = BackmappingWorkflow(
        aa_coords, basename, aa_param, cg_param, stereo, mapping_file
    )
    return backmapper


def test_loading_input_files(backmapper):
    """
    Test loading the input files for the backmapping workflow.
    """

    aa_param = "tests/data/aa_gmx/PSM.itp"
    cg_param = "tests/data/cg_martini/PSM.itp"

    # Check if the input files are loaded correctly
    assert nx.is_isomorphic(backmapper.aa, load_itp_as_network(aa_param))
    assert nx.is_isomorphic(backmapper.cg, load_itp_as_network(cg_param))
    assert backmapper.cg_bead_order == "NC3 PO4 AM1 AM2 T1A C2A C3A C1B C2B C3B C4B \n"


def test_lipid_class_detection(backmapper):
    backmapper.set_lipid_class()
    assert backmapper.lipid_class == LipidClass.SM


def test_region_assignment(backmapper):
    backmapper.set_lipid_class()
    backmapper.assign_regions()
    assert backmapper.cg_regions == {
        "head": ["1", "2"],
        "linker": ["3", "4"],
        "tailA": ["5", "6", "7"],
        "tailB": ["8", "9", "10", "11"],
    }


def test_allocate_head(backmapper):
    backmapper.set_lipid_class()
    backmapper.assign_regions()
    backmapper.allocate_head()
    assert backmapper.aa.nodes["1"]["map"] == ["1"]


def test_allocate_linker(backmapper):
    backmapper.set_lipid_class()
    backmapper.assign_regions()
    backmapper.allocate_linker()
    # C1S is mapped partway between PO4 and AM1
    assert backmapper.aa.nodes["25"]["map"] == ["3", "2"]
    # C2S is mapped partway between AM1 and AM2
    assert backmapper.aa.nodes["29"]["map"] == ["3", "4"]


def test_allocate_tail(backmapper):
    backmapper.set_lipid_class()
    backmapper.assign_regions()
    backmapper.allocate_both_tails()
    # C18S is mapped to the final TailA bead
    assert backmapper.aa.nodes["75"]["map"] == ["7"]
    # C216 is mapped to the final TailB bead
    assert backmapper.aa.nodes["124"]["map"] == ["11"]


def test_mapping_file_output(backmapper):
    backmapper.set_lipid_class()
    backmapper.assign_regions()
    backmapper.allocate_head()
    backmapper.allocate_linker()
    backmapper.allocate_both_tails()
    backmapper.catch_unmapped_atoms()
    backmapper.write_mapping_file()

    assert os.path.isfile(
        os.path.join(os.path.dirname(__file__), "data", "Mapping", "PSM.charmm36.map")
    )
