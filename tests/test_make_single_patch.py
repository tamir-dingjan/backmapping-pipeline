import os
import numpy as np
import mdtraj as md
from backmapping.patch import Patch

# Define patch selection config entries for this test
config_data = {
    "patches": {
        "patchdir": "data/patches",
        "traj": "data/cg_traj/cg.xtc",
        "top": "data/cg_traj/cg.gro",
        "frame_start": 0,
        "frame_end": 2,
        "frame_stride": 1,
        "resid": [1, 2],
        "size": 5,
    }
}


# TODO: rewrite to use fixture
def single_patch():
    patch_instance = Patch(
        patchdir=config_data["patches"]["patchdir"],
        outname="test_patch.gro",
        cg_traj=config_data["patches"]["traj"],
        cg_top=config_data["patches"]["top"],
        patch_frame=0,
        patch_resid=1,
        patch_size=config_data["patches"]["size"],
    )
    return patch_instance


def test_patch_instance_configs_are_passed_correctly(single_patch):

    assert single_patch.patchdir == config_data["patches"]["patchdir"]
    assert single_patch.outname == "test_patch.gro"
    assert single_patch.cg_traj == config_data["patches"]["traj"]
    assert single_patch.cg_top == config_data["patches"]["top"]
    assert single_patch.patch_frame == 0
    assert single_patch.patch_resid == 1
    assert single_patch.patch_size == config_data["patches"]["size"]
    assert single_patch.traj is not None
    np.testing.assert_equal(
        single_patch.core, np.asarray([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    )
    assert single_patch.excluded_beads is not None


def test_patch_lipid_selection(single_patch):
    assert single_patch.output is not None
    assert single_patch.output.n_frames == 1
    assert single_patch.output.n_atoms == 1035
    assert single_patch.output.n_residues == 90
    assert single_patch.output.topology.atom(1034).name == "C6B"
    assert single_patch.output.topology.atom(1034).residue.name == "KfS2"
    assert single_patch.output.topology.atom(1034).residue.resSeq == 4686
    assert single_patch.output.topology.atom(1034).residue.n_atoms == 13


def test_patch_saved_to_disk_is_correct(single_patch):
    # Check if the patch file is created
    disk_patch_path = os.path.join(
        os.path.dirname(__file__),
        config_data["patches"]["patchdir"],
        single_patch.outname,
    )
    assert os.path.isfile(disk_patch_path)

    # Load the saved patch file and check its contents
    disk_patch = md.load(disk_patch_path, top=disk_patch_path)
    assert disk_patch is not None
    assert disk_patch.n_frames == 1
    assert disk_patch.n_atoms == 1035
    assert disk_patch.n_residues == 90
    assert disk_patch.topology.atom(1034).name == "C6B"
    assert disk_patch.topology.atom(1034).residue.name == "KfS2"
    assert disk_patch.topology.atom(1034).residue.resSeq == 4686
    assert disk_patch.topology.atom(1034).residue.n_atoms == 13


if __name__ == "__main__":
    patch = single_patch()
    test_patch_instance_configs_are_passed_correctly(patch)

    patch.select_patch_residues()
    test_patch_lipid_selection(patch)

    patch.write_patch()
    test_patch_saved_to_disk_is_correct(patch)
