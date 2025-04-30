import os
import numpy as np
from backmapping.patch import PatchCoordinator

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


def get_patch_coordinator():
    """
    Create a PatchCoordinator instance with the specified configuration.
    """
    # Create a PatchCoordinator instance
    patch_coordinator = PatchCoordinator(config_data)
    return patch_coordinator


def test_patch_coordinator_initialization():
    patch_coordinator = get_patch_coordinator()

    # Check if the patch coordinator is initialized correctly
    # TODO: correct these assertions
    assert patch_coordinator.patchdir == config_data["patches"]["patchdir"]
    assert patch_coordinator.outname == "patch_0_1.pdb"
    assert patch_coordinator.cg_traj == config_data["patches"]["traj"]
    assert patch_coordinator.cg_top == config_data["patches"]["top"]
    assert patch_coordinator.patch_frame == 0
    assert patch_coordinator.patch_resid == 1
    assert patch_coordinator.patch_size == config_data["patches"]["size"]


if __name__ == "__main__":
    # Run the test functions if this script is executed directly
    test_patch_coordinator_initialization()
