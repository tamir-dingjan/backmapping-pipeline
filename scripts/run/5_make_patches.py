import argparse
from backmapping.config import load_config, check_fields
from backmapping.logger import logger
from backmapping.patch import PatchCoordinator


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
        "patchdir",
        "traj",
        "top",
        "frame_start",
        "frame_end",
        "frame_stride",
        "resid",
        "size",
    ]
    check_fields(config["patches"], required_paths)

    # Instantiate the PatchCoordinator
    logger.info(f"Making patches from trajectory: {config['patches']['traj']}")
    # logger.warning("This is a test")

    patch_coordinator = PatchCoordinator(config)

    patch_coordinator.copy_simulation_parameters()
    patch_coordinator.make_patches()
    patch_coordinator.make_patch_topologies()
    # patch_coordinator.backmap_all_patches()
    # patch_coordinator.correct_patch_stereoconformation()
    # patch_coordinator.prepare_patches_for_simulation()

    patch_coordinator.process_per_patch()


if __name__ == "__main__":
    logger.info(f"Starting {__file__}")
    main()
    logger.info(f"Finished {__file__}")
