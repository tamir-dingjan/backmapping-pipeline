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
    patch_coordinator = PatchCoordinator(config)
    logger.info("Loading patches...")
    patch_coordinator.load_patches()
    patch_coordinator.rebuild_bad_stereo()


if __name__ == "__main__":
    logger.info(f"Starting {__file__}")
    main()
    logger.info(f"Finished {__file__}")
