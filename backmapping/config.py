import yaml
import os
from backmapping.logger import logger

logger.debug("Reading config file")


def load_config(config_path="config.yaml"):
    try:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
            logger.setLevel(config["logging"]["level"])
            absolutize_paths(config, os.path.dirname(config_path))
        return config

    except Exception as e:
        logger.error(f"Error loading config file at {config_path}: {e}")


def check_fields(config, fields: list):
    logger.debug("Checking config fields")
    for field in fields:
        if not field in config.keys():
            logger.warning(f"Warning! Couldn't find field {field} in config file.")
        else:
            logger.debug(f"Found field {field}")


def convert_relative_to_absolute_path(path: str, root: str):
    if not os.path.exists(root):
        logger.error(f"Error! Couldn't find project root path: {root}")
    elif not os.path.exists(os.path.join(root, path)):
        logger.error(f"Error! Couldn't find path: {os.path.join(root,path)}")
    else:
        return os.path.abspath(os.path.join(root, path))


def absolutize_paths(config: dict, root: str):
    if not "filepaths" in config.keys():
        logger.error(f"Error! No filepaths defined in config file.")
    else:
        for label, path in config["filepaths"].items():
            if not os.path.isabs(path):
                config["filepaths"][label] = convert_relative_to_absolute_path(
                    path, root
                )
