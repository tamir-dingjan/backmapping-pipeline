import os
import tempfile
import yaml
from backmapping.config import load_config


def test_load_config_converts_relative_paths(tmp_path):
    config_data = {
        "logging": {"level": "DEBUG"},
        "filepaths": {
            "relative_path": "a/relative/path",
            "absolute_path": "/absolute/path",
        },
    }

    config_file = os.path.join(tmp_path, "test_config.yaml")

    # Create the relative path
    os.makedirs(os.path.join(tmp_path, "a/relative/path"))

    with open(config_file, "w") as f:
        yaml.dump(config_data, f)

    config = load_config(config_file)

    # relative path should now be absolute
    assert os.path.isabs(config["filepaths"]["relative_path"])
    assert config["filepaths"]["relative_path"] == os.path.abspath(
        os.path.join(tmp_path, "a/relative/path")
    )

    # absolute path should be unchanged
    assert config["filepaths"]["absolute_path"] == "/absolute/path"
