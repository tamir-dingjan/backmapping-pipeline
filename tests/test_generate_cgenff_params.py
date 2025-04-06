import os
import shutil
import yaml
import subprocess
from backmapping.config import load_config

# Define config entries for this test
config_data = {
    "logging": {"level": "DEBUG"},
    "filepaths": {
        "aa_coords": "data/aa_coords",
        "aa_cgenff_parameters": "data/aa_cgenff",
        "cgenff": "/opt/silcsbio.2024.1/cgenff/cgenff",
    },
}


def count_atom_lines(path):
    count = 0
    with open(path, "r") as f:
        for line in f.readlines():
            if "ATOM" in line:
                count += 1
    return count


def test_generate_cgenff_params_produces_output_files(tmp_path):
    # Populate temp dir with data
    os.makedirs(os.path.join(tmp_path, "data", "aa_coords"))
    os.makedirs(os.path.join(tmp_path, "data", "aa_cgenff"))

    shutil.copy(
        os.path.join(os.path.dirname(__file__), "data", "aa_coords", "PSM.mol2"),
        os.path.join(tmp_path, "data", "aa_coords"),
    )

    config_path = os.path.join(tmp_path, "config.yaml")
    with open(config_path, "w") as f:
        yaml.dump(config_data, f)

    # Preparation script path
    script_path = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "scripts",
        "prepare",
        "0_generate_cgenff_params.py",
    )

    # Run script
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        f"{os.path.dirname(os.path.dirname(__file__))}:{env['PYTHONPATH']}"
    )
    result = subprocess.run(
        args=["python", script_path, "--config", config_path],
        cwd=tmp_path,
        env=env,
        capture_output=True,
    )

    # Validate
    assert result.returncode == 0, f"Script failed with return code: {result.stderr}"
    assert os.path.isfile(os.path.join(tmp_path, "data", "aa_cgenff", "PSM.str"))
    assert count_atom_lines(
        os.path.join(os.path.dirname(__file__), "data", "aa_cgenff", "PSM.str")
    ) == count_atom_lines(os.path.join(tmp_path, "data", "aa_cgenff", "PSM.str"))
