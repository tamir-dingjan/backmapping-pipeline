import os
import numpy as np
import yaml
import shutil
import glob
import subprocess
from backmapping import kicker
from backmapping.io_utils import load_cg_trajectory, load_coordinates
from backmapping.logger import logger

EXCLUDED_BEADS_FILE = os.path.join(
    os.path.dirname(__file__), "excluded_bead_resnames.yaml"
)

MAX_STEREOCHEM_CORRECTION_ITER = 2


class Patch:

    # TODO: refactor to work with config dict rather than individual settings
    # This will let each patch know where backward.py and the schro.ve files are
    def __init__(
        self, config, patchdir, outname, patch_frame: int, patch_resid: int, traj=None, post_production=False
    ):
        logger.debug("Initializing Patch")
        self.config = config
        self.patchdir = os.path.abspath(patchdir)
        self.outname = outname
        self.outfile = os.path.join(
            self.patchdir,
            self.outname,
        )
        self.cg_traj = config["patches"]["traj"]
        self.cg_top = config["patches"]["top"]
        self.patch_frame = patch_frame
        self.patch_resid = patch_resid
        self.patch_size = config["patches"]["size"]
        self.traj = traj
        self.core = None
        self.excluded_beads = None
        self.lipids = None
        self.stereoconf_file = None
        self.stereochem_correction_iter = 0

        if not post_production:
            # Read in the excluded bead residue names
            # These are used for excluding beads during backmapping
            if not os.path.exists(EXCLUDED_BEADS_FILE):
                logger.error(f"Excluded beads file not found: {EXCLUDED_BEADS_FILE}")
            else:
                try:
                    with open(EXCLUDED_BEADS_FILE, "r") as f:
                        self.excluded_beads = yaml.safe_load(f)
                except yaml.YAMLError as e:
                    logger.error(f"Error loading excluded_beads.yaml: {e}")

            # Check if patchdir exists
            if not os.path.exists(self.patchdir):
                logger.error(f"Patch directory {self.patchdir} does not exist.")
                raise FileNotFoundError(f"Patch directory {self.patchdir} does not exist.")

            # Optionally load the trajectory if passed in
            if self.traj is None:
                logger.debug(
                    f"Loading trajectory {self.cg_traj} with topology {self.cg_top}."
                )
                self.traj = load_cg_trajectory(self.cg_traj, self.cg_top)
            if self.traj is None:
                logger.error(
                    f"Failed to load trajectory {self.cg_traj} with topology {self.cg_top}."
                )
                raise ValueError(
                    f"Failed to load trajectory {self.cg_traj} with topology {self.cg_top}."
                )

            # Validate the patch selection frame and residue are within the trajectory
            if self.patch_frame > len(self.traj):
                logger.error(
                    f"Patch selection frame {self.patch_frame} is out of bounds for trajectory with {len(self.traj)} frames."
                )
                raise ValueError(
                    f"Patch selection frame {self.patch_frame} is out of bounds for trajectory with {len(self.traj)} frames."
                )

            # Exclude the excluded beads from the trajectory
            # This is done to prevent the water residues from clobbering
            # lipid residue numbers
            include = f"not resname {self.excluded_beads["ion"]} and not resname {self.excluded_beads["water"]}"
            lipid_indices = self.traj.topology.select(include)
            self.traj = self.traj.atom_slice(lipid_indices)
            logger.debug(f"Excluded beads from trajectory: {self.excluded_beads}")
            logger.debug(f"Selected {len(lipid_indices)} lipid atoms from trajectory.")

            # Select the patch core residue
            self.core = self.traj.topology.select(f"resSeq {self.patch_resid}")

            if len(self.core) == 0:
                logger.error(
                    f"Patch selection residue {self.patch_resid} not found in trajectory."
                )
                raise ValueError(
                    f"Patch selection residue {self.patch_resid} not found in trajectory."
                )
            logger.debug(f"Selected core residue beads: {self.core}")

    def select_patch_residues(self):
        # TODO: Refactor to work on a single frame (self.patch_frame)
        # Keep a backup of the coordinates
        coord_backup = self.traj.xyz[self.patch_frame].copy()
        box = self.traj.unitcell_lengths[self.patch_frame]
        logger.debug(f"Box dimensions: {box}")

        # Shift the core resid COM to the center of the box
        # NOTE: If the core resid spans a PBC boundary, the COM will be incorrect
        # So, this selection uses a single atom from the core resid

        com = self.traj.xyz[self.patch_frame][self.core[0]]

        # Move the whole bilayer to place the core COM XY at the center of the box
        comshift = np.zeros_like(self.traj.xyz[self.patch_frame])
        comshift[:, 0] = (box[0] / 2) - com[0]
        comshift[:, 1] = (box[1] / 2) - com[1]
        self.traj.xyz[self.patch_frame] += comshift

        # Shift out-of-PBC atoms back into the box
        pbcshift = np.zeros_like(self.traj.xyz[self.patch_frame])
        for axis in [0, 1]:
            pbcshift[
                np.where(self.traj.xyz[self.patch_frame][:, axis] < 0)[0],
                axis,
            ] = box[axis]
            pbcshift[
                np.where(self.traj.xyz[self.patch_frame][:, axis] > box[axis])[0],
                axis,
            ] = (
                -1 * box[axis]
            )

        self.traj.xyz[self.patch_frame] += pbcshift

        # Select the patch residues
        # Start by taking every atom within the XY patch
        # This takes too many lipids. We will filter this to take only lipids whose COM is in the patch zone
        patch_mask = np.zeros_like(self.traj.xyz[self.patch_frame])
        patch_xmax = (
            np.mean(self.traj.xyz[self.patch_frame][self.core], axis=0)[0]
            + 0.5 * self.patch_size
        )
        patch_xmin = (
            np.mean(self.traj.xyz[self.patch_frame][self.core], axis=0)[0]
            - 0.5 * self.patch_size
        )
        patch_ymax = (
            np.mean(self.traj.xyz[self.patch_frame][self.core], axis=0)[1]
            + 0.5 * self.patch_size
        )
        patch_ymin = (
            np.mean(self.traj.xyz[self.patch_frame][self.core], axis=0)[1]
            - 0.5 * self.patch_size
        )

        patch_mask[
            np.where(
                (self.traj.xyz[self.patch_frame][:, 0] > patch_xmin)
                & (self.traj.xyz[self.patch_frame][:, 0] < patch_xmax)
                & (self.traj.xyz[self.patch_frame][:, 1] > patch_ymin)
                & (self.traj.xyz[self.patch_frame][:, 1] < patch_ymax)
            )
        ] = 1

        # Get the atoms indices within the mask
        # We are taking advantage of the fact that the indices of the patch_mask
        # match the indices of the atoms
        patch_atoms = list(
            set(self.traj.topology.atom(x).index for x in np.where(patch_mask == 1)[0])
        )

        # Expand to complete residues
        included_residues = list(
            set(self.traj.topology.atom(x).residue.resSeq for x in patch_atoms)
        )

        # Filter the residues to include by whether the lipid's COM is within the patch boundaries
        # Need to select residues using an ID that isn't shared with ions or water
        com_residues = []
        for res in included_residues:
            com = np.mean(
                self.traj.xyz[self.patch_frame][
                    self.traj.topology.select("resSeq %s" % res)
                ],
                axis=0,
            )
            if (
                (com[0] > patch_xmin)
                & (com[0] < patch_xmax)
                & (com[1] > patch_ymin)
                & (com[1] < patch_ymax)
            ):
                com_residues.append(res)

        patch_whole = self.traj.topology.select(
            "resSeq %s" % " ".join(map(str, com_residues))
        )

        # Extract the bilayer patch coordinates
        self.output = self.traj[self.patch_frame].atom_slice(patch_whole)
        logger.info(
            f"Selected {self.output.n_atoms} beads, {self.output.n_residues} residues."
        )

        # Restore the original trajectory coordinates
        # This is because the COM shift affects the source trajectory,
        # which needs to be reused for the next patch
        self.traj.xyz[self.patch_frame] = coord_backup

    def write_patch(self):
        self.output.save_gro(self.outfile)
        logger.info(f"Saved patch coordinates to {self.outfile}")

    def get_lipid_contents(self):
        # Determine the count and order of lipids to be written to the topology
        lipids = [res.name for res in self.output.topology.residues]

        # Only contiguous lipids of the same residue can be aggregated
        # I.e., if the patch contains lipids [A, B, C, B, B, D]
        # then the counts should be [1,1,1,2,1]
        contents_lipids = []
        contents_counts = []

        prior_lipid = ""
        count = 0

        for lipid in lipids:
            # If the lipid is the same as the previous one,
            # increment contents_counts
            if lipid == prior_lipid:
                count += 1
                contents_counts[-1] = count
            else:
                count = 1
                contents_lipids.append(lipid)
                contents_counts.append(count)
                prior_lipid = lipid

        return lipids, contents_lipids, contents_counts

    def run_backward(self):
        logger.debug(f"Running backward.py for {self.patchdir}")
        output = os.path.join(self.patchdir, "patch_aa.gro")
        args = [
            "python",
            self.config["filepaths"]["backward"],
            "-f",
            self.outfile,
            "-p",
            "patch.top",
            "-o",
            output,
            "-to",
            "charmm36",
            "-from",
            "martini",
            "-kick",
            "0",
        ]
        subprocess.run(args, cwd=self.patchdir)
        if not os.path.isfile(output):
            logger.error(f"Backward.py failed: {self.patchdir}")
            raise FileNotFoundError(f"Backward.py failed: {self.patchdir}")

    def remake_box_vectors(self):
        logger.debug(f"Remaking box vectors for {self.patchdir}")
        output = os.path.join(self.patchdir, "box.gro")
        args = [
            self.config["filepaths"]["gmx"],
            "editconf",
            "-f",
            "patch_aa.gro",
            "-o",
            output,
            "-box",
            str(self.config["patches"]["size"]),
            str(self.config["patches"]["size"]),
            str(self.config["patches"]["size"] * 2),
        ]
        subprocess.run(args, cwd=self.patchdir)
        if not os.path.isfile(output):
            logger.error(f"Box vectors generation failed: {self.patchdir}")
            raise FileNotFoundError(f"Box vectors generation failed: {self.patchdir}")

    def kick_overlapping_atoms(self):
        logger.debug(f"Kicking box file for: {self.patchdir}")

        box_file = os.path.join(self.patchdir, "box.gro")
        output = os.path.join(self.patchdir, "kicked.gro")

        kicker.run(box_file, kick_amount=0.02, radius=0.02)

        if not os.path.isfile(output):
            logger.error(f"Coordinate kicking failed: {self.patchdir}")
            raise FileNotFoundError(f"Coordinate kicking failed: {self.patchdir}")

    def minimise_in_vacuum(self, inname: str, outname: str):
        coord = os.path.join(self.patchdir, f"{inname}.gro")
        output = os.path.join(self.patchdir, f"{outname}.gro")
        if os.path.isfile(coord):
            logger.debug(f"Minimising in vacuum: {coord}")
            args = [
                self.config["filepaths"]["gmx"],
                "grompp",
                "-f",
                os.path.join(self.config["filepaths"]["mdp"]),
                "-c",
                coord,
                "-p",
                "patch.top",
                "-o",
                f"{outname}.tpr",
                "-maxwarn",
                "2",
            ]
            subprocess.run(args, cwd=self.patchdir)
            args = [
                self.config["filepaths"]["gmx"],
                "mdrun",
                "-deffnm",
                f"{outname}",
                "-ntomp",
                "5",
            ]
            subprocess.run(args, cwd=self.patchdir)
        if not os.path.isfile(output):
            logger.error(f"Minimisation in vacuum failed: {self.patchdir}")
            raise FileNotFoundError(f"Minimisation failed: {self.patchdir}")

    def has_correct_stereo(self, coordfile: str):
        # Check if the provided coordinate file has the correct stereoconformation
        coord = os.path.join(self.patchdir, f"{coordfile}.gro")
        output = os.path.join(self.patchdir, "stereo.check")
        execfile = os.path.join(self.patchdir, "run_stereo_check.sh")
        stereo_check = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "stereo/check.py")
        )
        stereo_lookup = os.path.abspath(
            os.path.join(self.config["filepaths"]["stereo"], "stereo.json")
        )
        with open(execfile, "w") as f:
            f.writelines(
                "\n".join(
                    [
                        "#!/usr/bin/bash",
                        # Try to use the same stereo.tpr from the first stereochemistry correction step
                        # This file only stores connectivity, not conformation
                        # f"{self.config["filepaths"]["gmx"]} grompp -f {self.config["filepaths"]["mdp"]} -c {coordfile} -p patch.top -o stereo_check.tpr -maxwarn 1",
                        f"source {self.config["filepaths"]["schro_venv"]}",
                        f"python {stereo_check} {coord} stereo.tpr {stereo_lookup}",
                    ]
                )
            )
        logger.debug(f"Checking stereoconformation for {self.patchdir}")
        subprocess.run(
            ["bash", execfile],
            cwd=self.patchdir,
        )

        if not os.path.isfile(output):
            logger.error(f"Checking stereoconformation failed: {self.patchdir}")
            raise FileNotFoundError(
                f"Checking stereoconformation failed: {self.patchdir}"
            )

        # Load the result of the check and return true or false
        with open(output, "r") as f:
            stereomatch = f.readlines()

        if "True" in stereomatch:
            return True
        else:
            return False

    def apply_correct_stereoisomers(self, coordfile: str):
        # NOTE: Must use a version of Gromacs which produces TPR files
        # readable by the schro.ve virtual environment
        coord = os.path.join(self.patchdir, f"{coordfile}.gro")
        output = os.path.join(self.patchdir, "stereo.gro")
        execfile = os.path.join(self.patchdir, "run_stereo_stamp.sh")
        stereo = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "stereo/stereo.py")
        )
        stereo_lookup = os.path.abspath(
            os.path.join(self.config["filepaths"]["stereo"], "stereo.json")
        )
        with open(execfile, "w") as f:
            f.writelines(
                "\n".join(
                    [
                        "#!/usr/bin/bash",
                        f"{self.config["filepaths"]["gmx"]} grompp -f {self.config["filepaths"]["mdp"]} -c {coord} -p patch.top -o stereo.tpr -maxwarn 1",
                        f"source {self.config["filepaths"]["schro_venv"]}",
                        f"python {stereo} {coord} stereo.tpr {stereo_lookup}",
                    ]
                )
            )
        logger.debug(f"Apply stereoconformation for {self.patchdir}")
        subprocess.run(
            ["bash", execfile],
            cwd=self.patchdir,
        )
        if not os.path.isfile(output):
            logger.error(f"Stereo stamp failed: {self.patchdir}")
            raise FileNotFoundError(f"Stereo stamp failed: {self.patchdir}")

    def restore_residue_numbers(self):
        logger.debug(f"Restoring residue numbers for {self.patchdir}")
        resnums_source = os.path.join(self.patchdir, "kicked.gro")
        resnums_coord = os.path.join(self.patchdir, "stereo.gro")
        resnums_dest = os.path.join(self.patchdir, "resnums.gro")
        resnums_output = []

        with open(resnums_source, "r") as source:
            linecount = sum(1 for line in source) - 1
            source.seek(0)

            with open(resnums_coord, "r") as coord:
                for line_i, (source_line, coord_line) in enumerate(
                    zip(source.readlines(), coord.readlines())
                ):
                    # Skip the first 2 lines and the final line
                    if (line_i == 0) or (line_i == 1) or (line_i == linecount):
                        resnums_output.append(coord_line)
                    else:
                        # Replace the residue number in the coordinate line
                        new_line = source_line[:5] + coord_line[5:]
                        resnums_output.append(new_line)
        with open(resnums_dest, "w") as dest:
            dest.writelines(resnums_output)
        if not os.path.isfile(resnums_dest):
            logger.error(f"Residue number restoration failed: {self.patchdir}")
            raise FileNotFoundError(
                f"Residue number restoration failed: {self.patchdir}"
            )

    def solvate(self):
        logger.debug(f"Solvating patch: {self.patchdir}")
        if self.stereoconf_file == None:
            logger.error("Stereoconformer file not set before solvation!")
            raise ValueError(f"Stereoconformer file not set before solvation")
        elif not os.path.isfile(self.stereoconf_file):
            logger.error(f"Couldn't find stereoconformer file: {self.stereoconf_file}")
            raise FileNotFoundError(
                f"Couldn't find stereoconformer file: {self.stereoconf_file}"
            )
        output = os.path.join(self.patchdir, "solv.gro")
        args = [
            self.config["filepaths"]["gmx"],
            "solvate",
            "-cp",
            f"{self.stereoconf_file}",
            "-p",
            "patch.top",
            "-o",
            output,
        ]
        logger.debug(f"Solvating with command:\n{args}\n")
        subprocess.run(args, cwd=self.patchdir)
        if not os.path.isfile(output):
            logger.error(f"Solvation failed: {self.patchdir}")
            raise FileNotFoundError(f"Solvation failed: {self.patchdir}")

    def delete_membrane_waters(self):
        logger.debug(f"Deleting membrane waters: {self.patchdir}")
        output = os.path.join(self.patchdir, "solv_fix.gro")
        args = [
            "perl",
            self.config["filepaths"]["water_deletor"],
            "-in",
            os.path.join(self.patchdir, "solv.gro"),
            "-out",
            output,
            "-ref",
            "O31",
            "-middle",
            "C216",
            "-nwater",
            "3",
        ]
        water_deletion = subprocess.run(
            args, cwd=self.patchdir, capture_output=True, text=True
        )
        remaining_waters = water_deletion.stdout.split("\n")[-6].split()[0]

        patch_topol = os.path.join(self.patchdir, "patch.top")
        system_topol = os.path.join(self.patchdir, "system.top")
        with open(patch_topol, "r") as source:
            with open(system_topol, "w") as dest:
                for line in source.readlines():
                    if "SOL" in line:
                        newline = line.split()[0] + "\t" + remaining_waters + "\n"
                    else:
                        newline = line
                    dest.writelines(newline)
                dest.write("\n")

        if not (os.path.isfile(system_topol) and os.path.isfile(output)):
            logger.error(f"Water deletion failed: {self.patchdir}")
            raise FileNotFoundError(f"Water deletion failed: {self.patchdir}")
        logger.info(f"New water count: {remaining_waters}")

    def add_ions(self):
        output = os.path.join(self.patchdir, "ions.gro")
        args = [
            self.config["filepaths"]["gmx"],
            "grompp",
            "-f",
            os.path.join(self.config["filepaths"]["mdp"]),
            "-c",
            os.path.join(self.patchdir, "solv_fix.gro"),
            "-p",
            os.path.join(self.patchdir, "system.top"),
            "-o",
            os.path.join(self.patchdir, "ions.tpr"),
            "-maxwarn",
            "1",
        ]
        subprocess.run(args, cwd=self.patchdir)
        if not os.path.isfile(os.path.join(self.patchdir, "ions.tpr")):
            logger.error(f"Ion generation failed: {self.patchdir}")
            raise FileNotFoundError(f"Ion generation failed: {self.patchdir}")

        execfile = os.path.join(self.patchdir, "make_ions.sh")
        with open(execfile, "w") as f:
            f.writelines(
                "\n".join(
                    [
                        "#!/usr/bin/bash",
                        f"{self.config["filepaths"]["gmx"]} genion -s ions.tpr -p system.top -o ions.gro -neutral -conc 0.15 << EOF",
                        "SOL",
                        "EOF",
                    ]
                )
            )
        subprocess.run(
            ["bash", execfile],
            cwd=self.patchdir,
        )
        if not os.path.isfile(output):
            logger.error(f"Ion insertion failed: {self.patchdir}")
            raise FileNotFoundError(f"Ion insertion failed: {self.patchdir}")

    def generate_index(self):
        output = os.path.join(self.patchdir, "index.ndx")

        # Make initial index file
        execfile = os.path.join(self.patchdir, "make_index.sh")
        with open(execfile, "w") as f:
            f.writelines(
                "\n".join(
                    [
                        "#!/usr/bin/bash",
                        f"{self.config["filepaths"]["gmx"]} make_ndx -f ions.gro -o index.ndx << EOF",
                        "q",
                        "EOF",
                    ]
                )
            )
        idx_groups = subprocess.run(
            ["bash", execfile], cwd=self.patchdir, capture_output=True, text=True
        )
        if not os.path.isfile(output):
            logger.error(f"Index generation failed: {self.patchdir}")
            raise FileNotFoundError(f"Index generation failed: {self.patchdir}")

        # Merge water and ions
        lipid_group = -1
        for line in idx_groups.stdout.split("\n"):
            if "Water_and_ions" in line:
                lipid_group = str(int(line.split()[0]) + 1)
                break
        if lipid_group == -1:
            logger.error(
                f"Failed to find water_and_ions group in index file: {self.patchdir}"
            )
            raise ValueError(
                f"Failed to find water_and_ions group in index file: {self.patchdir}"
            )

        # Make lipid group
        execfile = os.path.join(self.patchdir, "edit_index.sh")
        with open(execfile, "w") as f:
            f.writelines(
                "\n".join(
                    [
                        "#!/usr/bin/bash",
                        f"{self.config["filepaths"]["gmx"]} make_ndx -f ions.gro -n index.ndx -o index.ndx << EOF",
                        '!"Water_and_ions"',
                        f"name {lipid_group} lipid",
                        "q",
                        "EOF",
                    ]
                )
            )
        subprocess.run(
            ["bash", execfile],
            cwd=self.patchdir,
        )

    def minimise(self):
        logger.debug(f"Minimising patch: {self.patchdir}")
        output = os.path.join(self.patchdir, "min.tpr")
        args = [
            self.config["filepaths"]["gmx"],
            "grompp",
            "-f",
            os.path.join(self.config["filepaths"]["mdp"]),
            "-c",
            os.path.join(self.patchdir, "ions.gro"),
            "-p",
            os.path.join(self.patchdir, "system.top"),
            "-o",
            os.path.join(self.patchdir, "min.tpr"),
            "-n",
            os.path.join(self.patchdir, "index.ndx"),
            "-maxwarn",
            "0",
        ]
        subprocess.run(args, cwd=self.patchdir)
        if not os.path.isfile(output):
            logger.error(f"Minimisation preparation failed: {self.patchdir}")
            raise FileNotFoundError(f"Minimisation preparation failed: {self.patchdir}")

        args = [
            self.config["filepaths"]["gmx"],
            "mdrun",
            "-deffnm",
            "min",
            "-ntomp",
            "5",
        ]
        subprocess.run(args, cwd=self.patchdir)
        if not os.path.isfile(os.path.join(self.patchdir, "min.gro")):
            logger.error(f"Minimisation failed: {self.patchdir}")
            raise FileNotFoundError(f"Minimisation failed: {self.patchdir}")

    def extract_minimised_lipids(self, coordfile: str, outfile: str):
        coord = os.path.join(self.patchdir, f"{coordfile}.gro")
        topol = os.path.join(self.patchdir, f"{coordfile}.tpr")
        out = os.path.join(self.patchdir, f"{outfile}.gro")
        
        if not (os.path.exists(coord) and os.path.exists(topol)):
            raise FileNotFoundError(f"Could not find coordinate file ({coordfile}.gro) or topology file ({coordfile}.tpr)")

        u = load_coordinates(coord, topol)
        u.select_atoms("not resname CL NA SOL").write(out)
        logger.debug(f"Extracted minimised lipids to file: {out}")

    def desolvate_topology(self):
        dry_topology = []
        with open(os.path.join(self.patchdir, "patch.top"), "r") as infile:
            for line in infile.readlines():
                if not "SOL" in line:
                    dry_topology.append(line)
        with open(os.path.join(self.patchdir, "patch.top"), "w") as outfile:
            outfile.writelines(dry_topology)


class PatchCoordinator:

    def __init__(self, config):
        self.config = config
        self.patches = []

        # Process the patch selection settings to determine how many patches to make
        self.frame_range = range(
            self.config["patches"]["frame_start"],
            self.config["patches"]["frame_end"],
            self.config["patches"]["frame_stride"],
        )
        self.resid_range = self.config["patches"]["resid"]

    def copy_simulation_parameters(self):
        # Set up the patch forcefield, lipid topology/parameter files, and minimisation settings
        self.topol = os.path.abspath(
            os.path.join(self.config["patches"]["patchdir"], "topol")
        )
        logger.info(f"Copying parameters to {self.topol}")
        if not os.path.exists(self.topol):
            os.makedirs(self.topol)

        # Charmm forcefield
        self.charmm36_ff = os.path.join(
            self.topol, os.path.basename(self.config["filepaths"]["charmm36_ff"])
        )
        shutil.copytree(
            self.config["filepaths"]["charmm36_ff"],
            self.charmm36_ff,
            dirs_exist_ok=True,
        )

        # Additional cgenff params
        self.additional_params = os.path.join(self.topol, "additional_params.prm")
        shutil.copy2(
            os.path.join(
                self.config["filepaths"]["aa_gmx_parameters"], "additional_params.prm"
            ),
            self.additional_params,
        )

        # Lipid topology files
        lipid_topology_files = glob.glob(
            os.path.join(self.config["filepaths"]["aa_gmx_parameters"], "*.itp")
        )
        for file in lipid_topology_files:
            shutil.copy2(file, os.path.join(self.topol, os.path.basename(file)))

        # Minimisation settings
        shutil.copy2(
            self.config["filepaths"]["mdp"],
            self.topol,
        )

    def make_patches(self):
        traj = load_cg_trajectory(
            self.config["patches"]["traj"], self.config["patches"]["top"]
        )
        for frame in self.frame_range:
            for resid in self.resid_range:
                # Make a patch directory if it doesn't exist
                patchdir = os.path.join(
                    self.config["patches"]["patchdir"], f"patch_{frame}_{resid}"
                )
                if not os.path.exists(patchdir):
                    os.makedirs(patchdir)

                patch = Patch(
                    config=self.config,
                    patchdir=patchdir,
                    outname=f"patch_{frame}_{resid}.gro",
                    patch_frame=frame,
                    patch_resid=resid,
                    traj=traj,
                )

                patch.select_patch_residues()
                patch.write_patch()

                self.patches.append(patch)

    def make_patch_topologies(self):
        for patch in self.patches:
            lipids, contents_lipids, contents_counts = patch.get_lipid_contents()

            topology_file = os.path.join(
                patch.patchdir,
                "patch.top",
            )

            with open(topology_file, "w") as f:
                f.write("#include ")
                f.write(f'"{self.charmm36_ff}/forcefield.itp"')
                f.write("\n")

                f.write("#include ")
                f.write(f'"{self.additional_params}"')
                f.write("\n")

                seen = set()
                for lipid in lipids:
                    if not lipid in seen:
                        f.write("#include ")
                        f.write('"' + os.path.join(self.topol, f"{lipid}.itp") + '"')
                        f.write("\n")
                        seen.add(lipid)
                f.write("\n")

                # Water and ions
                f.write("; Include water topology\n")
                f.write("#include ")
                f.write(f'"{self.charmm36_ff}/tip3p.itp"\n')
                f.write("#ifdef POSRES_WATER\n")
                f.write("; Position restraint for each water oxygen \n")
                f.write("[ position_restraints ] \n")
                f.write("; i funct       fcx        fcy        fcz \n")
                f.write("1    1       1000       1000       1000 \n")
                f.write("#endif \n")
                f.write("; Include topology for ions \n")
                f.write("#include ")
                f.write(f'"{self.charmm36_ff}/ions.itp"\n')

                # Lipid contents
                f.write("[ system ]\n")
                f.write(f"{os.path.basename(patch.patchdir)}\n")
                f.write("[ molecules ]\n")

                f.writelines(
                    [
                        lipid + "\t" + str(count) + "\n"
                        for (lipid, count) in zip(contents_lipids, contents_counts)
                    ]
                )

                patch.lipids = list(seen)

    def backmap_all_patches(self):
        for patch in self.patches:
            patch.run_backward()
            patch.remake_box_vectors()
            patch.kick_overlapping_atoms()
            patch.minimise_in_vacuum("kicked", "minvac")
            # The residue numbers are reset to begin from
            # 1 in minvac.gro. This should not be a problem,
            # but they can be copied back over from kicked.gro:
            # patch.restore_residue_numbers()

    def correct_patch_stereoconformation_single(self, patch, coordfile: str):
        # This method is in the PatchCoordinator rather than the Patch class
        # because it tells the Patch what to do.
        patch.apply_correct_stereoisomers(coordfile)
        # Produces "stereo.gro"

        patch.minimise_in_vacuum(
            "stereo", f"stereo_min_{patch.stereochem_correction_iter}"
        )
        # Produces "stereo_min_0.gro"

        while not (
            patch.has_correct_stereo(f"stereo_min_{patch.stereochem_correction_iter}")
            or patch.stereochem_correction_iter > MAX_STEREOCHEM_CORRECTION_ITER
        ):
            patch.apply_correct_stereoisomers(
                f"stereo_min_{patch.stereochem_correction_iter}"
            )
            # Overwrites "stereo.gro"

            patch.stereochem_correction_iter += 1
            patch.minimise_in_vacuum(
                "stereo", f"stereo_min_{patch.stereochem_correction_iter}"
            )
        patch.stereoconf_file = os.path.abspath(
            os.path.join(
                patch.patchdir, f"stereo_min_{patch.stereochem_correction_iter}.gro"
            )
        )

    def correct_patch_stereoconformation(self):
        # Continue attempting to correct stereochemistry until MAX_STEREOCHEM_CORRECTION_ITER
        # is reached
        try:
            for patch in self.patches:
                self.correct_patch_stereoconformation_single(patch)

        except Exception as e:
            logger.error(e)

    def prepare_patches_for_simulation(self):
        try:
            for patch in self.patches:
                patch.solvate()
                patch.delete_membrane_waters()
                patch.add_ions()
                patch.generate_index()
                patch.minimise()
        except Exception as e:
            logger.error(e)

    def load_patches(self, post_production=False):
        traj = load_cg_trajectory(
            self.config["patches"]["traj"], self.config["patches"]["top"]
        )
        for frame in self.frame_range:
            for resid in self.resid_range:

                patchdir = os.path.join(
                    self.config["patches"]["patchdir"], f"patch_{frame}_{resid}"
                )

                patch = Patch(
                    config=self.config,
                    patchdir=patchdir,
                    outname=f"patch_{frame}_{resid}.gro",
                    patch_frame=frame,
                    patch_resid=resid,
                    traj=traj,
                    post_production=post_production
                )

                self.patches.append(patch)

    def process_per_patch(self):
        for patch in self.patches:
            logger.info(f"Processing patch: {patch.patchdir}")
            try:
                patch.run_backward()
                patch.remake_box_vectors()
                patch.kick_overlapping_atoms()
                patch.minimise_in_vacuum("kicked", "minvac")
                self.correct_patch_stereoconformation_single(patch, "minvac")
                patch.solvate()
                patch.delete_membrane_waters()
                patch.add_ions()
                patch.generate_index()
                patch.minimise()
                # Re-check correct stereo after minimisation in solvent
                # this requires first preparing a lipid-only coordinate file from the minimised system
                patch.extract_minimised_lipids("min", "min_lipids")
                iter_min_solvent = 0
                while not (
                    patch.has_correct_stereo("min_lipids")
                    or (iter_min_solvent > MAX_STEREOCHEM_CORRECTION_ITER)
                ):
                    logger.info(f"Rebuilding patch")
                    # If incorrect stereo, correct the patch stereo from the minimised system
                    # Iteratively correct stereoconformation of extracted lipids
                    # Reset the correction iterator first to force overwriting
                    # the previous vacumm-minimised coordinate files
                    patch.stereochem_correction_iter = 0
                    # Desolvate the patch.top to allow minimisation in vacuum
                    patch.desolvate_topology()
                    self.correct_patch_stereoconformation_single(patch, "min_lipids")
                    # This will update the patch.stereoconf_file, so we can redo the following
                    # steps with new coordinates
                    patch.solvate()
                    patch.delete_membrane_waters()
                    patch.add_ions()
                    patch.generate_index()
                    patch.minimise()
                    patch.extract_minimised_lipids("min", "min_lipids")
                    iter_min_solvent += 1

            except Exception as e:
                logger.error(e)

    def check_stereo(self):
        for patch in self.patches:
            logger.info(f"Checking stereo for patch: {patch.patchdir}")
            try:
                patch.extract_minimised_lipids("min", "min_lipids")
                result = patch.has_correct_stereo("min_lipids")
                logger.info(f"Stereocheck result: {result}")
            except Exception as e:
                logger.error(e)
                
    def check_stereo_post_production(self):
        for patch in self.patches:
            logger.info(f"Checking stereo for md0.gro file in patch: {patch.patchdir}")
            try:
                patch.extract_minimised_lipids("md0", "md0_lipids")
                result = patch.has_correct_stereo("md0_lipids")
                logger.info(f"Stereocheck result: {result}")
            except Exception as e:
                logger.error(e)

    def rebuild_bad_stereo(self):
        for patch in self.patches:
            logger.info(f"Checking stereo for patch: {patch.patchdir}")
            try:
                patch.extract_minimised_lipids("min", "min_lipids")
                iter_min_solvent = 0
                while not (
                    patch.has_correct_stereo("min_lipids")
                    or (iter_min_solvent > MAX_STEREOCHEM_CORRECTION_ITER)
                ):
                    logger.info(f"Rebuilding patch")
                    patch.run_backward()
                    patch.remake_box_vectors()
                    patch.kick_overlapping_atoms()
                    # The patch.top file contains a solvent entry, so remove this before attempting to minimise in vacuum
                    patch.desolvate_topology()
                    patch.minimise_in_vacuum("kicked", "minvac")
                    # Reset the correction iterator first to force overwriting
                    # the previous vacumm-minimised coordinate files
                    patch.stereochem_correction_iter = 0
                    self.correct_patch_stereoconformation_single(patch, "minvac")
                    patch.solvate()
                    patch.delete_membrane_waters()
                    patch.add_ions()
                    patch.generate_index()
                    patch.minimise()
                    patch.extract_minimised_lipids("min", "min_lipids")
                    iter_min_solvent += 1

            except Exception as e:
                logger.error(e)
