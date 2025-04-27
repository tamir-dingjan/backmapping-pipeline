import os
import numpy as np
import yaml
import shutil
import glob
import subprocess
from backmapping.io_utils import load_cg_trajectory
from backmapping.logger import logger

EXCLUDED_BEADS_FILE = os.path.join(
    os.path.dirname(__file__), "excluded_bead_resnames.yaml"
)


class Patch:

    # TODO: refactor to work with config dict rather than individual settings
    # This will let each patch know where backward.py and the schro.ve files are
    def __init__(
        self,
        patchdir,
        outname,
        cg_traj,
        cg_top,
        patch_frame: int,
        patch_resid: int,
        patch_size: float,
    ):
        logger.debug("Initializing Patch")
        self.patchdir = patchdir
        self.outname = outname
        self.cg_traj = cg_traj
        self.cg_top = cg_top
        self.patch_frame = patch_frame
        self.patch_resid = patch_resid
        self.patch_size = patch_size
        self.traj = None
        self.core = None
        self.excluded_beads = None
        self.lipids = None

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

        self.traj = load_cg_trajectory(self.cg_traj, self.cg_top)
        if self.traj is None:
            logger.error(
                f"Failed to load trajectory {self.cg_traj} with topology {self.cg_top}."
            )
            raise ValueError(
                f"Failed to load trajectory {self.cg_traj} with topology {self.cg_top}."
            )
        logger.debug(f"Loaded trajectory {self.cg_traj} with topology {self.cg_top}.")

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

        # Save out the bilayer patch coordinates
        self.output = self.traj[self.patch_frame].atom_slice(patch_whole)
        logger.debug(
            f"Selected {self.output.n_atoms} beads, {self.output.n_residues} residues."
        )

    def write_patch(self):
        self.output.save_gro(
            os.path.join(
                self.patchdir,
                self.outname,
            )
        )
        logger.debug(f"Saved patch coordinates to {self.patchdir}/{self.outname}")

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
        pass

    def remake_box_vectors(self):
        pass

    def kick_overlapping_atoms(self):
        pass

    def minimise_in_vacuum(self):
        pass

    def apply_correct_stereoisomers(self):
        pass

    def restore_residue_numbers(self):
        pass

    def solvate(self):
        pass

    def delete_membrane_waters(self):
        pass

    def add_ions(self):
        pass

    def generate_index(self):
        pass

    def minimise(self):
        pass


class PatchCoordinator:

    def __init__(self, config):
        self.config = config
        self.patches = []

        # Process the patch selection settings to determine how many patches to make
        self.frame_range = range(
            self.config["frame_start"],
            self.config["frame_end"],
            self.config["frame_stride"],
        )
        self.resid_range = self.config["patches"]["resid"]

    def copy_simulation_parameters(self):
        # Set up the patch forcefield, lipid topology/parameter files, and minimisation settings
        self.topol = os.path.join(self.config["patches"]["patchdir"], "topol")
        logger.info(f"Copying parameters to {self.topol}")
        if not os.path.exists(self.topol):
            os.makedirs(self.topol)

        # Charmm forcefield
        self.charmm36_ff = os.path.join(
            self.topol, os.path.basename(self.config["filepaths"]["charmm36_ff"])
        )
        shutil.copy2(self.config["filepaths"]["charmm36_ff"], self.charmm36_ff)

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
        for frame in self.frame_range:
            for resid in self.resid_range:
                # Make a patch directory if it doesn't exist
                patchdir = os.path.join(
                    self.config["patches"]["patchdir"], f"patch_{frame}_{resid}"
                )
                if not os.path.exists(patchdir):
                    os.makedirs(patchdir)

                patch = Patch(
                    patchdir=patchdir,
                    outname=f"patch_{frame}_{resid}.gro",
                    cg_traj=self.config["patches"]["traj"],
                    cg_top=self.config["patches"]["top"],
                    patch_frame=frame,
                    patch_resid=resid,
                    patch_size=self.config["patches"]["size"],
                )

                patch.select_patch_residues()
                patch.write_patch()

                # TODO: Refactor so that the CG trajectory is loaded once
                # by the PatchCoordinator and passed to each patch instance
                # instead of loading it for each patch
                del patch.traj
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
                f.write(f"{patch.patchdir}\n")
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
            patch.minimise_in_vacuum()
            patch.apply_correct_stereoisomers()
            patch.restore_residue_numbers()

    def prepare_patches_for_simulation(self):
        for patch in self.patches:
            patch.solvate()
            patch.delete_membrane_waters()
            patch.add_ions()
            patch.generate_index()
            patch.minimise()
