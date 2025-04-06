import argparse
import os.path
import subprocess
import glob
import sys
import json
import MDAnalysis as mda
from io import StringIO
from rdkit import Chem
from schrodinger import adapter
from schrodinger.rdkit import rdkit_adapter
from schrodinger.structutils import analyze
from backmapping.config import load_config, check_fields
from backmapping.logger import logger, check_if_file_exists


def make_tpr_file(basename: str, wdir, mdp, gmx):
    box_args = [
        f"{gmx}",
        "editconf",
        "-f",
        f"{basename}_ini.pdb",
        "-o",
        f"{basename}_box.gro",
        "-d",
        "1",
        "-c",
        "-bt",
        "cubic",
    ]
    subprocess.run(box_args, cwd=wdir)

    if not check_if_file_exists(os.path.join(wdir, f"{basename}_box.gro")):
        sys.exit()

    if not check_if_file_exists(os.path.join(wdir, f"{basename}.top")):
        sys.exit()

    tpr_args = [
        f"{gmx}",
        "grompp",
        "-f",
        f"{mdp}",
        "-c",
        f"{basename}_box.gro",
        "-p",
        f"{basename}.top",
        "-o",
        f"{basename}.tpr",
        "-maxwarn",
        "1",
    ]
    subprocess.run(tpr_args, cwd=wdir)

    if not check_if_file_exists(os.path.join(wdir, f"{basename}.tpr")):
        sys.exit()

    pdb_args = [
        f"{gmx}",
        "editconf",
        "-f",
        f"{basename}.tpr",
        "-o",
        f"{basename}_box.pdb",
        "-conect",
    ]
    subprocess.run(pdb_args, cwd=wdir)

    if not check_if_file_exists(os.path.join(wdir, f"{basename}_box.pdb")):
        sys.exit()

    return os.path.join(wdir, f"{basename}.tpr")


def get_PDBBlock(u, atom_group, frame: int):
    stream = StringIO()
    stream_value = None

    with mda.coordinates.PDB.PDBWriter(stream) as f:
        for ts in u.trajectory[[frame]]:
            f.write(atom_group)
            f._write_pdb_bonds()
            f.END()
        stream_value = stream.getvalue()
    if stream_value == None:
        logger.error("Failed to get PDB block.")
        return
    return stream_value


def get_unsat_stereo(shmol):
    m_unsat = rdkit_adapter.to_rdkit(shmol, implicitH=False)

    unsat_stereo = {}
    i = -1
    for bond in m_unsat.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin = (
                m_unsat.GetAtomWithIdx(bond.GetBeginAtomIdx())
                .GetPDBResidueInfo()
                .GetName()
            )
            end = (
                m_unsat.GetAtomWithIdx(bond.GetEndAtomIdx())
                .GetPDBResidueInfo()
                .GetName()
            )
            i += 1
            unsat_stereo[i] = {
                "begin": begin,
                "end": end,
                "stereo": bond.GetStereo().name[6:],
            }
    return unsat_stereo


def get_amide_stereo(shmol):
    m_amide = rdkit_adapter.to_rdkit(shmol, implicitH=False)

    amide_stereo = {}

    amide_bonds = m_amide.GetSubstructMatches(Chem.MolFromSmarts("[#1]-[#7]-[#6]=[#8]"))

    if len(amide_bonds) == 0:
        return amide_stereo

    i = -1
    for bond in amide_bonds:
        # Get the orientation of this amide bond
        dihedral_angle = Chem.rdMolTransforms.GetDihedralDeg(
            m_amide.GetConformer(0), *bond
        )
        # Convert angles to a label
        if (dihedral_angle > -90) and (dihedral_angle < 90):
            orient = "cis"
        else:
            orient = "trans"

        # Record the orientation
        # Atom ordering in the detected bond respects the order in the SMARTS query above
        i += 1
        amide_stereo[i] = {
            "H": bond[0],
            "N": bond[1],
            "C": bond[2],
            "O": bond[3],
            "stereo": orient,
        }

    return amide_stereo


def get_stereo_labels(filepath, config):
    """
    Get stereoisomer lables for this molecule.
    """
    logger.debug(f"Getting stereo labels for file: {filepath}")

    # Split off the "_ini.pdb" part of the filename to get the usable basename
    basename = os.path.splitext(os.path.basename(filepath))[0]
    if not "_" in basename:
        logger.error(f"Couldn't get basename for file: {filepath}")
        return
    else:
        basename = basename.split("_")[0]

    tpr = make_tpr_file(
        basename,
        config["filepaths"]["aa_gmx_parameters"],
        config["filepaths"]["mdp"],
        config["filepaths"]["gmx"],
    )

    u = mda.Universe(tpr, filepath)

    PDBBlock = get_PDBBlock(u, atom_group=u.select_atoms("all"), frame=0)

    if PDBBlock == None:
        logger.error(f"Couldn't make PDBBlock for: {basename}")
        return

    m_san = Chem.MolFromPDBBlock(
        molBlock=PDBBlock, sanitize=True, removeHs=False, proximityBonding=False
    )

    shmol = rdkit_adapter.from_rdkit(m_san)

    # Schrodinger automatically adds hydrogens to carbons with fewer than 4 bonded neighbors
    # Need to delete these to correctly record double bond stereochemistry
    added_hydrogens = []
    for a in shmol.atom:
        if a.pdbres == "UNK ":
            added_hydrogens.append(a)

    shmol.deleteAtoms(added_hydrogens)

    # Checking for a valid Lewis structure causes Schrodinger to recognize the double bond
    analyze.has_valid_lewis_structure(shmol)
    logger.info("SMILES: ", adapter.to_smiles(shmol, stereo=True))

    detected_stereo = analyze.get_chiral_atoms(shmol)

    # Map atom indices back to atom names
    chiral_stereo = {}
    for idx, stereo in detected_stereo.items():
        atom_name = m_san.GetAtomWithIdx(idx - 1).GetPDBResidueInfo().GetName()
        chiral_stereo[atom_name] = stereo

    # Get double bond stereoisomer
    unsat_stereo = get_unsat_stereo(shmol)

    # Get amide dihedral conformation
    amide_stereo = get_amide_stereo(shmol)

    stereo_reference = {
        basename: {
            "chiral": chiral_stereo,
            "unsat": unsat_stereo,
            "amide": amide_stereo,
        }
    }

    logger.info("Stereoisomer record: ", stereo_reference)

    with open(
        os.path.join(config["filepaths"]["aa_gmx_parameters"], f"{basename}.json"),
        "w",
    ) as f:
        json.dump(stereo_reference, f, indent=4)


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
        "aa_gmx_parameters",
        "mdp",
        "gmx",
    ]
    check_fields(config["filepaths"], required_paths)

    logger.info("Generating stereoisomer reference from coordinate files")
    coord_files = glob.glob(
        os.path.join(config["filepaths"]["aa_gmx_parameters"], "*ini.pdb")
    )

    for file_path in coord_files:
        get_stereo_labels(file_path, config)

    # Combine the stereoisomer records into a single file
    stereo = {}
    stereo_files = glob.glob(
        os.path.join(config["filepaths"]["aa_gmx_parameters"], "*json")
    )
    for stereo_file in stereo_files:
        with open(stereo_file, "r") as f:
            stereo.update(json.load(f))

    with open(os.path.join(config["filepaths"]["stereo"], "stereo.json"), "w") as f:
        json.dump(stereo, f, indent=4)


if __name__ == "__main__":
    logger.info(f"Starting {__file__}")
    main()
    logger.info(f"Finished {__file__}")
