import numpy as np
from io import StringIO
import MDAnalysis as mda
from rdkit import Chem
from schrodinger import adapter
from schrodinger.rdkit import rdkit_adapter
from schrodinger.structutils import analyze
from backmapping.logger import logger


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
            # Store these with a str to match the lookup read from json
            unsat_stereo[str(i)] = {
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
        # Store these with a str to match the lookup read from json
        amide_stereo[str(i)] = {
            "H": bond[0],
            "N": bond[1],
            "C": bond[2],
            "O": bond[3],
            "stereo": orient,
        }

    return amide_stereo


def get_stereo(PDBBlock):
    if PDBBlock is None:
        logger.error("No PDBBlock provided.")
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
    logger.info(f"SMILES: {adapter.to_smiles(shmol, stereo=True)}")

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

    stereoconformer = {
        "chiral": chiral_stereo,
        "unsat": unsat_stereo,
        "amide": amide_stereo,
    }
    return stereoconformer


def is_same_bonded_stereo(a, b):
    """
    The unsaturation stereochemistry blocks are written with a per-bond dictionary, e.g.,:
        "0": {
            "begin": " O13",
            "end": " P  ",
            "stereo": "NONE"
        },
        "1": {
            "begin": " C5S",
            "end": " C4S",
            "stereo": "E"
        },
        "2": {
            "begin": " O22",
            "end": " C21",
            "stereo": "NONE"
        }

    Likewise, the amide stereochemistry blocks:
        "0": {
            "H": 78,
            "N": 27,
            "C": 79,
            "O": 80,
            "stereo": "trans"
        }

    However, the order of these bonds may not be the same in patches and in the lookup.
    So, to compare the two, we need to compare each bond independently.

    """

    for i in a:
        match_bond = False
        for j in b:
            if i == j:
                match_bond = True
        if not match_bond:
            return False

    return True


def get_plane_equation_from_points(points):
    if points.shape != (3, 3):
        logger.error("Can't compute plane from input of shape: " + str(points.shape))
        return

    x1, y1, z1 = points[0]
    x2, y2, z2 = points[1]
    x3, y3, z3 = points[2]
    a1 = x2 - x1
    b1 = y2 - y1
    c1 = z2 - z1
    a2 = x3 - x1
    b2 = y3 - y1
    c2 = z3 - z1
    a = b1 * c2 - b2 * c1
    b = a2 * c1 - a1 * c2
    c = a1 * b2 - b1 * a2
    d = -a * x1 - b * y1 - c * z1
    return a, b, c, d


def mirror_point(a, b, c, d, point):
    x1 = point[0]
    y1 = point[1]
    z1 = point[2]
    k = (-a * x1 - b * y1 - c * z1 - d) / float((a * a + b * b + c * c))
    x2 = a * k + x1
    y2 = b * k + y1
    z2 = c * k + z1
    x3 = 2 * x2 - x1
    y3 = 2 * y2 - y1
    z3 = 2 * z2 - z1

    return np.asarray([x3, y3, z3])
