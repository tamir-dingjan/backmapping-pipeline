import os
import glob
from rdkit import Chem
from backmapping import substructure
from backmapping.io_utils import load_structure_from_pdb_file
from backmapping.lipidclass import LipidClass


def test_substructure_detection_for_whole_lipid():
    coord_path = os.path.join(
        os.path.dirname(__file__), "data", "aa_gmx", "PSM_box.pdb"
    )
    structure = load_structure_from_pdb_file(coord_path)

    detected_lipid_class = substructure.detect_lipid_class(structure)

    assert detected_lipid_class == LipidClass.SM


def test_substructure_detection_for_PC_headgroup():
    coord_path = os.path.join(
        os.path.dirname(__file__), "data", "lipid_class_substructures", "PC.pdb"
    )
    structure = load_structure_from_pdb_file(coord_path)

    detected_lipid_class = substructure.detect_lipid_class(structure)

    assert detected_lipid_class == LipidClass.PC


def test_substructure_detection_for_PE_headgroup():
    coord_path = os.path.join(
        os.path.dirname(__file__), "data", "lipid_class_substructures", "PE.pdb"
    )
    structure = load_structure_from_pdb_file(coord_path)

    detected_lipid_class = substructure.detect_lipid_class(structure)

    assert detected_lipid_class == LipidClass.PE


def test_substructure_detection_for_PA_headgroup():
    coord_path = os.path.join(
        os.path.dirname(__file__), "data", "lipid_class_substructures", "PA.pdb"
    )
    structure = load_structure_from_pdb_file(coord_path)

    detected_lipid_class = substructure.detect_lipid_class(structure)

    assert detected_lipid_class == LipidClass.PA


def test_substructure_detection_for_PG_headgroup():
    coord_path = os.path.join(
        os.path.dirname(__file__), "data", "lipid_class_substructures", "PG.pdb"
    )
    structure = load_structure_from_pdb_file(coord_path)

    detected_lipid_class = substructure.detect_lipid_class(structure)

    assert detected_lipid_class == LipidClass.PG


def test_substructure_detection_for_PI_headgroup():
    coord_path = os.path.join(
        os.path.dirname(__file__), "data", "lipid_class_substructures", "PI.pdb"
    )
    structure = load_structure_from_pdb_file(coord_path)

    detected_lipid_class = substructure.detect_lipid_class(structure)

    assert detected_lipid_class == LipidClass.PI


def test_substructure_detection_for_PS_headgroup():
    coord_path = os.path.join(
        os.path.dirname(__file__), "data", "lipid_class_substructures", "PS.pdb"
    )
    structure = load_structure_from_pdb_file(coord_path)

    detected_lipid_class = substructure.detect_lipid_class(structure)

    assert detected_lipid_class == LipidClass.PS


def test_substructure_detection_for_SM_headgroup():
    coord_path = os.path.join(
        os.path.dirname(__file__), "data", "lipid_class_substructures", "SM.pdb"
    )
    structure = load_structure_from_pdb_file(coord_path)

    detected_lipid_class = substructure.detect_lipid_class(structure)

    assert detected_lipid_class == LipidClass.SM


def test_substructure_detection_for_Cer_headgroup():
    coord_path = os.path.join(
        os.path.dirname(__file__), "data", "lipid_class_substructures", "Cer.pdb"
    )
    structure = load_structure_from_pdb_file(coord_path)

    detected_lipid_class = substructure.detect_lipid_class(structure)

    assert detected_lipid_class == LipidClass.Cer


def test_substructure_detection_for_deoxyCer_headgroup():
    coord_path = os.path.join(
        os.path.dirname(__file__), "data", "lipid_class_substructures", "DeoxyCer.pdb"
    )
    structure = load_structure_from_pdb_file(coord_path)

    detected_lipid_class = substructure.detect_lipid_class(structure)

    assert detected_lipid_class == LipidClass.DeoxyCer


def test_substructure_detection_for_HexCer_headgroup():
    coord_path = os.path.join(
        os.path.dirname(__file__), "data", "lipid_class_substructures", "HexCer.pdb"
    )
    structure = load_structure_from_pdb_file(coord_path)

    detected_lipid_class = substructure.detect_lipid_class(structure)

    assert detected_lipid_class == LipidClass.HexCer


def test_substructure_detection_for_sterol():
    coord_path = os.path.join(
        os.path.dirname(__file__), "data", "lipid_class_substructures", "sterol.pdb"
    )
    structure = load_structure_from_pdb_file(coord_path)

    detected_lipid_class = substructure.detect_lipid_class(structure)

    assert detected_lipid_class == LipidClass.Sterol


if __name__ == "__main__":
    test_substructure_detection_for_sterol()
    test_substructure_detection_for_deoxyCer_headgroup()
    test_substructure_detection_for_Cer_headgroup()
    test_substructure_detection_for_HexCer_headgroup()
    test_substructure_detection_for_PC_headgroup()
    test_substructure_detection_for_PE_headgroup()
    test_substructure_detection_for_PA_headgroup()
    test_substructure_detection_for_PG_headgroup()
    test_substructure_detection_for_PI_headgroup()
    test_substructure_detection_for_PS_headgroup()
    test_substructure_detection_for_SM_headgroup()
    test_substructure_detection_for_whole_lipid()
