from rdkit import Chem
from backmapping.lipidclass import LipidClass
from backmapping.logger import logger


# Define substructure matching for each lipid class defined in the workflow
def get_substructure_for_lipid_class(LipidClass):
    # NOTE: The substructure matching uses SMARTS patterns without bond valences
    # because the PDB file format doesn't include that information.
    # Note that the hydrogen atoms included in these patterns include implicit hydrogens
    # that are implied by RDKit as it interprets the connectivity (i.e., PH)
    match LipidClass:
        case LipidClass.PC:
            return "[H]C([H])(OC)C([H])(OC)C([H])([H])O[PH](O)(O)OC([H])([H])C([H])([H])[N+](C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H]"
        case LipidClass.PE:
            return "[H]C([H])(OC)C([H])(OC)C([H])([H])O[PH](O)(O)OC([H])([H])C([H])([H])[N+]([H])([H])[H]"
        case LipidClass.PI:
            return "[H]OC1([H])C([H])(O[H])C([H])(O[H])C([H])(O[PH](O)(O)OC([H])([H])C([H])(OC)C([H])([H])OC)C([H])(O[H])C1([H])O[H]"
        case LipidClass.PA:
            return "[H]O[PH](O)(O)OC([H])([H])C([H])(OC)C([H])([H])OC"
        case LipidClass.PS:
            return "[H]C([H])(OC)C([H])(OC)C([H])([H])O[PH](O)(O)OC([H])([H])C([H])(C(O)O)[N+]([H])([H])[H]"
        case LipidClass.PG:
            return "[H]OC([H])([H])C([H])(O[H])C([H])([H])O[PH](O)(O)OC([H])([H])C([H])(OC)C([H])([H])OC"
        case LipidClass.Cer:
            return "[H]OC([H])([H])C([H])(N([H])C(O)C([H])[H])C([H])(C)O[H]"
        case LipidClass.SM:
            return "NCCOPOCCNCO"
        case LipidClass.HexCer:
            return "[H]CC([H])(O[H])C([H])(N([H])C(O)C([H])[H])C([H])([H])OC1([H])OC([H])(C([H])([H])O[H])C([H])(O[H])C([H])(O[H])C1([H])O[H]"
        case _:
            return ""


def detect_lipid_class(structure: Chem.Mol):
    # Try each defined lipid class's substructure pattern
    # If we don't get a unique match, complain
    matched_classes = []
    for lipid_class in LipidClass:

        pattern = get_substructure_for_lipid_class(lipid_class)

        matches = structure.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        if matches != ():
            logger.debug(f"Match for {lipid_class}")
            matched_classes.append(lipid_class)
        else:
            continue

    if len(matched_classes) == 0:
        logger.error("No lipid class matched")
        return None
    if len(matched_classes) != 1:
        logger.error(f"Matched multiple lipid classes: {matched_classes}")
        return None
    else:
        return matched_classes[0]
