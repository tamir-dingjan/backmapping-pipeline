import numpy as np
from backmapping.logger import logger
from backmapping.stereo.utils import get_plane_equation_from_points, mirror_point


class InvalidPatchStereo(Exception):
    """Raised when a patch doesn't contain the same stereo records as the lookup"""

    def __init__(self, message="Invalid patch stereoconformation detected"):
        self.message = f"{message}"
        super().__init__(self.message)


def validate_patch_stereoconformation(stereo, stereo_lookup):
    # Check that all chiral centers, unsaturations, and amide bonds
    # defined in the stereo lookup are also detected in the patch
    if sorted(stereo["chiral"].keys()) != sorted(stereo_lookup["chiral"].keys()):
        logger.warning("Chiral records are invalid for this patch")
        return False

    # Check that the same number of unsaturated bonds are listed
    if sorted(stereo["unsat"].keys()) != sorted(stereo_lookup["unsat"].keys()):
        logger.warning("Unsat records are invalid for this patch")
        return False

    # Check that each bond in the lookup is defined in the patch
    for unsat_lookup in stereo_lookup["unsat"].values():
        valid = False
        for unsat_patch in stereo["unsat"].values():
            if (unsat_lookup["begin"] == unsat_patch["begin"]) and (
                unsat_lookup["end"] == unsat_patch["end"]
            ):
                valid = True
                continue
            # Also consider cases where the bond is recorded in the reverse order
            elif (unsat_lookup["begin"] == unsat_patch["end"]) and (
                unsat_lookup["end"] == unsat_patch["begin"]
            ):
                valid = True
                continue
        if not valid:
            logger.warning(f"Unsat record not found in lookup: {unsat_patch}")
            return False

    # Check that each amide bond in the lookup is defined in the patch
    if stereo["amide"].keys() != stereo_lookup["amide"].keys():
        logger.warning("Amide records are invalid for this patch")
        return False

    return True


def get_nonmatching_chiral_centers(stereo, stereo_lookup):
    mismatched_chiral_centers = []

    for chiral_name in stereo_lookup["chiral"].keys():
        if stereo["chiral"][chiral_name] != stereo_lookup["chiral"][chiral_name]:
            mismatched_chiral_centers.append(chiral_name)

    return mismatched_chiral_centers


def get_chiral_proton(chiral_atom, u):
    bonded_protons = u.select_atoms(f"element H and bonded (index {chiral_atom.index})")
    if len(bonded_protons) != 1:
        logger.debug(f"Couldn't identify a chiral proton for atom {chiral_atom.name}")
    return bonded_protons


def flip_chiral_proton(chiral_atom, chiral_proton, heavy_nbors, u, positions):
    # Flip the chiral center and the attached proton in the plane defined by
    # the chiral center's other bonded neighbors
    # However, if this chiral center is in a carbohydrate, then moving the chiral center
    # may change the carbohydrate ring conformer, which we do not want. In this case
    # we can just swap the positions of the proton and the substituent hydroxyl oxygen atoms.
    check_single_oxygen_nbor = sorted([x.element for x in heavy_nbors]) == [
        "C",
        "C",
        "O",
    ]
    check_hydroxyl_group = False

    if check_single_oxygen_nbor:
        bonded_oxygens = u.select_atoms(
            "element O and bonded (index %s)" % chiral_atom.index
        )
        check_hydroxyl_group = sorted(
            [x.element for x in bonded_oxygens[0].bonded_atoms]
        ) == ["C", "H"]

    if check_hydroxyl_group:
        logger.debug("Chiral center has a hydroxyl group.")
        # Swap positions of the hydroxyl oxygen and the chiral proton

        hydroxyl_oxygen = bonded_oxygens[0]
        logger.debug(
            "Swapping positions of: %s <-> %s"
            % (chiral_proton.name, hydroxyl_oxygen.name)
        )
        swap_vec = positions[chiral_proton.index] - positions[hydroxyl_oxygen.index]
        shifted_hydroxyl_position = positions[hydroxyl_oxygen.index] + swap_vec
        for hydroxyl_proton in u.select_atoms(
            "element H and bonded (index %s)" % hydroxyl_oxygen.index
        ):
            positions[hydroxyl_proton.index] = (
                positions[hydroxyl_proton.index] + swap_vec
            )
        positions[chiral_proton.index] = positions[hydroxyl_oxygen.index]
        positions[hydroxyl_oxygen.index] = shifted_hydroxyl_position

    else:
        # Flip the chiral center and chiral proton in the plane defined
        # by the heavy bonded neighbors
        logger.debug("Stamping: %s & %s" % (chiral_proton.name, chiral_atom.name))
        nbor_coords = u.atoms.positions[heavy_nbors.indices]
        a, b, c, d = get_plane_equation_from_points(nbor_coords)

        positions[chiral_proton.index] = mirror_point(
            a, b, c, d, chiral_proton.position
        )
        positions[chiral_atom.index] = mirror_point(a, b, c, d, chiral_atom.position)
    return positions


def flip_chiral_methyl(chiral_atom, heavy_nbors, u, positions):
    methyl = []
    non_methyl = []
    for nbor in heavy_nbors:
        if sorted([x.element for x in nbor.bonded_atoms]) == [
            "C",
            "H",
            "H",
            "H",
        ]:
            methyl.append(nbor)
        else:
            non_methyl.append(nbor)
    if (len(methyl) != 1) and (len(non_methyl) != 3):
        logger.error(
            "Did not identify a single methyl and 3 non-methyl neighbors for chiral center: {chiral_atom.name}"
        )
        return

    # Reflect the methyl group and chiral center in the plane defined by the chiral center's other bonded neighbors
    # Note that non_methyl is a list not an AtomGroup here
    logger.debug(f"Stamping chiral center: {chiral_atom.name}")

    non_methyl_coords = u.atoms.positions[[x.index for x in non_methyl]]

    a, b, c, d = get_plane_equation_from_points(non_methyl_coords)

    positions[chiral_atom.index] = mirror_point(a, b, c, d, chiral_atom.position)

    methyl_carbon = methyl[0]
    methyl_protons = u.select_atoms(
        "element H and bonded (index %s )" % methyl_carbon.index
    )
    logger.debug(
        "Stamping methyl group: %s"
        % (methyl_carbon.name + " " + " ".join(methyl_protons.names))
    )
    positions[methyl_carbon.index] = mirror_point(a, b, c, d, methyl_carbon.position)
    for proton in methyl_protons:
        positions[proton.index] = mirror_point(a, b, c, d, proton.position)
    return positions


def correct_chiral_center(chiral_name, u, res, positions):
    chiral_atom = u.select_atoms(f"resid {str(res)} and name {chiral_name}")
    if len(chiral_atom) != 1:
        logger.warning(
            f"Could not identify a unique chiral center with name: {chiral_name}"
        )
    else:
        chiral_atom = chiral_atom[0]

    # Check if there is an attached proton
    bonded_protons = get_chiral_proton(chiral_atom, u)
    heavy_nbors = u.select_atoms(
        f"(not element H) and bonded (index {chiral_atom.index})"
    )

    if len(bonded_protons) == 1:
        logger.debug(f"Found proton bonded to chiral center: {bonded_protons.names[0]}")
        chiral_proton = bonded_protons[0]
        positions = flip_chiral_proton(
            chiral_atom, chiral_proton, heavy_nbors, u, positions
        )

    elif len(bonded_protons) == 0:
        # This chiral center only has bonds to carbon atoms
        # Process this as part of a sterol
        positions = flip_chiral_methyl(chiral_atom)

    else:
        # Chiral center has multiple protons - apparently...
        logger.error(f"Found multiple protons at chiral center: {chiral_atom.name}")
    return positions


def get_nonmatching_unsat(stereo, stereo_lookup):
    # The order of unsaturated bond stereo records in the lookup and in the patch
    # may differ. So check all-to-all
    mismatches = []
    for unsat_patch in stereo["unsat"].values():
        matched = False
        for unsat_lookup in stereo_lookup["unsat"].values():
            if unsat_patch == unsat_lookup:
                matched = True
            # Also consider the case where the begin and end atoms are labelled
            # in reverse - as long as the stereoconformer matches, this is still OK
            elif (
                (unsat_patch["begin"] == unsat_lookup["end"])
                and (unsat_patch["end"] == unsat_lookup["begin"])
                and (unsat_patch["stereo"] == unsat_lookup["stereo"])
            ):
                matched = True
        if not matched:
            mismatches.append(unsat_patch)
    return mismatches


def correct_unsat(bond, u, res, positions):
    begin_name = bond["begin"]
    end_name = bond["end"]
    logger.debug(f"Unmatched bond stereoconfonformer: {begin_name} - {end_name}")

    begin_atom = u.select_atoms(f"resid {res} and name {begin_name}")
    end_atom = u.select_atoms(f"resid {res} and name {end_name}")

    if (len(begin_atom) != 1) and (len(end_atom) != 1):
        logger.warning(f"Could not identify bond atoms: {begin_name} - {end_name}")
        return

    begin_atom = begin_atom[0]
    end_atom = end_atom[0]

    # Find a substituent proton
    for atom1, atom2 in zip([begin_atom, end_atom], [end_atom, begin_atom]):
        proton = u.select_atoms(f"element H and bonded (index {atom1.index})")
        substituent = u.select_atoms(
            f"not (element H) and not (index {atom2.index}) and bonded (index {atom1.index})"
        )
        if (len(proton) == 1) and (len(substituent) == 1):
            break
    if len(proton) != 1:
        logger.warning("Could not find unsaturated proton to swap with")
        return
    elif len(substituent) != 1:
        logger.warning("Could not find non-proton substituent to swap with")
        return

    proton = proton[0]
    substituent = substituent[0]

    logger.debug(
        f"Swapping positions of proton ({proton.name}) and substituent ({substituent.name})"
    )

    swap_vec = positions[proton.index] - positions[substituent.index]
    shifted_substituent_position = positions[substituent.index] + swap_vec
    for substituent_proton in u.select_atoms(
        f"element H and bonded (index {substituent.index})"
    ):
        positions[substituent_proton.index] = (
            positions[substituent_proton.index] + swap_vec
        )
    positions[proton.index] = positions[substituent.index]
    positions[substituent.index] = shifted_substituent_position
    return positions


def get_nonmatching_amide(stereo, stereo_lookup):
    mismatches = []
    for amide_patch in stereo["amide"].values():
        matched = False
        for amide_lookup in stereo_lookup["amide"].values():
            if amide_patch == amide_lookup:
                matched = True
        if not matched:
            mismatches.append(amide_patch)
    return mismatches


def correct_amide(amide, u, res, positions):
    # Correct by reflecting the N and attached proton in a plane defined by
    # the two neighbouring carbon atoms (C1, C2) and the point midway between them.
    # The midway point is on a straight line between C1 and C2, so
    # a single plane cannot be defined by these points.
    # To allow a single plane definition, we have to displace the midway point (X)
    # in a direction orthogonal to the C1-C2 vector and orthogonal to the
    # X-N vector. This will define the plane as normal to the N-H bond, which is
    # useful for reflecting the N and H atoms.

    # The amide atoms are stored in the stereo lookup by their element and 0-based
    # index within the lipid. This index does not match the index within the patch,
    # so use this number to select within the lipid atom group.
    nitrogen = u.select_atoms(f"resid {res}")[amide["N"]]
    carbon_1 = u.select_atoms(f"resid {res}")[amide["C"]]
    for carbon in u.select_atoms(f"element C and bonded (index {nitrogen.index})"):
        if carbon.index == carbon_1.index:
            continue
        carbon_2 = carbon
    proton = u.select_atoms(f"element H and bonded (index {nitrogen.index})")[0]

    # Construct the midway point
    c2_c1_vec = positions[carbon_1.index] - positions[carbon_2.index]
    midway = positions[carbon_2.index] + 0.5 * c2_c1_vec

    # Find the vector orthogonal to C1-C2 and orthogonal to X-N
    n_midway_vec = midway - positions[nitrogen.index]

    orthogonal_vec = np.cross(c2_c1_vec, n_midway_vec)

    # Displace the midway point along the orthogonal vector
    # It doesn't matter whether in the +ve or -ve direction
    # since we only need this point for the plane construction
    displaced_midway = midway + orthogonal_vec

    # Construct the plane

    plane_coords = np.asarray(
        [positions[carbon_1.index], positions[carbon_2.index], displaced_midway],
    )

    a, b, c, d = get_plane_equation_from_points(plane_coords)

    positions[nitrogen.index] = mirror_point(a, b, c, d, nitrogen.position)
    positions[proton.index] = mirror_point(a, b, c, d, proton.position)
    return positions


def correct_stereo(stereo, stereo_lookup, u, res, positions):

    valid_patch = validate_patch_stereoconformation(stereo, stereo_lookup)
    if not valid_patch:
        logger.warning(
            "Patch stereoconformation is not valid. Cannot correct stereoconformation."
        )
        raise InvalidPatchStereo

    # Corect amides first, because the reflection from cis to trans amide
    # makes a larger change when beginning from a minimised conformation
    # Correct amide conformations
    if not stereo["amide"] == {}:
        mismatched_amide = get_nonmatching_amide(stereo, stereo_lookup)
        if mismatched_amide == []:
            logger.debug("Amide conformations match")
        else:
            for amide in mismatched_amide:
                positions = correct_amide(amide, u, res, positions)

    # Correct individual chiral centers
    mismatched_chiral_centers = get_nonmatching_chiral_centers(stereo, stereo_lookup)
    if mismatched_chiral_centers == []:
        logger.debug("Chiral centers match")
    else:
        for chiral_name in mismatched_chiral_centers:
            positions = correct_chiral_center(chiral_name, u, res, positions)

    # Correct unsaturations
    mismatched_unsat = get_nonmatching_unsat(stereo, stereo_lookup)
    if mismatched_unsat == []:
        logger.debug("Unsaturation stereoconfigurations match")
    else:
        for bond in mismatched_unsat:
            positions = correct_unsat(bond, u, res, positions)

    return positions
