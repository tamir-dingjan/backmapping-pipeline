from backmapping.logger import logger


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


def get_nonmatching_chiral_centers(stereo, stereo_lookup):
    mismatched_chiral_centers = []

    for chiral_name in stereo_lookup["chiral"].keys():
        if stereo["chiral"][chiral_name] != stereo_lookup["chiral"][chiral_name]:
            mismatched_chiral_centers.append(chiral_name)

    return mismatched_chiral_centers
