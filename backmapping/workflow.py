from datetime import datetime
import os
import numpy as np
import networkx as nx
import fnmatch
from rdkit import Chem
from collections import defaultdict
from backmapping.logger import logger
from backmapping.lipidclass import LipidClass, unsupported_lipid_classes
from backmapping.substructure import detect_lipid_class
from backmapping import io_utils
from backmapping.utils import *


class BackmappingWorkflow:

    def __init__(self, coord, basename, aa, cg, stereo, outfile):
        logger.debug("Instantiating BackmappingWorkflow")
        self.basename = basename
        self.aa = None
        self.cg = None
        self.stereo = None
        self.outfile = outfile
        self.outfile_lines = []
        self.aa_regions = None
        self.cg_regions = None
        self.lipid_class = LipidClass
        self.bead_c1 = ""
        self.bead_c2 = ""
        self.bead_c3 = ""
        self.ceramide_nitrogen = ""
        self.seen_nodes = []
        self.is_supported_lipid_class = None

        try:
            # Load the lipid coordinates, topology files, and stereo reference
            if os.path.splitext(coord)[1] == ".mol2":
                self.structure = io_utils.load_structure_from_mol2_file(coord)
            elif os.path.splitext(coord)[1] == ".pdb":
                self.structure = io_utils.load_structure_from_pdb_file(coord)

            self.aa = io_utils.load_itp_as_network(aa)
            self.cg = io_utils.load_itp_as_network(cg)
            self.cg_bead_order = io_utils.get_bead_order_from_itp(cg)
            self.stereo = io_utils.load_stereo_reference(stereo)
        except Exception as e:
            logger.error(
                f"Couldn't instantiate the workflow for the following input: \ncoord: {coord}\nall-atom topology: {aa}\ncoarse-grained topology: {cg}"
            )
            raise (e)

    # Set lipid class
    def set_lipid_class(self):
        self.lipid_class = detect_lipid_class(self.structure)
        if self.lipid_class == None:
            msg = "Error! Couldn't detect lipid class. Check the input structure."
            logger.error(msg)
            raise Exception(msg)

        # Note if we have an unsupported lipid class (currently limited to sterols)
        if self.lipid_class in unsupported_lipid_classes:
            logger.warning("Detected unsupported lipid class: {self.lipid_class}.")
            self.is_supported_lipid_class = False
        else:
            self.is_supported_lipid_class = True

        # Additionally, if the lipid class is a sphingolipid, set the ceramide nitrogen
        if self.lipid_class in (
            LipidClass.SM,
            LipidClass.Cer,
            LipidClass.HexCer,
            LipidClass.DeoxyCer,
        ):
            # Set the ceramide nitrogen identity
            # TODO: Edit all references to ceramide nitrogen to use this variable
            self.ceramide_nitrogen = get_ceramide_nitrogen(self.aa)

            # If HexCer, add in the bonds for the head group beads and set the bead labels
            if self.lipid_class == LipidClass.HexCer:

                # Add in the hexose bonds for hexosylceramides
                self.bead_c1 = select_bead_by_name("C1", self.cg)
                self.bead_c2 = select_bead_by_name("C2", self.cg)
                self.bead_c3 = select_bead_by_name("C3", self.cg)

                self.cg.add_edge(self.bead_c1[0], self.bead_c2[0])
                self.cg.add_edge(self.bead_c1[0], self.bead_c3[0])
                self.cg.add_edge(self.bead_c3[0], self.bead_c2[0])

    # Assign regions to atoms and beads
    def assign_regions(self):
        self.assign_aa_regions()
        self.assign_cg_regions()

    def assign_aa_regions(self):
        """
        Adds region properties for topology atoms.

        Tail regions: atoms which are C or H with only C or H neighbors. Hydrogens
        must also have no heteroatoms among their 2nd neighbors.
        """

        region_lookup = defaultdict(list)
        network = self.aa

        if self.lipid_class in (
            LipidClass.SM,
            LipidClass.Cer,
            LipidClass.DeoxyCer,
            LipidClass.HexCer,
        ):
            ceramide_c2 = None
        else:
            tailA_dist = 6
            tailB_dist = 7

        # Assign tail region if the atom is not hetero and all neighbors are not hetero
        for n in network.nodes:
            if (not network.nodes[n]["hetero"]) and (
                not has_heteroatom_neighbors(network, n)
            ):
                network.nodes[n]["region"] = "generic_tail"
                region_lookup["generic_tail"].append(n)
                # The above step will also mark head group aliphatic hydrogens as
                # tail. Catch these
                if (network.nodes[n]["element"] == "H") and bool(
                    sum(
                        [
                            has_heteroatom_neighbors(network, nbor)
                            for nbor in network.neighbors(n)
                        ]
                    )
                ):
                    network.nodes[n]["region"] = "unknown"
                    region_lookup["generic_tail"].remove(n)

            else:
                network.nodes[n]["region"] = "unknown"

        # If we have a deoxyceramide, we need to make sure C1 is labelled as "unknown" instead of "generic_tail"
        if self.lipid_class in (
            LipidClass.SM,
            LipidClass.Cer,
            LipidClass.HexCer,
            LipidClass.DeoxyCer,
        ):
            # Using path termini will select C1, C3, and the alpha-carbon of the acyl chain
            for target in get_path_termini(
                "element", ["C", "C"], get_ceramide_nitrogen(network), network
            ):
                # If this target has three hydrogens, remove it and them from the generic tail group
                if [
                    network.nodes[n]["element"] for n in network.neighbors(target)
                ].count("H") == 3:
                    if target in region_lookup["generic_tail"]:
                        region_lookup["generic_tail"].remove(target)
                        network.nodes[target]["region"] = "unknown"
                        # Also remove attached hydrogens
                        for nbor in network.neighbors(target):
                            if (network.nodes[nbor]["element"] == "H") & (
                                nbor in region_lookup["generic_tail"]
                            ):
                                region_lookup["generic_tail"].remove(nbor)
                                network.nodes[nbor]["region"] = "unknown"

        # Assign head and linker carbons.
        # The headgroup phosphate divides the linker and head regions.
        # All lipids except for ceramide have a phosphate in the head group.
        # So P and neighboring atoms are in the head region.
        for n in network.nodes:
            if (network.nodes[n]["element"] == "P") or (
                "P" in [network.nodes[nbor]["element"] for nbor in network.neighbors(n)]
            ):
                network.nodes[n]["region"] = "head"
                region_lookup["head"].append(n)
                if network.nodes[n]["element"] == "P":
                    region_lookup["phosphate"].append(n)

        # Assign hexosylceramide head group boundary C11 carbon
        if self.lipid_class == LipidClass.HexCer:
            # Find this carbon
            target = get_path_termini(
                "element", ["C", "C", "O", "C"], get_ceramide_nitrogen(network), network
            )
            if type(target) != str:
                raise Exception(
                    "Error! Couldn't find the hexosylceramide ether-linked ring carbon. Stopping."
                )
            else:
                network.nodes[target]["region"] = "head"
                region_lookup["head"].append(target)

        # Assign the linker region
        # Defined as being between the tail and head.
        # Need to distinguish between tail carbonyl carbon and glycerol.
        # Also, for PG lipids, need to distinguish between the two glycerol groups.
        # So, build from the tail
        for n in network.nodes:
            if (network.nodes[n]["region"] == "unknown") and (
                "generic_tail"
                in [network.nodes[nbor]["region"] for nbor in network.neighbors(n)]
            ):
                # In sphingolipids and ether-linked lipids, this atom will have only a single oxygen neighbor
                # In glycerophospholipids, this atom will have the carbonyl and ether oxygen neighbors
                network.nodes[n]["region"] = "linker"
                region_lookup["linker"].append(n)

        # Expand out the linker region between the tail and head
        while check_unassigned_linker_atoms(network):
            for n in network.nodes:
                if check_prelinker_atom(network, n):
                    network.nodes[n]["region"] = "linker"
                    region_lookup["linker"].append(n)

        # Assign remaining unknown atoms as the head group
        # For all lipids containing a P atom, this will be the remaining atoms
        # on the other side of the phosphate from the acyl chains.
        # For ceramides, which have no P atom, all the unknown atoms will already
        # be labelled as linker.
        for n in network.nodes:
            if network.nodes[n]["region"] == "unknown":
                network.nodes[n]["region"] = "head"
                region_lookup["head"].append(n)

        # Assign tailA and tailB
        # tailA is the tail connected to glycerol sn-2 position, closer to the phosphate.
        # So, begin by identifying this tail.

        while check_unassigned_generic_tail(network):
            for n in region_lookup["generic_tail"]:

                # Propagate existing tail labels
                if "tailA" in [
                    network.nodes[nbor]["region"] for nbor in network.neighbors(n)
                ]:
                    network.nodes[n]["region"] = "tailA"
                    region_lookup["tailA"].append(n)

                elif "tailB" in [
                    network.nodes[nbor]["region"] for nbor in network.neighbors(n)
                ]:
                    network.nodes[n]["region"] = "tailB"
                    region_lookup["tailB"].append(n)

                # If we're not a neighbor to tailA or tailB, we need to consider applying these labels de novo
                else:
                    # For glycerolipids, use distance to phosphate - 6-bond distance for tailA, 7 for tailB
                    if not self.lipid_class in (
                        LipidClass.SM,
                        LipidClass.Cer,
                        LipidClass.HexCer,
                        LipidClass.DeoxyCer,
                    ):

                        dist = nx.algorithms.shortest_path_length(
                            network, source=n, target=region_lookup["phosphate"][0]
                        )

                        if "linker" in [
                            network.nodes[nbor]["region"]
                            for nbor in network.neighbors(n)
                        ]:
                            if dist == tailA_dist:
                                network.nodes[n]["region"] = "tailA"
                                region_lookup["tailA"].append(n)

                            elif dist == tailB_dist:
                                network.nodes[n]["region"] = "tailB"
                                region_lookup["tailB"].append(n)
                    else:
                        # For sphingolipids, use the comparison between C2 and N distances.
                        # This more complex definition covers the multiple types of sphingolipids
                        # with varying linker region structure: sphingomyelins, ceramides, and phytoceramides.# If we reach here, we're dealing with a sphingolipid.
                        # So calculate the distance between the N2 amide nitrogen and the tail atoms
                        ceramide_nitrogen = get_ceramide_nitrogen(network)

                        # Find the C2 atom
                        for nbor in network.neighbors(ceramide_nitrogen):
                            # If we have an oxygen neighbor, or this is a proton, it's not C2
                            if (network.nodes[nbor]["element"] == "H") or (
                                "O"
                                in [
                                    network.nodes[i]["element"]
                                    for i in network.neighbors(nbor)
                                ]
                            ):
                                continue
                            else:
                                ceramide_c2 = nbor

                        # If we couldn't identify a C2, we have a problem
                        if ceramide_c2 == None:
                            raise Exception(
                                "Error! Couldn't find C2 in this sphingolipid! Stopping."
                            )

                        # Compare the N and C2 distances
                        # If dist_C2 < dist_N, the atom is on tailA
                        # If dist_N < dist_C2, the atom is on tailB
                        dist_n = nx.algorithms.shortest_path_length(
                            network, source=n, target=ceramide_nitrogen
                        )
                        dist_c2 = nx.algorithms.shortest_path_length(
                            network, source=n, target=ceramide_c2
                        )

                        if dist_c2 < dist_n:
                            network.nodes[n]["region"] = "tailA"
                            region_lookup["tailA"].append(n)
                        else:
                            network.nodes[n]["region"] = "tailB"
                            region_lookup["tailB"].append(n)

        self.aa_regions = region_lookup

    def assign_cg_regions(self):

        region_lookup = defaultdict(list)

        # Select regions based on bead name matching
        for n in self.cg.nodes:

            if fnmatch.fnmatch(self.cg.nodes[n]["name"], "[CDT][123456]A"):
                self.cg.nodes[n]["region"] = "tailA"
                region_lookup["tailA"].append(n)

            elif fnmatch.fnmatch(self.cg.nodes[n]["name"], "[CDT][123456]B"):
                self.cg.nodes[n]["region"] = "tailB"
                region_lookup["tailB"].append(n)

            elif fnmatch.fnmatch(self.cg.nodes[n]["name"], "[GA][LM][12]"):
                self.cg.nodes[n]["region"] = "linker"
                region_lookup["linker"].append(n)

            else:
                self.cg.nodes[n]["region"] = "head"
                region_lookup["head"].append(n)

        self.cg_regions = region_lookup

    # Allocate atoms to beads within the head group region
    def allocate_head(self):
        # For PC, PE, PG, and PS, use a single phosphate-attached bead to map to all
        # non-phosphate head group atoms.
        # This mapping also applies to sphingmyelin head groups.
        # The linking atoms are split partially between the phosphate and the remainder.
        #
        #  P --- O --- CH2 --- CH2 --- X --->
        #
        # PO4   PO4    PO4     PO4    XXX
        #       PO4    XXX     XXX
        #       PO4            XXX
        #       XXX            XXX
        #
        # The assignment is different for PA (all head group is PO4 bead):
        #
        #
        #  P --- O --- H
        #
        # PO4   PO4   PO4
        #
        #
        # The case for PI lipids is more complex:
        #
        #
        #                 C1      C3
        #                 C3
        #
        #                CHOH --- CHOH
        #               /          \
        #  P --- O --- CH           CHOH    C2
        #               \          /        C3
        #                CHOH --- CHOH
        #
        # PO4   PO4   C1  C1      C2
        #        C1       C2
        #
        #
        # If this lipid contains a phosphate in the head group,
        # the phosphate P and O atoms are mapped to the PO4 bead
        phosphate_bead = select_bead_by_name("PO4", self.cg)
        if "phosphate" in self.aa_regions.keys():
            # TODO: Is atom index 0 of the "phosphate" region always the phosphorous atom?
            self.aa.nodes[self.aa_regions["phosphate"][0]]["map"] = phosphate_bead

            for n in get_deadend_nbors(self.aa_regions["phosphate"][0], self.aa):
                self.aa.nodes[n]["map"] = phosphate_bead

        if self.lipid_class == LipidClass.PI:
            self.allocate_head_pi()
        elif self.lipid_class == LipidClass.PA:
            self.allocate_head_pa()
        elif self.lipid_class == LipidClass.HexCer:
            self.allocate_head_hexcer()
        elif (
            self.lipid_class == LipidClass.Cer
            or self.lipid_class == LipidClass.DeoxyCer
        ):
            # Ceramide lipids have no head group region, so no need to map any beads
            pass
        else:
            # One of the other phospholipid classes
            # Get the first CH2 atom by following the path from the linker through the phosphate
            aa_target = get_path_termini(
                "element",
                ["O", "P", "O", "C"],
                get_single_region_nbor("linker", "head", self.aa),
                self.aa,
            )

            # Select the X bead attached to the PO4
            # Note that this won't select the PI C1 bead, since that is not a dead-end neighbor
            bead_map = get_deadend_nbors(
                get_single_region_nbor("head", "linker", self.cg), self.cg
            )

            # Assign P-O-(C)-C
            self.aa.nodes[aa_target]["map"] = bead_map + phosphate_bead
            for n in get_deadend_nbors(aa_target, self.aa):
                self.aa.nodes[n]["map"] = bead_map + phosphate_bead

            # Assign the neighboring atoms
            for n in self.aa.neighbors(aa_target):
                # Assign P-(O)-C-C
                if self.aa.nodes[n]["element"] == "O":
                    self.aa.nodes[n]["map"] = (3 * phosphate_bead) + bead_map

                # Assign P-O-C-(C)
                elif self.aa.nodes[n]["element"] == "C":
                    for i in [n] + get_deadend_nbors(n, self.aa):
                        self.aa.nodes[i]["map"] = (3 * bead_map) + phosphate_bead

            # The remaining atoms in the head group are 1:1 mapped to the X bead
            for n in self.aa_regions["head"]:
                if self.aa.nodes[n].get("map", []) == []:
                    self.aa.nodes[n]["map"] = bead_map

    def allocate_head_pi(self):
        # Get the CH atom
        aa_target = get_path_termini(
            "element",
            ["O", "P", "O", "C"],
            get_single_region_nbor("linker", "head", self.aa),
            self.aa,
        )

        # Get the C1, C2, and C3 beads
        # NOTE - the CG element typing on PI lipids is messed up because the head group beads are P1 beads
        # Maybe just do this based on bead names?
        self.bead_c1 = select_bead_by_name("C1", self.cg)
        self.bead_c2 = select_bead_by_name("C2", self.cg)
        self.bead_c3 = select_bead_by_name("C3", self.cg)

        # Assign single-bead mappings
        # CH
        self.aa.nodes[aa_target]["map"] = self.bead_c1

        # Assign the C2 and C3 parts of the ring - these are symmetrical, so use a path
        # to identify the para CHOH and then map its neighbors
        for n_i, n in enumerate(
            [
                nbor
                for nbor in self.aa.neighbors(
                    get_path_termini("element", ["C", "C", "C"], aa_target, self.aa)
                )
                if self.aa.nodes[nbor]["element"] == "C"
            ]
        ):
            self.aa.nodes[n]["map"] = [self.bead_c2, self.bead_c3][n_i]

        # Assign the intermediate mappings - any atom without a mapping takes a combination of its neighbor's mappings
        atoms_to_map = self.aa_regions["head"].copy()
        # Remove the phosphoether oxygen connecting to the linker region
        atoms_to_map.remove(get_single_region_nbor("head", "linker", self.aa))

        while (
            np.asarray(
                [self.aa.nodes[x].get("map", None) for x in atoms_to_map],
                dtype=object,
            )
            == None
        ).any():
            for n in self.aa_regions["head"]:
                if self.aa.nodes[n].get("map", None) == None:
                    # Assign a map combined from neighbor mappings
                    nbor_maps = []
                    for i in self.aa.neighbors(n):
                        if not (self.aa.nodes[i].get("map", None) == None):
                            nbor_maps.extend(self.aa.nodes[i].get("map"))
                    # If all the neighbors are un-mapped, move to a different atom
                    if len(nbor_maps) == 0:
                        continue
                    else:
                        self.aa.nodes[n]["map"] = list(set(nbor_maps))

    def allocate_head_pa(self):
        # All the head group atoms are the PO4 bead in PA lipids
        for n in self.aa_regions["head"]:
            if self.aa.nodes[n].get("map", []) == []:
                self.aa.nodes[n]["map"] = [
                    get_single_region_nbor("head", "linker", self.cg)
                ]

    def allocate_head_hexcer(self):
        # Map the hexose of the hexosylceramides.
        # Start with 1:1 bead maps
        # C11 is bead C3
        aa_target = get_single_region_nbor("head", "linker", self.aa)
        self.aa.nodes[aa_target]["map"] = self.bead_c3

        # C8 is bead C1
        aa_target = get_path_termini(
            "element",
            ["C", "C"],
            get_single_region_nbor("head", "linker", self.aa),
            self.aa,
        )
        self.aa.nodes[aa_target]["map"] = self.bead_c1

        # C9 is bead C2
        aa_target = get_path_termini(
            "element",
            ["C", "C", "C", "C"],
            get_single_region_nbor("head", "linker", self.aa),
            self.aa,
        )
        self.aa.nodes[aa_target]["map"] = self.bead_c2

        # Interpolate mappings for in-between atoms
        # Any atom without a mapping takes a combination of its neighbor's mappings
        atoms_to_map = self.aa_regions["head"].copy()

        while (
            np.asarray(
                [self.aa.nodes[x].get("map", None) for x in atoms_to_map],
                dtype=object,
            )
            == None
        ).any():
            for n in self.aa_regions["head"]:
                if self.aa.nodes[n].get("map", None) == None:
                    # Assign a map combined from neighbor mappings
                    nbor_maps = []
                    for i in self.aa.neighbors(n):
                        if not (self.aa.nodes[i].get("map", None) == None):
                            nbor_maps.extend(self.aa.nodes[i].get("map"))
                    # If all the neighbors are un-mapped, move to a different atom
                    if len(nbor_maps) == 0:
                        continue
                    else:
                        self.aa.nodes[n]["map"] = list(set(nbor_maps))

    # Allocate atoms to beads within the linker region
    def allocate_linker(self):
        # Sphingolipids have a different linker region structure so handle these separately
        if self.lipid_class in (
            LipidClass.SM,
            LipidClass.Cer,
            LipidClass.HexCer,
            LipidClass.DeoxyCer,
        ):
            if self.lipid_class == LipidClass.SM:
                self.allocate_linker_sm()
            elif self.lipid_class in (
                LipidClass.Cer,
                LipidClass.HexCer,
                LipidClass.DeoxyCer,
            ):
                self.allocate_linker_cer()
            # Fill in assignments for dead-end neighbors
            # This catches the amide carbonyl oxygen, in both sphingomyelin and ceramide.
            for n in self.aa_regions["linker"]:
                if self.aa.nodes[n].get("map", []) == []:
                    continue
                else:
                    for nbor in get_deadend_nbors(n, self.aa):
                        self.aa.nodes[nbor]["map"] = self.aa.nodes[n]["map"]
        else:
            # Non-sphingolipids
            # The phosphoether linking C1 and P is 2/3 PO4 and 1/3 GL1
            aa_target = get_single_region_nbor("head", "linker", self.aa)
            self.aa.nodes[aa_target]["map"] = [
                get_single_region_nbor("head", "linker", self.cg),
                get_single_region_nbor("head", "linker", self.cg),
                get_single_region_nbor("linker", "head", self.cg),
            ]

            # C1 is 2/3 GL1, 1/3 PO4
            aa_target = get_single_region_nbor("linker", "head", self.aa)

            bead_map = [
                get_single_region_nbor("linker", "head", self.cg),
                get_single_region_nbor("linker", "head", self.cg),
                get_single_region_nbor("head", "linker", self.cg),
            ]

            self.aa.nodes[aa_target]["map"] = bead_map

            for n in get_deadend_nbors(aa_target, self.aa):
                self.aa.nodes[n]["map"] = bead_map

            # C2 is 2/3 GL1, 1/3 GL2
            aa_target = get_path_termini(
                "element",
                ["C"],
                get_single_region_nbor("linker", "head", self.aa),
                self.aa,
            )

            bead_map = [
                get_single_region_nbor("linker", "head", self.cg),
                get_single_region_nbor("linker", "head", self.cg),
                get_single_region_nbor("linker", "tailB", self.cg),
            ]

            self.aa.nodes[aa_target]["map"] = bead_map

            for n in get_deadend_nbors(aa_target, self.aa):
                self.aa.nodes[n]["map"] = bead_map

            # C3 is GL2 GL2 GL2 PO4
            aa_target = get_path_termini(
                "element",
                ["C", "C"],
                get_single_region_nbor("linker", "head", self.aa),
                self.aa,
            )

            bead_map = [
                get_single_region_nbor("linker", "tailB", self.cg),
                get_single_region_nbor("linker", "tailB", self.cg),
                get_single_region_nbor("linker", "tailB", self.cg),
                get_single_region_nbor("head", "linker", self.cg),
            ]

            self.aa.nodes[aa_target]["map"] = bead_map

            for n in get_deadend_nbors(aa_target, self.aa):
                self.aa.nodes[n]["map"] = bead_map

            # O21 is GL1 GL1 GL2 C1A
            # Add in the tail bead assignments later
            aa_target = get_path_termini(
                "element",
                ["C", "O"],
                get_single_region_nbor("linker", "head", self.aa),
                self.aa,
            )

            bead_map = [
                get_single_region_nbor("linker", "tailA", self.cg),
                get_single_region_nbor("linker", "tailA", self.cg),
                get_single_region_nbor("linker", "tailB", self.cg),
            ]

            self.aa.nodes[aa_target]["map"] = bead_map

            # C21 is GL1 C1A
            self.aa.nodes[get_single_region_nbor("linker", "tailA", self.aa)]["map"] = (
                get_region_nbors("linker", "tailA", self.cg)
            )

            # O22 is GL1
            # NOTE this atom isn't present in ether-linked lipids. If it doesn't exist,
            # need to extend the carbon mapping to the attached hydrogens
            # to next head group atom - or can we do this in the tail assign?
            aa_target = get_path_termini(
                "element",
                ["C", "O", "C", "O"],
                get_single_region_nbor("linker", "head", self.aa),
                self.aa,
            )
            if not (aa_target == set()):
                bead_map = get_region_nbors("linker", "head", self.cg)
                self.aa.nodes[aa_target]["map"] = bead_map

            # O31 is GL2
            aa_target = get_path_termini(
                "element",
                ["C", "C", "O"],
                get_single_region_nbor("linker", "head", self.aa),
                self.aa,
            )
            self.aa.nodes[aa_target]["map"] = get_single_region_nbor(
                "linker", "tailB", self.cg
            )

            # C31 is GL2 GL1 C1B
            aa_target = get_single_region_nbor("linker", "tailB", self.aa)

            bead_map = [
                get_single_region_nbor("linker", "tailB", self.cg),
                get_single_region_nbor("linker", "tailA", self.cg),
            ]

            self.aa.nodes[aa_target]["map"] = bead_map

            # O32 is GL2
            # NOTE this atom isn't present in ether-linked lipids.  If it doesn't exist, move
            # to next head group atom
            aa_target = get_path_termini(
                "element",
                ["C", "C", "O", "C", "O"],
                get_single_region_nbor("linker", "head", self.aa),
                self.aa,
            )
            if not (aa_target == set()):
                cg_target = get_single_region_nbor("linker", "tailB", self.cg)
                self.aa.nodes[aa_target]["map"] = cg_target

    def allocate_linker_sm(self):
        # The phosphoether O1 is 2/3 PO4 and 1/3 AM1
        aa_target = get_single_region_nbor("head", "linker", self.aa)
        self.aa.nodes[aa_target]["map"] = [
            get_single_region_nbor("head", "linker", self.cg),
            get_single_region_nbor("head", "linker", self.cg),
            get_single_region_nbor("linker", "head", self.cg),
        ]

        # C1 is AM1 / PO4
        aa_target = get_single_region_nbor("linker", "head", self.aa)
        self.aa.nodes[aa_target]["map"] = [
            get_single_region_nbor("linker", "head", self.cg),
            get_single_region_nbor("head", "linker", self.cg),
        ]

        # C2 is AM1 / AM2
        aa_target = get_path_termini(
            "element",
            ["C"],
            get_single_region_nbor("linker", "head", self.aa),
            self.aa,
        )
        self.aa.nodes[aa_target]["map"] = [
            get_single_region_nbor("linker", "head", self.cg),
            get_single_region_nbor("linker", "tailB", self.cg),
        ]

        # N2 is AM2
        aa_target = get_path_termini(
            "element",
            ["C", "N"],
            get_single_region_nbor("linker", "head", self.aa),
            self.aa,
        )
        self.aa.nodes[aa_target]["map"] = [
            get_single_region_nbor("linker", "tailB", self.cg)
        ]

        # The amide carbonyl is AM2 / C1A
        # Add in the tail bead assignment later
        aa_target = get_path_termini(
            "element",
            ["C", "N", "C"],
            get_single_region_nbor("linker", "head", self.aa),
            self.aa,
        )
        self.aa.nodes[aa_target]["map"] = [
            get_single_region_nbor("linker", "tailB", self.cg)
        ]

        # C3 is AM1
        aa_target = get_path_termini(
            "element",
            ["C", "C"],
            get_single_region_nbor("linker", "head", self.aa),
            self.aa,
        )
        self.aa.nodes[aa_target]["map"] = [
            get_single_region_nbor("linker", "tailA", self.cg)
        ]

        # O3 is also AM1
        aa_target = get_path_termini(
            "element",
            ["C", "C", "O"],
            get_single_region_nbor("linker", "head", self.aa),
            self.aa,
        )
        self.aa.nodes[aa_target]["map"] = [
            get_single_region_nbor("linker", "tailA", self.cg)
        ]

    def allocate_linker_cer(self):
        # Ceramides need to be assigned from the linker region nitrogen atom rather than the head group
        source = get_path_termini(
            "element", ["C", "H"], self.ceramide_nitrogen, self.aa
        )

        bead_am1 = select_bead_by_name("AM1", self.cg)
        bead_am2 = select_bead_by_name("AM2", self.cg)

        # C1, and C3 are AM1
        aa_target = get_path_termini("element", ["C", "C"], source, self.aa)
        if (len(aa_target) != 2) or (type(aa_target) != set):
            logger.error(
                "Error! This lipid looks like a ceramide, but couldn't find two linker region carbons. Stopping."
            )
        for target in aa_target:
            self.aa.nodes[target]["map"] = bead_am1

        # O1 and O3 are AM1 / AM2
        # If there is no O1, (e.g., --C1 lipids), just map the single oxygen atom
        aa_target = get_path_termini("element", ["C", "C", "O"], source, self.aa)

        if type(aa_target) == str:
            # If we got a string back from get_path_termini, we have only a single number
            self.aa.nodes[aa_target]["map"] = bead_am1 + bead_am2
        elif (type(aa_target) == set) and (len(aa_target) != 2):
            logger.error(
                "Error! This lipid looks like a ceramide, but couldn't find two linker region oxygens. Stopping."
            )
        else:
            for target in aa_target:
                self.aa.nodes[target]["map"] = bead_am1 + bead_am2

        # C2 is 2/3 AM1, 1/3 AM2
        aa_target = get_path_termini("element", ["C"], source, self.aa)
        self.aa.nodes[aa_target]["map"] = bead_am1 * 2 + bead_am2

        # N2 is AM2
        self.aa.nodes[self.ceramide_nitrogen]["map"] = bead_am2

        # Amide carbonyl is also AM2
        aa_target = get_path_termini("element", ["C", "N", "C"], source, self.aa)
        self.aa.nodes[aa_target]["map"] = bead_am2

    # Allocate atoms to beads within the tail regions
    def allocate_both_tails(self):
        self.allocate_tail("tailA")
        self.allocate_tail("tailB")

    def allocate_tail(self, tail_label):
        frac_pos = get_fractional_positions(tail_label, self.aa)

        # Scale fractional position by the number of beads
        scaled_pos = {}
        bead_lookup = get_tail_bead_lookup(tail_label, self.cg)
        num_beads = len(bead_lookup)

        for n in frac_pos.keys():
            scaled_pos[n] = frac_pos[n] * (num_beads - 1)

        # Assign bead by integer, and remainder gives weighting to next bead
        bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]

        for n in scaled_pos.keys():
            # Integer gives starting bead
            pos = np.modf(scaled_pos[n])
            bead_start = bead_lookup[int(pos[-1])]

            # Remainder gives end bead influence
            # Catch final bead case - the last position has no influence from a next bead
            # For this case, the integer portion is the index of the final bead
            if int(pos[-1]) == (num_beads - 1):
                bead_end = ""
            else:
                bead_end = bead_lookup[int(pos[-1]) + 1]

            end_influence = np.digitize(pos[0], bins) - 1
            start_influence = (len(bins) - 1) - end_influence

            bead_map = start_influence * [bead_start] + end_influence * [bead_end]

            # If all the beads are the same reduce to a single mapped bead
            if len(np.unique(bead_map)) == 1:
                bead_map = [bead_map[0]]

            # Add to existing bead map since the first carbon has mappings from the linker region
            self.aa.nodes[n]["map"] = self.aa.nodes[n].get("map", []) + bead_map

        # For each carbon of the tails, assign the same bead map to its hydrogens
        for n in self.aa_regions[tail_label]:
            # Move on if this atom isn't a carbon - since we can't try to assign
            # hydrogens before their carbons have mappings
            if self.aa.nodes[n]["element"] == "C":
                for i in [n] + get_deadend_nbors(n, self.aa):
                    if self.aa.nodes[i]["element"] == "H":
                        # print(i,n)
                        self.aa.nodes[i]["map"] = self.aa.nodes[n]["map"]

    # Final catch for unmapped atoms
    def catch_unmapped_atoms(self):
        # Extend the neighbour's map to the un-mapped atom
        for n in self.aa.nodes:
            # Previous version - only for hydrogens
            # if (self.aa.nodes[n].get('map', []) == []) & (self.aa.nodes[n]['element'] == 'H'):
            if self.aa.nodes[n].get("map", []) == []:
                # print(n)
                # Get the mapping of the attached atom
                mapping = self.aa.nodes[list(self.aa.neighbors(n))[0]]["map"]
                self.aa.nodes[n]["map"] = mapping

    # Write the mapping file
    def write_mapping_file(self):
        self.write_header()
        self.write_atom_to_bead_map()
        self.write_head_atom_spacing()
        self.write_linker_atom_spacing()
        self.write_tail_atom_spacing("tailA")
        self.write_tail_atom_spacing("tailB")
        self.write_chiral_centers()
        self.write_double_bond_stereo()
        self.write_amide_trans()

        # Write to disk
        with open(self.outfile, "w") as f:
            f.write("\n".join(self.outfile_lines))

    def write_header(self):
        header_lines = [
            "; Created by backmapping-pipeline at: "
            + datetime.now().strftime("%d-%m-%Y %H:%M:%S"),
            "",
            "[ molecule ]",
            self.basename,
            "",
            "[ martini ]",
            self.cg_bead_order,
            "[ mapping ]",
            "charmm36",
            "",
            "[ atoms ]",
        ]
        self.outfile_lines.extend(header_lines)

    def write_atom_to_bead_map(self):
        atom_to_bead_map_lines = []
        for n in self.aa:
            line = (
                "{:>5}".format(int(n))
                + "{:>6}".format(self.aa.nodes[n]["name"])
                + "   "
                + " ".join(
                    [self.cg.nodes[x]["name"] for x in self.aa.nodes[n].get("map", [])]
                )
            )
            atom_to_bead_map_lines.append(line)
        self.outfile_lines.extend(atom_to_bead_map_lines)
        self.outfile_lines.append("")

    def write_head_atom_spacing(self):
        spacing_lines = []
        # For lipids which have a head group region,
        # begin spacing the head group from the atom
        # which borders the linker region. In glycerophospholipids,
        # this is the phosphate atom.
        if "head" in self.aa_regions.keys():
            source = get_single_region_nbor("head", "linker", self.aa)
            nbors = [get_single_region_nbor("linker", "head", self.aa)] + list(
                self.aa.neighbors(get_single_region_nbor("linker", "head", self.aa))
            )
            nbors.remove(source)

            spacing_lines.append("[ out ]")

            spacing_lines.append(
                " ".join([self.aa.nodes[n]["name"] for n in [source] + nbors])
            )

            self.seen_nodes = [source] + nbors

            all_paths = []
            for target in self.aa_regions["head"]:
                all_paths = all_paths + [
                    path for path in nx.all_simple_paths(self.aa, source, target)
                ]

            all_paths.sort(key=len)

            # The PI and HexCer sugar ring substituents are spaced out using
            # [out] tags only, while other lipid class head groups are
            # placed first with [tetra] tags, followed by [out] tags

            if self.lipid_class in (LipidClass.PI, LipidClass.HexCer):
                # Use [out] tags to space the ring substituents out
                for path in all_paths:
                    for n_i, n in enumerate(path):
                        if n in self.seen_nodes:
                            continue
                        else:
                            # If this atom is a ring carbon, don't modify its position
                            if self.aa.nodes[n]["element"] == "C":
                                continue
                            else:
                                spacing = [
                                    self.aa.nodes[n]["name"]
                                    for n in get_atom_spacing(n, path[n_i - 1], self.aa)
                                ]
                                tag = "[ out ]"
                                spacing_lines.append(tag)
                                spacing_lines.append(" ".join(spacing))
                                self.seen_nodes.append(n)
            else:
                for path in all_paths:
                    for n_i, n in enumerate(path):
                        if n in self.seen_nodes:
                            continue
                        else:
                            # Define the atoms needed for a tetrahedral placement, namely:
                            # 1. The anchor atom - becomes a vertex of the tetrahedron
                            # 2. The central atom (won't be moved) - this is n
                            # 3. Up to three additional atoms joined to the center to be placed on vertices of the tetrahedron
                            anchor = path[n_i - 1]
                            center = n

                            # Remove the anchor from the list of neighbors to place
                            nbors_to_place = [
                                nbor
                                for nbor in self.aa.neighbors(n)
                                if not (nbor == anchor)
                            ]

                            # Depending on the number of atoms to be placed, we can use an out or tetra tag
                            # If the center has no non-anchor neighbors, use an out tag
                            if len(nbors_to_place) == 0:
                                spacing = [
                                    self.aa.nodes[n]["name"]
                                    for n in get_atom_spacing(n, path[n_i - 1], self.aa)
                                ]
                                tag = "[ out ]"

                            else:
                                # If we have neighbors that are not the anchor, we need a tetrahedral placement
                                spacing_atoms = [anchor, center] + nbors_to_place
                                spacing = [
                                    self.aa.nodes[i]["name"] for i in spacing_atoms
                                ]
                                tag = "[ tetra ]"

                            # Adding just the center to the record of seen nodes means that
                            # dead-end atoms are first placed using a tetra flag, and then
                            # subsequently re-placed using an out tag.
                            # This might be beneficial, but let's see...
                            spacing_lines.append(tag)
                            spacing_lines.append(" ".join(spacing))
                            self.seen_nodes.append(n)
        self.outfile_lines.extend(spacing_lines)

    def write_linker_atom_spacing(self):
        spacing_lines = []
        spacing_lines.append("[ out ]")
        # Ceramide lipids do not have a head group region, so set the source
        # for tracing paths throughout the linker region to the nitrogen atom
        # For other lipids, set the source to the head region atom
        # bordering the linker region.
        if self.lipid_class in (LipidClass.Cer, LipidClass.DeoxyCer):
            source = self.ceramide_nitrogen
            self.seen_nodes = [self.ceramide_nitrogen]
        else:
            # Position the region neighbors first, and don't move these later on
            # Head group border
            spacing = [
                self.aa.nodes[n]["name"]
                for n in get_atom_spacing(
                    get_single_region_nbor("linker", "head", self.aa),
                    get_single_region_nbor("head", "linker", self.aa),
                    self.aa,
                )
            ]
            self.seen_nodes.append(get_single_region_nbor("linker", "head", self.aa))
            spacing_lines.append(" ".join(spacing))

            # Start traversing bonded paths within the linker region from the head region nbor
            # The head region has been spaced out at this point, so let's start the paths
            # from the already-placed head region nbor.
            source = get_single_region_nbor("head", "linker", self.aa)

        # Write spacing lines for the linker region atoms bordering the tails
        for tail in ["tailA", "tailB"]:
            spacing = [
                self.aa.nodes[n]["name"]
                for n in get_atom_spacing(
                    get_single_region_nbor("linker", tail, self.aa),
                    get_single_region_nbor(tail, "linker", self.aa),
                    self.aa,
                )
            ]
            self.seen_nodes.append(get_single_region_nbor("linker", tail, self.aa))
            spacing_lines.append(" ".join(spacing))

        # Define paths to all linker region atoms from the source selected above
        all_paths = []
        for target in self.aa_regions["linker"]:
            all_paths = all_paths + [
                path for path in nx.all_simple_paths(self.aa, source, target)
            ]

        all_paths.sort(key=len)

        for path in all_paths:
            for n_i, n in enumerate(path):
                if n in self.seen_nodes:
                    continue
                else:
                    # Define the atoms needed for a tetrahedral placement, namely:
                    # 1. The anchor atom - becomes a vertex of the tetrahedron
                    # 2. The central atom (won't be moved) - this is n
                    # 3. Up to three additional atoms joined to the center to be placed on vertices of the tetrahedron
                    anchor = path[n_i - 1]
                    center = n

                    # Remove the anchor from the list of neighbors to place
                    # We also need to avoid re-placing the region border atoms placed
                    # above. So, exclude these atoms from the nbors
                    nbors_to_place = [
                        nbor
                        for nbor in self.aa.neighbors(n)
                        if not ((nbor == anchor) or (nbor in self.seen_nodes))
                    ]

                    # Depending on the number of atoms to be placed, we can use an out or tetra tag
                    # If the center has no non-anchor neighbors, use an out tag
                    if len(nbors_to_place) == 0:
                        spacing = [
                            self.aa.nodes[n]["name"]
                            for n in get_atom_spacing(n, path[n_i - 1], self.aa)
                        ]
                        tag = "[ out ]"

                    else:
                        # If we have neighbors that are not the anchor, we need a tetrahedral placement

                        spacing_atoms = [anchor, center] + nbors_to_place
                        spacing = [self.aa.nodes[i]["name"] for i in spacing_atoms]
                        tag = "[ tetra ]"

                    # Adding just the center to the record of seen nodes means that
                    # dead-end atoms are first placed using a tetra flag, and then
                    # subsequently re-placed using an out tag.
                    # This might be beneficial, but let's see...
                    spacing_lines.append(tag)
                    spacing_lines.append(" ".join(spacing))
                    self.seen_nodes.append(n)
        self.outfile_lines.extend(spacing_lines)

    def write_tail_atom_spacing(self, tail_label: str):
        spacing_lines = []
        source = get_single_region_nbor(tail_label, "linker", self.aa)
        self.seen_nodes = [source]  # TODO: Clobbering - check this is OK

        # Add all tail carbon atoms to the list of seen atoms
        for n in self.aa_regions[tail_label]:
            if self.aa.nodes[n]["element"] == "C":
                self.seen_nodes.append(n)

        all_paths = []
        for target in self.aa_regions[tail_label]:
            all_paths = all_paths + [
                path for path in nx.all_simple_paths(self.aa, source, target)
            ]

        all_paths.sort(key=len)

        for path in all_paths:
            for n_i, n in enumerate(path):
                if n in self.seen_nodes:
                    continue
                else:
                    anchor = path[n_i - 1]
                    center = n
                    nbors_to_place = [
                        nbor
                        for nbor in self.aa.neighbors(n)
                        if not ((nbor == anchor) or (nbor in self.seen_nodes))
                    ]

                    # Depending on the number of atoms to be placed, we can use an out or tetra tag
                    # If the center has no non-anchor neighbors, use an out tag
                    if len(nbors_to_place) == 0:
                        spacing = [
                            self.aa.nodes[n]["name"]
                            for n in get_atom_spacing(n, path[n_i - 1], self.aa)
                        ]
                        tag = "[ out ]"

                    else:
                        # If we have neighbors that are not the anchor, we need a tetrahedral placement
                        spacing_atoms = [anchor, center] + nbors_to_place
                        spacing = [self.aa.nodes[i]["name"] for i in spacing_atoms]
                        tag = "[ tetra ]"

                    # Adding just the center to the record of seen nodes means that
                    # dead-end atoms are first placed using a tetra flag, and then
                    # subsequently re-placed using an out tag.
                    # This might be beneficial, but let's see...
                    spacing_lines.append(tag)
                    spacing_lines.append(" ".join(spacing))
                    self.seen_nodes.append(n)
        self.outfile_lines.extend(spacing_lines)

    def write_chiral_centers(self):
        spacing_lines = ["[ chiral ]"]
        # Check if we have a chiral center defined in the stereo reference
        # Most lipids have at least one
        if not "chiral" in self.stereo[self.basename].keys():
            logger.debug("No chiral centers found in stereo reference")
            return

        chiral_centers = self.stereo[self.basename]["chiral"].keys()
        logger.debug(f"Found {len(chiral_centers)} chiral centers in stereo reference")
        for chiral_center in chiral_centers:
            chiral_center_node = select_bead_by_name(chiral_center.strip(), self.aa)[0]
            nbors = self.aa.neighbors(chiral_center_node)
            heavy_nbors = []
            for nbor in nbors:
                # Store the names of the hydrogen atom to be placed and the heavy atoms
                # Strip whitespace for agreement with the ITP-derived names
                if self.aa.nodes[nbor]["element"] == "H":
                    hydrogen = self.aa.nodes[nbor]["name"].strip()
                else:
                    heavy_nbors.append(self.aa.nodes[nbor]["name"].strip())
            spacing_line = " ".join([hydrogen, chiral_center.strip()] + heavy_nbors)
            spacing_lines.append(spacing_line)
        self.outfile_lines.extend(spacing_lines)

    def write_double_bond_stereo(self):
        spacing_lines = []
        # Check if we have a double bond defined in the stereo reference
        if "unsat" in self.stereo[self.basename].keys():
            logger.debug("Found double bond in stereo reference")
        else:
            logger.debug("No double bond found in stereo reference")
            return

        # Check if the unsaturation has defined stereochemistry
        for unsat_bond in self.stereo[self.basename]["unsat"].values():
            if unsat_bond["stereo"] == "NONE":
                logger.debug(
                    f"No stereochemistry defined for double bond: {unsat_bond}"
                )
                continue
            elif unsat_bond["stereo"] == "E":
                spacing_lines.append("[ trans ]")
            elif unsat_bond["stereo"] == "Z":
                spacing_lines.append("[ cis ]")

            # Get the atom names to write spacing
            begin = unsat_bond["begin"]
            end = unsat_bond["end"]

            # Also need to get the names for the attached hydrogen atoms
            # Strip whitespace for agreement with the ITP-derived names
            hydrogen_names = []
            for carbon_name in [begin, end]:
                found_atoms = select_bead_by_name(carbon_name.strip(), self.aa)
                if len(found_atoms) != 1:
                    # Should only ever match a single atom
                    msg = f"Error! Found {len(found_atoms)} atoms matching {carbon_name.strip()}"
                    logger.error(msg)
                    raise ValueError(msg)
                hydrogen_name = get_deadend_nbors(found_atoms[0], self.aa)[0]
                # Double bonds by definition have only a single hydrogen at each carbon
                # So there will only be a single deadend neighbor
                hydrogen_names.append(hydrogen_name)

            # Write spacing line for the double bond in both the forward and backward direction
            # This is to space out both of the hydrogen atoms
            bond_names = [
                self.aa.nodes[hydrogen_names[0]]["name"],
                begin.strip(),
                end.strip(),
                self.aa.nodes[hydrogen_names[1]]["name"],
            ]
            spacing_lines.append(" ".join(bond_names))
            bond_names.reverse()
            spacing_lines.append(" ".join(bond_names))
        self.outfile_lines.extend(spacing_lines)

    def write_amide_trans(self):
        spacing_lines = []
        # An amide entry is always defined, but may not be populated
        # Check if we have an amide bond defined in the stereo reference
        if self.stereo[self.basename]["amide"] == {}:
            logger.debug("No amide bond found in stereo reference")
            return
        else:
            logger.debug("Found amide bond/s in stereo reference")

        # Get atom names for the amide atoms
        spacing_lines.append("[ trans ]")
        for amide_bond in self.stereo[self.basename]["amide"].values():
            h = amide_bond["H"]
            n = amide_bond["N"]
            c = amide_bond["C"]
            o = amide_bond["O"]
            # Compensate for the fact the ITP file is 1-indexed but the stereo reference is 0-indexed
            spacing_lines.append(
                " ".join(
                    [self.aa.nodes[str(atom + 1)]["name"] for atom in [h, n, c, o]]
                )
            )
        self.outfile_lines.extend(spacing_lines)
