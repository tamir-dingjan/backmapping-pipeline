import numpy as np
import networkx as nx
import fnmatch
from collections import defaultdict


def has_heteroatom_neighbors(network, node):
    """
    Method to check if the supplied node has heteroatom neighbors in the
    given network.

    The neighbors are heteroatoms if they are anything other than carbon
    or hydrogen. The element is determined by the first character of the
    atom name during ITP loading.

    Parameters
    ----------
    network : networkx network
        The topology network
    node : str
        The integer ID of the node being checked

    Returns
    -------
    bool
        True if the node has a heteroatom neighbor, otherwise False.

    """

    return bool(sum([network.nodes[n]["hetero"] for n in network.neighbors(node)]))


def check_prelinker_atom(network, n):
    """
    This check detects if there are remaining atoms that should be considered
    linker atoms - parts of glycerol or the fatty acyl carbonyl.

    The detection is based on atoms being labelled unknown with linker neighbors.

    Parameters
    ----------
    network : networkx network

    n : node ID


    Returns
    -------
    bool


    """
    if (network.nodes[n]["region"] == "unknown") and (
        "linker" in [network.nodes[nbor]["region"] for nbor in network.neighbors(n)]
    ):
        return True
    else:
        return False


def check_unassigned_linker_atoms(network):
    for n in network.nodes:
        if check_prelinker_atom(network, n):
            return True
    else:
        return False


def check_generic_tail(network, n):
    if network.nodes[n]["region"] == "generic_tail":
        return True
    else:
        return False


def check_unassigned_generic_tail(network):
    for n in network.nodes:
        if check_generic_tail(network, n):
            return True
    else:
        return False


def get_deadend_nbors(node, network):
    """
    Returns a list of neighbors for the given node which are only connected to that node.

    Parameters
    ----------
    node : str
        node ID.
    network : networkx network
        Bond topology network.

    Returns
    -------
    list
        List of neighboring nodes which are dead-end neighbors.

    """

    return [
        nbor
        for nbor in network.neighbors(node)
        if list(network.neighbors(nbor)) == [node]
    ]


def get_path_termini(attribute, path, source, network, exclusions=None):
    """
    Find the terminal nodes that match a path of specified attribute values.


    Parameters
    ----------
    attribute : str
        The node attribute specifying the path value.
    path : list
        The attribute values that define the path.
    source : str
        Node ID to begin path.
    network : networkx network
        Topology network.
    exclusions: list
        List of nodes to exclude from matching.

    Returns
    -------
    str, set
        The set of node IDs that match the given attribute path from the source.
        If there is only a single matching ID, returns a string.

    """
    # termini = [] # list approach
    termini = set()  # set approach

    if exclusions == None:
        exclusions = []

    # This is how the termini is handed back up the recursion tree
    if len(path) == 0:
        return source

    else:
        # Walk a step forward on the path
        # 1. Shorten the path by one step
        # 2. Add the source to the list of nodes seen
        # 3. Find new termini starting from the matching neighbor
        # print(source, list(network.neighbors(source)), exclusions)
        for nbor in network.neighbors(source):
            if (network.nodes[nbor][attribute] == path[0]) and (
                set(exclusions).isdisjoint({nbor})
            ):
                shorter_path = path.copy()
                shorter_path.pop(0)
                exclusions.append(source)
                # termini = termini + get_path_termini(attribute, shorter_path, nbor, network, exclusions)
                result = get_path_termini(
                    attribute, shorter_path, nbor, network, exclusions
                )
                # print("result:", result)
                # if result == []: # list approach
                if result == set():  # set approach
                    continue
                else:
                    # termini.extend(result) # list approach

                    # set approach
                    if type(result) == set:
                        termini.update(result)
                    else:
                        termini.add(result)
                    # set approach

    # Reduce list to unique termini
    # Note, if termini is a string, it might be a multi-digit integer.
    # Directly casting that as a set will break the digits apart.
    # So, we need to catch this case.
    # Try to catch this by using termini as a set and adding new results with set.add(),
    # which preserves multi-character strings

    # termini = list(set(termini))
    # print("termini:", termini, type(termini))
    if len(termini) == 1:
        return list(termini)[0]
    else:
        return termini


def write_with_newline(line, file):
    file.write(line + "\n")


def get_single_region_nbor(source_region, nbor_region, network):
    """
    Picks a single node from the source region which neighbors a node in the
    specified neighbor region.

    This routine returns the first node that satisfies the criteria,
    so do not use if you have regions which share more than one boundary edge.

    Parameters
    ----------
    source_region : str
        Source region from which region neighbors are picked.
    nbor_region : str
        Neighboring region to use as criteria for node picking.
    network : networkx network
        The topology network.

    Returns
    -------
    nbor : str
        The node ID from source_region which neighbors nbor_region

    """

    for n in network.nodes:
        if network.nodes[n]["region"] == source_region:
            for nbor in network.neighbors(n):
                if network.nodes[nbor]["region"] == nbor_region:
                    return n


def get_fractional_positions(tail, network):
    pos = {}
    frac_pos = {}
    count = 0
    reached_tail_end = False
    seen_nodes = []

    # Beginning from the start of the tail, track all the tail carbons
    n = get_single_region_nbor("linker", tail, network)

    while not reached_tail_end:

        for nbor in network.neighbors(n):
            # Move to next node if it's a new tail carbon atom
            if (
                (network.nodes[nbor]["region"] == tail)
                and (network.nodes[nbor]["element"] == "C")
                and not (nbor in seen_nodes)
            ):
                next_node = nbor
                seen_nodes.append(n)
                pos[n] = count
                count += 1

        # If we couldn't update the next_node with a valid neighbor, that's the end of the tail
        if next_node == n:
            reached_tail_end = True
            pos[n] = count
            # count += 1
        else:
            n = next_node

    # Make positions fractional
    for node in pos.keys():
        frac_pos[node] = pos[node] / count

    return frac_pos


def get_tail_bead_lookup(tail, network):
    n = get_single_region_nbor(tail, "linker", network)
    bead_lookup = [n]
    seen_nodes = []

    # How many tail beads do we need to sort?
    tail_beads = [n for n in network if network.nodes[n]["region"] == tail]

    while len(tail_beads) > len(bead_lookup):
        for nbor in network.neighbors(n):
            if (network.nodes[nbor]["region"] == tail) and not (nbor in seen_nodes):
                seen_nodes.append(n)
                bead_lookup.append(nbor)
        n = nbor

    return bead_lookup


def get_atom_spacing(n, anchor, network):
    spacers = list(network.neighbors(anchor))
    spacers.remove(n)
    return [n, anchor] + spacers


def select_bead_by_name(name, network):
    """
    Returns a list of bead IDs which match the given name.

    Parameters
    ----------
    name : str
        The desired bead name.
    network : networkx network
        The bond topology network to search.

    Returns
    -------
    list
        Matching bead IDs.

    """
    return [n for n in network.nodes if network.nodes[n]["name"] == name]


def get_ceramide_nitrogen(network):
    """
    Returns the node ID of the first non-choline nitrogen atom found.
    Raises an exception if fewer or more than one nitrogen atom is found.

    Parameters
    ----------
    network : networkx network
        The all-atom topology network.

    Raises
    ------
    Exception
        The error raised if fewer or more than a single nitrogen atom is found.

    Returns
    -------
    str
        The node ID of the identified nitrogen atom.

    """
    nitrogens = []

    for i in network.nodes:
        if network.nodes[i]["element"] == "N":
            check_choline = sum(
                [check_methyl_carbon(network, nbor) for nbor in network.neighbors(i)]
            )
            if not check_choline:
                nitrogens.append(i)

    if len(nitrogens) != 1:
        raise Exception(
            "Error! There is more than one nitrogen atom in this ceramide! Stopping."
        )
    else:
        return nitrogens[0]


def check_methyl_carbon(network, i):
    """
    Returns True if the provided node is a carbon atom bonded to three hydrogens.

    Parameters
    ----------
    network : networkx topology
        All-atom networkx topology.
    i : str
        An atom in the network.

    Returns
    -------
    bool.
        True if a carbon with three hydrogens, otherwise False.

    """
    hydrogen_count = 0

    for nbor in network.neighbors(i):
        if network.nodes[nbor]["element"] == "H":
            hydrogen_count += 1

    if hydrogen_count == 3:
        return True
    else:
        return False


def get_region_nbors(source_region, nbor_region, network):
    """
    Picks nodes from the source region which neighbor nodes in the
    specified neighbor region.

    Parameters
    ----------
    source_region : str
        Source region from which region neighbors are picked.
    nbor_region : str
        Neighboring region to use as criteria for node picking.
    network : networkx network
        The topology network.

    Returns
    -------
    nbors : list
        The list of node IDs from source_region which neighbor
        nbor_region

    """
    nbors = []
    for n in network.nodes:
        if network.nodes[n]["region"] == source_region:
            for nbor in network.neighbors(n):
                if network.nodes[nbor]["region"] == nbor_region:
                    nbors.append(n)
    return nbors


def get_double_bonds(tail_label, network, aa_regions):
    """
    Returns a list of the double bonds for [cis] spacing.
    The bonds returned contain four atoms: the double-bonded carbons
    and the hydrogen atoms attached to them.

    Parameters
    ----------
    tail_label : TYPE
        DESCRIPTION.
    network : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    counted = []
    pairs = []

    # Get carbons with only a single deadend neighbor
    unsats = [
        n for n in aa_regions[tail_label] if (len(get_deadend_nbors(n, network)) == 1)
    ]

    # Find neighboring carbons
    for i in unsats:
        for j in unsats:
            if i == j or j in counted:
                continue
            # Are these two neighbors?
            if j in network.neighbors(i):
                counted.append(i)
                counted.append(j)

                # Get the double bond hydrogen atoms for the [cis] flag
                a = get_deadend_nbors(i, network)[0]
                d = get_deadend_nbors(j, network)[0]

                # Get the double bond neighboring carbon atoms
                # for x in network.neighbors(i):
                #     if (x == j) or (x in get_deadend_nbors(i, network)):
                #         continue
                #     else:
                #         a = x

                # for x in network.neighbors(j):
                #     if (x == i) or (x in get_deadend_nbors(j, network)):
                #         continue
                #     else:
                #         d = x

                pairs.append([a, i, j, d])

    return pairs
