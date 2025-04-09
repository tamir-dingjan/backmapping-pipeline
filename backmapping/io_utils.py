from rdkit import Chem
import networkx as nx
import fnmatch


def load_structure_from_file(path: str):
    structure = Chem.rdmolfiles.MolFromPDBFile(path, sanitize=True, removeHs=False)
    return structure


def load_itp_as_network(filepath):
    with open(filepath, "r") as itp:
        # Construct a new network
        topology = nx.Graph()

        current_section = None
        heteroflag = None

        for line in itp.readlines():
            if (line[0] == ";") or (line[0] == "#"):
                continue

            elif len(line.split()) == 0:
                # Blank lines mark the end of a section
                current_section = None

            elif current_section == "atoms":
                # Add this atom as a network node
                atom = line.split()
                # Define heteroatom flag based on the first character of the atomname
                if (atom[1][0] == "H") or (atom[1][0] == "C"):
                    heteroflag = False
                else:
                    heteroflag = True

                topology.add_node(
                    atom[0],
                    element=atom[1][0],
                    atomtype=atom[1],
                    name=atom[4],
                    hetero=heteroflag,
                )

            elif current_section == "bonds":
                # Add this bond as an edge
                bond = line.split()
                topology.add_edge(bond[0], bond[1])

            elif fnmatch.fnmatch(line, "*[[]*atoms*[]]*"):
                current_section = "atoms"

            elif fnmatch.fnmatch(line, "*[[]*bonds*[]]*"):
                current_section = "bonds"

    return topology
