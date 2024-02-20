"""
Modified version of ligand_approx.py from the graph_search_ligand_approximation repository. 
"""

import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
from collections import defaultdict
import optparse
import re

import sys
sys.path.append('/nfs/polizzi/shared/programs/chemistry/graph-search-ligand-approximation/graph_matcher/')

from constants import fragment_graph_to_amino_acid_names, names_to_smile_map
from graph_search import match_graphs, convert_atom_nameidx_to_idx, is_attached_to_pi_system


def matrix_from_smile(smile):
    """Given a smile string, loads RDkit molecule & constructs an atom adjacency matrix and lists of atom properties"""
    template = AllChem.MolFromSmiles(smile)
    template = Chem.AddHs(template)  # Comment this out to only use heavy atoms.
    AllChem.EmbedMolecule(template)
    AllChem.UFFOptimizeMolecule(template)
    template.GetConformer()

    coords = []  # Feature matrix of coordinates for each atom.
    identities = []  # Feature matrix of identities for each atom.
    formal_charges = []
    # Basic adjacency matrix for atoms, edge weights are number of bonds (1.5 == aromatic).
    adj_matrix = np.zeros((template.GetNumAtoms(), template.GetNumAtoms()))

    for x in template.GetBonds():
        bidx, eidx = (x.GetBeginAtomIdx(), x.GetEndAtomIdx())
        adj_matrix[bidx, eidx] = x.GetBondTypeAsDouble()
        adj_matrix[eidx, bidx] = x.GetBondTypeAsDouble()

    for a in range(template.GetNumAtoms()):
        positions = template.GetConformer().GetAtomPosition(a)
        coords.append([positions.x, positions.y, positions.z])
        identities.append(template.GetAtomWithIdx(a).GetSymbol())
        formal_charges.append(template.GetAtomWithIdx(a).GetFormalCharge())

    return adj_matrix, coords, identities, formal_charges


def parse_pdb_rdkit_molecule(mol):
    """Given a rdkit molecule created from a pdb file, constructs an adjacency matrix for
    graph representation as well as useful information. Unlike parse_ligand_matrix_from_mol2,
    the pdb name for each index is stored in the molecule, We need these to use the molecule in COMBS."""
    coords = []  # Feature matrix of coordinates for each atom.
    identities = []  # Feature matrix of identities for each atom.
    pdb_indices = []
    formal_charges = []
    # Basic adjacency matrix for atoms, edge weights are number of bonds (1.5 == aromatic).
    adj_matrix = np.zeros((mol.GetNumAtoms(), mol.GetNumAtoms()))

    for x in mol.GetBonds():
        bidx, eidx = (x.GetBeginAtomIdx(), x.GetEndAtomIdx())
        adj_matrix[bidx, eidx] = x.GetBondTypeAsDouble()
        adj_matrix[eidx, bidx] = x.GetBondTypeAsDouble()

    for a in range(mol.GetNumAtoms()):
        # positions = mol.GetConformer().GetAtomPosition(a)
        # coords.append([positions.x, positions.y, positions.z])
        identities.append(mol.GetAtomWithIdx(a).GetSymbol())
        pdb_indices.append(mol.GetAtomWithIdx(a).GetPDBResidueInfo().GetName())
        formal_charges.append(mol.GetAtomWithIdx(1).GetFormalCharge())

    return adj_matrix, coords, identities, formal_charges, pdb_indices


def get_3_letter_code_from_mol2(path_to_mol2):
    """Very roughly parses the 3-letter code from the TRIPOS formatted mol2 file used as an input."""
    output = ""
    can_print = False
    # Loops over all the lines in the file and gets the 3-letter code from first line after tripos substructure line.
    with open(path_to_mol2, "r") as f:
        for i in f.readlines():
            if can_print:
                output += i.strip()
            if "@<TRIPOS>SUBSTRUCTURE" in i:
                can_print = True
    last_line = output.strip().split()
    pattern = re.compile("^[A-Z0-9]{3}$")
    matches = [y[0] for y in [pattern.match(x) for x in last_line] if y is not None]

    # Takes the first Capital-Alpha-Numeric 3-letter code from line under @<TRIPOS>SUBSTRUCTURE, this might not be right.
    return matches[0]

def get_3_letter_code_from_pdb(path_to_pdb):
    lines = open(path_to_pdb, 'r').read().splitlines()
    for l in lines:
        if l.startswith('HETATM'):
            code = l.split()[3]
            break 
    return code

def correct_pdb_ligand_name(pdb_path, ligand_name):
    """Corrects the three-letter-code output by RDkit (UNL) for the actual 3-letter ligand code."""
    lines = None
    with open(pdb_path, "r") as f:
        lines = [line for line in f.readlines()]
        
    if lines:
        with open(pdb_path, "w") as f:
            for line in lines:
                f.write(line.replace("UNL", ligand_name))


def parse_ligand_matrix_from_mol2(mol2_path, sdf_path, ligand_pdb_output):
    """
    Constructs an adjacency matrix from a mol2 file and presumed sdf file with same name,
    ligand as well as writes the ligand to a pdb file to be used by COMBS downstream if necessary
    """

    # If an sdf path isn't provided, try loading from the same path as mol2 with same name.
    # This doesn't really matter unless loading from mol2 fails.
    if not sdf_path:
        sdf_path = mol2_path.replace(".mol2", ".sdf")

    # Create an output pdb path using the mol2 file name if one isn't already provided.
    if not ligand_pdb_output:
        ligand_pdb_output = mol2_path.replace(".mol2", ".pdb")

    # Find the 3-letter code for the ligand.
    ligand_name = get_3_letter_code_from_mol2(mol2_path)

    # Load the ligand from a mol2 file, if this fails try loading from sdf, quits if both fail.
    mol = Chem.MolFromMol2File(mol2_path, removeHs=False)
    if mol is None:
        print(f"Couldn't load {mol2_path} file. Trying to load from {sdf_path} file.")
        mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
        if mol is None:
            print(f"Also couldn't load from SDF file.")
            raise RuntimeError("mol2 and sdf files could not be loaded by rdkit")

    # Dump ligand to pdb file which I think can be used as the input to COMBS.
    # Flavor 4 ensures bond orders are written into the pdb, correct UNL residues in generated pdb file.
    Chem.MolToPDBFile(mol, ligand_pdb_output, flavor=4) 
    correct_pdb_ligand_name(ligand_pdb_output, ligand_name)

    # Reload the ligand from the PDB file so that the atoms retain PDB file name data used by COMBS.
    template = Chem.MolFromPDBFile(ligand_pdb_output, removeHs=False)
    # AllChem.EmbedMolecule(template)
    # AllChem.UFFOptimizeMolecule(template)
    # template.GetConformer()

    # Return the adjacency matrix and node information as well is 3-letter lig name
    adj_matrix, coords, identities, formal_charges, pdb_indices = parse_pdb_rdkit_molecule(template)
    return adj_matrix, coords, identities, formal_charges, pdb_indices, ligand_name

def parse_ligand_matrix_from_pdb(pdb_path):
    ligand_name = get_3_letter_code_from_pdb(pdb_path)
    template = Chem.MolFromPDBFile(pdb_path, removeHs=False)
    # Return the adjacency matrix and node information as well is 3-letter lig name
    adj_matrix, coords, identities, formal_charges, pdb_indices = parse_pdb_rdkit_molecule(template)
    return adj_matrix, coords, identities, formal_charges, pdb_indices, ligand_name

def fetch_ligand_smile(three_letter_code):
    """
    If you have a csv dataframe of smiles known to the RCSB in your working directory, you can look up the
    ligand smile from its 3-letter code with this function.

    Function is unused at the moment.
    """
    expo_dict = pd.read_csv('Components-smiles-stereo-oe.smi', sep="\t", header=None, names=["SMILES", "ID", "Name"])
    expo_dict.set_index("ID", inplace=True)
    sub_smiles = expo_dict["SMILES"][three_letter_code]
    return sub_smiles


def draw_graph_from_adj_matrix(mtx, name, node_identities=None, size=12):
    """Plots a network representation of a graph's adjacency matrix for inspection."""
    use_cmap = False
    node_colors = []
    if node_identities:
        color_map = {"C": "grey", "N": "b", "O": "r", "H": "lightgrey", "S": "orange", "Cl": "green", "P": "purple"}
        for i in node_identities:
            try:
                node_colors.append(color_map[i])
            except KeyError:
                node_colors.append("teal")
        use_cmap = True
    G = nx.Graph(mtx)
    pos = nx.kamada_kawai_layout(G)
    fig, ax = plt.subplots(1, 1)
    if use_cmap:
        nx.draw_networkx(G, pos=pos, ax=ax, with_labels=True, node_size=300, font_size=10, node_color=node_colors)
        nx.draw_networkx_edge_labels(G, pos=pos, ax=ax, font_size=5) # Draws edge weights (poorly) on the lines
    else:
        nx.draw_networkx(G, pos=pos, ax=ax, with_labels=True, node_size=150, font_size=5)
    plt.title(name)
    fig.set_size_inches(size, size)
    plt.show()


# TODO: move graph generation to graph_search file for better compartmentalization?
def initialize_graph(mtx, coords, labels, formal_charges, name, pdb_indices=None, plot_graph=False):
    """Creates a networkx graph from the adjacency matrix and lists of node data. Also plots the graph if asked."""
    G = nx.Graph(mtx, name=name)

    if pdb_indices:
        # for atom, coord, node, index, charge in zip(labels, coords, G.nodes(), pdb_indices, formal_charges):
        for atom, node, index, charge in zip(labels, G.nodes(), pdb_indices, formal_charges):
            G.nodes[node]["atom"] = atom
            # G.nodes[node]["xyz"] = coord
            G.nodes[node]["pdb_idx"] = index.strip()
            G.nodes[node]['charge'] = charge
    else:
        # for atom, coord, node, charge in zip(labels, coords, G.nodes(), formal_charges):
        for atom, node, charge in zip(labels, G.nodes(), formal_charges):
            G.nodes[node]["atom"] = atom
            # G.nodes[node]["xyz"] = coord
            G.nodes[node]['charge'] = charge

    for node in G.nodes():
        # Precompute these since computing them on the graphs/subgraphs we match with can lead to undefined behavior
        G.nodes[node]['atom_neighbors'] = [G.nodes[x]['atom'] for x in list(G.neighbors(node))]
        G.nodes[node]['is_attached_pi'] = is_attached_to_pi_system(f"{node}", G)

    if plot_graph:
        draw_graph_from_adj_matrix(mtx, name, labels)
    return G


def sort_mapping_for_output(match_mapping):
    """
    COMBS requires that the order of the atoms in the lignad.txt for each chemical group is consistent.
    We have to sort them before outputting since we are using sets to generate them which shuffles the order.
    """
    sorted_match_list = [(x, match_mapping[x]) for x in sorted(match_mapping.keys())]
    return sorted_match_list


def write_matches_to_ligand_txt_file(f, all_match_lists, ligand_code, combs_group, ligand_graph, cg_incidence_counter, graph_coverage, verbose):
    f.write("# " + combs_group.upper() + "\n")
    for match_list in all_match_lists:
        for match_mapping in match_list:
            if verbose:
                print(match_mapping)
            match_mapping = sort_mapping_for_output(match_mapping)
            cg_incidence_counter[combs_group] += 1
            for mapped_amino_acid in fragment_graph_to_amino_acid_names[combs_group]:
                for frag_atm, lig_atm in match_mapping:
                    # Record matched atoms for summary statistics.
                    if "!" not in lig_atm:
                        graph_coverage[lig_atm.replace("*", "").replace("~", "").replace("$", "")] += 1

                    # Write line in ligand.txt
                    f.write(
                        f"{ligand_code} {ligand_graph.nodes[convert_atom_nameidx_to_idx(lig_atm)]['pdb_idx']} "
                        f"{mapped_amino_acid} "
                        f"{fragment_graph_to_amino_acid_names[combs_group][mapped_amino_acid][frag_atm]} "
                        f"{combs_group} {cg_incidence_counter[combs_group]} 1"
                    )
                    if "!" in lig_atm:
                        f.write(" no_rmsd")
                    if "*" in lig_atm:
                        f.write(" is_not_donor")
                    if "~" in lig_atm:
                        f.write(" is_not_acceptor")
                    if "$" in lig_atm:
                        f.write(" is_acceptor")
                    f.write("\n")
            f.write("\n")
        if verbose:
            print()
        f.write("\n")


def read_lig_txt(path_to_txt):
    """
    This function is copied directly from COMBS.
    """
    dtype_dict = {'lig_resname': str, 'lig_name': str, 'resname': str, 'name': str, 'CG_type': str, 'CG_group': int, 'CG_ligand_coverage': int, 'rmsd': bool, 'is_acceptor': bool, 'is_donor': bool, 'is_not_acceptor': bool, 'is_not_donor': bool}
    data = []
    with open(path_to_txt, 'r') as infile:
        for line in infile:
            if line[0] == '#':
                continue
            try:
                spl = line.strip().split()
                if len(spl) >= 7 and len(spl) <= 9:
                    new_spl = np.zeros(12, dtype=object)
                    new_spl[:7] = spl[:7]
                    if 'no_rmsd' not in line:
                        new_spl[7] = True
                    if 'is_acceptor' in line:
                        new_spl[8] = True
                    if 'is_donor' in line:
                        new_spl[9] = True
                    if 'is_not_acceptor' in line:
                        new_spl[10] = True
                    if 'is_not_donor' in line:
                        new_spl[11] = True
                    spl.append(True)
                else:
                    continue
                data.append(new_spl)
            except Exception:
                pass
    df = pd.DataFrame(data, columns=['lig_resname', 'lig_name', 'resname', 'name', 'CG_type',
                                     'CG_group', 'CG_ligand_coverage', 'rmsd', 'is_acceptor',
                                     'is_donor', 'is_not_acceptor', 'is_not_donor'])
    return df.astype(dtype=dtype_dict)


def main(pdb_path, mol2_path, sdf_path, lig_pdb_out, lig_txt_out, write_df_pickle=False, verbose=False, plot_graphs=False, hydrophobics=False, hydrophilics=False):
    # Create matrix representation of ligand graph.
    if pdb_path: 
        ligand_data = parse_ligand_matrix_from_pdb(pdb_path)
    else:
        ligand_data = parse_ligand_matrix_from_mol2(mol2_path, sdf_path, lig_pdb_out)
    pdb_indices = ligand_data[4]
    ligand_code = ligand_data[5]

    # Create Networkx Ligand Graph
    ligand_graph = initialize_graph(*ligand_data[:4], ligand_code, pdb_indices=pdb_indices, plot_graph=plot_graphs)

    # Write ligand.txt output
    graph_coverage = {f"{ligand_graph.nodes[x]['atom']}{x}": 0 for x in ligand_graph.nodes()}
    with open(lig_txt_out, "w") as f:
        cg_incidence_counter = defaultdict(int)

        # Select just the combs groups on interest if hydrophobic or hydrophilic flag is set.
        if hydrophobics:
            selected = {x: y for x, y in names_to_smile_map.items() if x in ["ch3", "pro", "isopropyl", "ph", "csc"]}
        elif hydrophilics:
            selected = {x: y for x, y in names_to_smile_map.items() if x not in ["ch3", "pro", "isopropyl", "ph", "csc"]}
        else:
            selected = names_to_smile_map

        # Loop over all (selected) combs_groups.
        for combs_group in selected:
            if verbose:
                print(combs_group.upper())

            # Compute combs_group to ligand matches.
            cg_graph = initialize_graph(*matrix_from_smile(names_to_smile_map[combs_group]), combs_group, plot_graph=plot_graphs)
            all_match_lists = match_graphs(ligand_graph, cg_graph)

            # Write this combs group's matches to text file.
            write_matches_to_ligand_txt_file(
                f, all_match_lists, ligand_code, combs_group, ligand_graph,
                cg_incidence_counter, graph_coverage, verbose
            )

        # Write some summary statistics to the end of the file.
        f.write("\n")
        # f.write("# >Coverage Depth Per Atom: " + graph_coverage.__repr__() + "\n")
        f.write("# >Fraction Heavy Atoms Matched: " + str(len({i:j for i,j in graph_coverage.items() if j != 0 and "H" not in i}) / len({i: j for i,j in graph_coverage.items() if "H" not in i})) + "\n")
        f.write("# >Average Coverage Per Heavy Atom: " + str(sum([j for i,j in graph_coverage.items() if "H" not in i]) / len([j for i,j in graph_coverage.items() if "H" not in i])) + "\n")
        if verbose:
            print("# >Coverage Depth Per Atom: " + graph_coverage.__repr__() + "\n")
            print("# >Fraction Heavy Atoms Matched: " + str(len({i:j for i,j in graph_coverage.items() if j != 0 and "H" not in i}) / len({i: j for i,j in graph_coverage.items() if "H" not in i})) + "\n")
            print("# >Average Coverage Per Heavy Atom: " + str(sum([j for i,j in graph_coverage.items() if "H" not in i]) / len([j for i,j in graph_coverage.items() if "H" not in i])) + "\n")

    # Read in ligand.txt to pandas dataframe and dump it to a .pkl file.
    if write_df_pickle:
        pd_data_frame = read_lig_txt(lig_txt_out)
        pd_data_frame.to_pickle(lig_txt_out[:-4] + ".pkl")


if __name__ == "__main__":
    parser = optparse.OptionParser(description='Generate a ligand.txt file from ligand mol2 file.')
    parser.add_option("--pdb_path", dest='pdb', default="",help="path to pdb file to generate ligand.txt for")
    parser.add_option("--mol2_path", dest="mol2", default="", help="path to mol2 file to convert to pdb and generate ligand.txt for")
    parser.add_option("--sdf_path", dest="sdf", default="", help="(optional) path to sdf file load from if mol2 fails, if unspecified will attempt to load from file with same path as mol2 but ending in .sdf, to prevent redundant loading can input a path to an empty file.")
    parser.add_option("--ligand_pdb_out", dest="lig_pdb_out", default="", help="(optional) path to ligand pdb file that is output while processing the sdf input, defaults to same path as mol2 file but ending in .pdb")
    parser.add_option("--ligand_txt_out", dest="lig_txt_out", default="ligand.txt", help="(optional) path to ligand.txt file output, defaults to ligand.txt in current directory.")
    parser.add_option("-p", "--write_df_pickle", dest="write_df_pickle", action="store_true", default=False, help="(optional) Flag if set writes Pandas dataframe version of ligand.txt to .pkl file of the same name as ligand txt path")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="(optional) Flag if set, prints matches to standard output, default: silenced printing.")
    parser.add_option("-g", "--plot_graphs", dest="plot", action="store_true", default=False, help="(optional) Flag if set, plots networkx graph representation of ligand and fragments. Useful for debugging.")
    parser.add_option("--only_hydrophobics", dest="hydrophobics", action="store_true", default=False, help="(optional) Flag if set, uses only hydrophobic fragments (isopropyl, pro, ch3, csc, ph)")
    parser.add_option("--only_hydrophilics", dest="hydrophilics", action="store_true", default=False, help="(optional) Flag if set, uses only hydrophilic fragments (all fragments other than isopropyl, pro, ch3, csc, ph)")
    (options, args) = parser.parse_args()

    pdb_path = options.pdb_path
    lig_txt_out = options.ligand_txt_out

    mol2_path = options.mol2
    sdf_path = options.sdf
    lig_pdb_out = options.lig_pdb_out
    write_df_pickle = options.write_df_pickle
    verbose = options.verbose
    plot_graphs = options.plot
    hydrophobics = options.hydrophobics
    hydrophilics = options.hydrophilics

    if hydrophobics and hydrophilics:
        print("cannot set both --only_hydrophobics and --only_hydrophilics as these flags are mutually exclusive.")
        exit(1)
    
    main(pdb_path, mol2_path, sdf_path, lig_pdb_out, lig_txt_out, write_df_pickle, verbose, plot_graphs, hydrophobics, hydrophilics)
