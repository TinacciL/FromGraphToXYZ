from mol_lib import MolFromGraphs, mol_graph_image
from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
name = 'C9H5' + '+'
atoms = ['H','H','H','H','H','C','C','C','C','C','C','C','C','C']
bonds = [1,1,1,1,1,1,3,1,3,1,3,1,2]
edges = [[0,5],[1,5],[2,5],[3,13],[4,13],[5,6],[6,7],[7,8],[8,9],[9,10],[10,11],[11,12],[12,13]]
connectivity = []
for i,item in enumerate(edges):
    connectivity.append((edges[i][0],edges[i][1], {'bond': bonds[i]}))
#connectivity = [(edges[0][0],edges[0][1], {'bond': bonds[0]}),(edges[1][0],edges[1][1], {'bond': bonds[1]}),(edges[2][0],edges[2][1], {'bond': bonds[2]}),(edges[3][0],edges[3][1], {'bond': bonds[3]}), (edges[4][0],edges[4][1], {'bond': bonds[4]}), (edges[5][0],edges[5][1], {'bond': bonds[5]})]
path = '/Users/tinaccil/Documents/GitHub/tmp/ions/miss_cations/'
n_dim = len(atoms)
G = nx.Graph()
for i,atom in enumerate(atoms):
    G.add_node(i)
    tmp_attr = {'atom': atom}
    G.nodes[i].update(tmp_attr.copy())
G.add_edges_from(connectivity)
tmp_mol = MolFromGraphs(G)
#mol_graph_image(G)
#rint(Chem.MolToSmiles(tmp_mol))
AllChem.EmbedMolecule(tmp_mol)
AllChem.UFFOptimizeMolecule(tmp_mol)
with open(path + name + '.xyz', "w+") as file_mol:
    file_mol.write(Chem.MolToXYZBlock(tmp_mol))
os.system('ase gui ' + path + name + '.xyz')
