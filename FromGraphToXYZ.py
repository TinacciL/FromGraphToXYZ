from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
from ase import io, neighborlist, atoms

#Function that encode the molecules from graph object to XYZ format
def MolFromGraphs(G):
    # Function based by https://stackoverflow.com/questions/51195392/smiles-from-graph
	'''
	Function that takes as input the networkx graph (each node must have an atom property H,N,C,O etc) and return the mol (rdkit) object
	the function dont discriminate between different type of bond, it care only about the connectivity.
	'''
	adjacency_matrix = nx.convert_matrix.to_numpy_matrix(G,weight='bond').tolist()
	node_list = []
	for i,node in enumerate(G):
		node_list.append(G.nodes[i]['atom'])
	#Create empty editable mol object
	mol = Chem.RWMol()
	#Add atoms to mol and keep track of index
	node_to_idx = {}
	for i in range(len(node_list)):
		a = Chem.Atom(node_list[i])
		molIdx = mol.AddAtom(a)
		node_to_idx[i] = molIdx
	#Add bonds between adjacent atoms
	for ix, row in enumerate(adjacency_matrix):
		for iy, bond in enumerate(row):
			#Only traverse half the matrix
			if iy >= ix:
				break
			#Add relevant bond type (there are many more of these)
			if bond == 0:
				continue
			elif bond == 1:
				bond_type = Chem.rdchem.BondType.SINGLE
				mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
			elif bond == 2:
			    bond_type = Chem.rdchem.BondType.DOUBLE
			    mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
			elif bond == 3:
			    bond_type = Chem.rdchem.BondType.TRIPLE
			    mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
	mol.UpdatePropertyCache(strict=False)
    #Add the possibility in N bearing species cations to have valence 4 for the N center
	for at in mol.GetAtoms():
		if at.GetAtomicNum() == 7 and at.GetExplicitValence()==4 and at.GetFormalCharge()==0:
			at.SetFormalCharge(1)     
	Chem.SanitizeMol(mol)
	return mol

#INPUT - Change the above line in order to save the computed .xyz in the choose path    
path = '/Users/aaa/Documents/FromGraphToXYZ/'
#INPUT - Name of the molecules for the .xyz file
name = 'CH4'
#INPUT - Ordinated list of atoms in the choose molecule
atoms = ['H','H','H','H','C']
#INPUT - Ordinated list of bond in the choose molecule (1 = sigma, 2 = pi, 3 = 2 * pi)
bonds = [1,1,1,1]
#INPUT - Ordinated order of Connettivity between list of atom
edges = [[0,4],[1,4],[2,4],[3,4]]

#Upload all the info in graph object
connectivity = []
for i,item in enumerate(edges):
    connectivity.append((edges[i][0],edges[i][1], {'bond': bonds[i]}))
n_dim = len(atoms)
G = nx.Graph()
for i,atom in enumerate(atoms):
    G.add_node(i)
    tmp_attr = {'atom': atom}
    G.nodes[i].update(tmp_attr.copy())
G.add_edges_from(connectivity)
#Encode the molecules from graph object to XYZ format
tmp_mol = MolFromGraphs(G)
#Read the molecules in RdKit
AllChem.EmbedMolecule(tmp_mol)
#Optimized the guessed structure with UFF in RdKit
AllChem.UFFOptimizeMolecule(tmp_mol)
#Save the .XYZ file of the optimized structure
with open(path + name + '.xyz', "w+") as file_mol:
    file_mol.write(Chem.MolToXYZBlock(tmp_mol))
