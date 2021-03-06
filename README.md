# FromGraphToXYZ
This is a Python 3.8 script that permits to start from a molecule graph and computed a XYZ configuration.

## Installation

This program used a python3 interface, to run this code you must install on your machine this list of packages:

* ```networkx```
* ```rdkit```
* ```ASE```

## Script Structure

0 - Import the required packages

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
from ase import io, neighborlist, atoms
```

1 - Function that encode the molecules from graph object to .XYZ format, and you will call it in the main of the code.

```python
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
```
2 - Code main: computing the .XYZ optimized structure at UFF level from molecular graph object.

   2.1 - Input the parameters in order to create an Molecular graph.
```python
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
```

   2.2 - Computed the info into molecular graph object

```python
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
```
   2.3 - Encode the molecules from graph object to XYZ format
```python
tmp_mol = MolFromGraphs(G)
```
   2.4 - Read the molecules in RdKit
```python
AllChem.EmbedMolecule(tmp_mol)
```
   2.5 - Optimized the guessed structure with UFF in RdKit
```python
AllChem.UFFOptimizeMolecule(tmp_mol)
```
   2.6 - Save the .XYZ file of the optimized structure
```python
with open(path + name + '.xyz', "w+") as file_mol:
    file_mol.write(Chem.MolToXYZBlock(tmp_mol))
```
## Associated publication
Please reffer to the following publication to cite our work or retrieve the info:

[Structures and Properties of Known and Postulated Interstellar Cations, L.Tinacci et al 2021 ApJS 256 35](https://doi.org/10.3847/1538-4365/ac194c) 

(https://doi.org/10.3847/1538-4365/ac194c)

## Acknowledgments
This project has received funding within the European Union???s Horizon 2020 research and innovation programme from the Marie Sklodowska-Curie for the project ???Astro-Chemical Origins??? (ACO), grant agreement No 811312.
