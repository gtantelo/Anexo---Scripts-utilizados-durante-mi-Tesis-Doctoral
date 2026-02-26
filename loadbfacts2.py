from pymol import cmd, stored, math
	
def loadBfacts (mol,startaa=1,source="newBfactors.txt", visual="Y"):
	"""
	Replaces B-factors with a list of values contained in a plain txt file
	
	usage: loadBfacts mol, [startaa, [source, [visual]]]
 
	mol = any object selection (within one single object though)
	startaa = number of first amino acid in 'new B-factors' file (default=1)
	source = name of the file containing new B-factor values (default=newBfactors.txt)
	visual = redraws structure as cartoon_putty and displays bar with min/max values (default=Y)
 
	example: loadBfacts 1LVM and chain A
	"""
	obj=cmd.get_object_list(mol)[0]
	counter=int(startaa)

	if visual=="Y":
		cmd.spectrum("b","rainbow", "%s and n. CA " %mol)
		cmd.ramp_new("count", obj, [0, 1], "rainbow")
		cmd.recolor()

cmd.extend("loadBfacts", loadBfacts);