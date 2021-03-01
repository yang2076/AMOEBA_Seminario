usage = '''
	======================================
	assign atom types to an input molecule
	python3 typingAMOEBAplus.py mol.xyz  
	======================================
'''

from pybel import *
import sys

def typing_xyz(xyz):
	if ".xyz" in xyz:
		fileformat = "xyz"
	if ".sdf" in xyz:
		fileformat = "sdf"
	for mol in readfile(fileformat,xyz):
		print(mol)
		matchDict = {}
		matchList = []
		commentsDict = {} 
		natoms = len(mol.atoms)
		for i in range(1, natoms+1, 1):
			matchList.append(i)
		matchDict = dict.fromkeys(matchList, 0)
		commentsDict = dict.fromkeys(matchList, 0)
		lines = open("amoebaplusCHON").readlines()
		for line in lines:
			if "#" not in line[0] and len(line)>10:
				data = line.split()
				myStr = data[0]
				classNum = data[1]
				className = data[2]
				comment = line.split("# ")[1][0:-1]
				smarts = Smarts(myStr)
				match = smarts.findall(mol)
				if match:
					for i in range(len(match)):	
						#matchDict[match[i][0]] = classNum
						matchDict[match[i][0]] = className
						commentsDict[match[i][0]] = comment	
		for i in range(1, natoms+1, 1):
			print("%3s %3s   %s"%(i, matchDict[i], commentsDict[i]))
	return

if len(sys.argv) != 2:
	print(usage)
	sys.exit()
else:
	typing_xyz(sys.argv[1])
