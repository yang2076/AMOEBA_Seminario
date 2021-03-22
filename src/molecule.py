import os, sys
import numpy as np
from chemFileConvert import *
from pybel import *

class molecules:

  def __init__(self, mol_dir = '../dir', mol_list = 'filelist_all.txt', given_class=True, given_txyz=False, H_comb=True, C_sp3_H_comb=True, style='number'):
    self.path = os.getcwd()
    self.mol_dir = os.path.join(self.path, mol_dir)
    self.mol_list = os.path.join(self.mol_dir, mol_list)
    self.given_class = given_class
    self.gen_types = {}
    self.mol_gen_types = {}

    self.type_conclusion = 'typing.txt' # conclude all the atomic types and correponding relations of the whole databse in one file 

    self.given_txyz = given_txyz
    self.keywords = []

    self.style = style # decide what kind of style you want to show the parameter
    # style = 'number' e.g. bond          13        14
    # style = 'word'   e.g. bond          HW        OW

    # Count and record molecules
    f = open(self.mol_list)
    lines = f.readlines()
    f.close()
    self.mol_num = 0
    self.mol_name = {}
    for line in lines:
      terms = line.split()
      if(len(terms) != 0):
        self.mol_num += 1
        self.mol_name[self.mol_num] = terms[0]
        words = " Guess=INDO MaxDisk=100GB freq"
        if(len(terms) == 2 and terms[1] == 'i'):
          self.keywords.append(words + " SCRF=(PCM)")
        else:
          self.keywords.append(words)

    # Convert to xyz
    count = 0
    for i in self.mol_name.values():
      mol_path = os.path.join(self.mol_dir, i)
      os.chdir(mol_path)
      LOG2XYZ(i + '.log')
      LOG2COM(i + '.log', extraKeyword = self.keywords[count], write_new_file=False)

    count = 0
    for i in self.mol_name.values():
      mol_path = os.path.join(self.mol_dir, i)
      os.chdir(mol_path)
      LOG2XYZ(i + '.log')
      LOG2COM(i + '.log', extraKeyword = self.keywords[count], write_new_file=False)
      conn = open("connect",'w')
      txyz = open(i + '.txyz')
      lines_1 = txyz.readlines()
      for line1 in lines_1:
        terms = line1.split()
        num = len(terms)
        if(num == 7):
          conn.write("%9s %4s\n" %(terms[5],terms[6]))
        elif(num == 8):
          conn.write("%9s %4s %4s\n" %(terms[5],terms[6],terms[7]))
        elif(num == 9):
          conn.write("%9s %4s %4s %4s\n" %(terms[5],terms[6],terms[7],terms[8]))
        elif(num == 10):
          conn.write("%9s %4s %4s %4s %4s\n" %(terms[5],terms[6],terms[7],terms[8],terms[9]))
        else:
          conn.write('\n')
      conn.close()
      COM2TXYZ(i + '.com','connect')
      count += 1
    os.chdir(self.path)

    # Combining all the H as the same, so do the ordinary sp3 C (C33, C32, C31, C30)
    comb = self.types_equi(H_comb, C_sp3_H_comb)
    if(comb[-1] == 3):
      self.H_comb = comb[0]
      self.C_sp3_H_comb = comb[1]
    elif(comb[-1] == 2):
      self.H_comb = []
      self.C_sp3_H_comb = comb[0]
    elif(comb[-1] == 1):
      self.H_comb = comb[0]
      self.C_sp3_H_comb = []
    else:
      self.H_comb = []
      self.C_sp3_H_comb = []

  def typing_xyz(self):
    # Read the general type
    lines = open("amoebaplus_type").readlines()
    types = {}
    i = 0
    for line in lines:
      if("#" not in line[0] and len(line) > 10):
        i += 1
        data = line.split()
        myStr = data[0]

        global_index = data[1]
        classNum = int(data[2])
        className = data[3]

        # Decide if the atom class is automatically generated in order without any combination
        if(not self.given_class):
          if(className not in types):
            classNum = i
            types[className] = classNum
          else:
            classNum = types[className]

        comment = line.split("# ")[1][0:-1]
        smarts = Smarts(myStr)
        self.gen_types[global_index] = (classNum, className, comment, smarts)


    # Convert to type
    conclusion = open(self.type_conclusion, 'w')
    count = 1
    for i in self.mol_name.values():
      mol_path = os.path.join(self.mol_dir, i)
      os.chdir(mol_path)
      f_type = open("%s.type" %(i), 'w')
      # Read xyz
      if(not self.given_txyz):
        for mol in readfile("xyz", i + '.xyz'):
          matchDict = {}
          self.mol_gen_types[i] = []
          natoms = len(mol.atoms)

          for j1, j2 in self.gen_types.items():
            smarts_0 = j2[3]
            match = smarts_0.findall(mol)
            if(match):
              for k in range(len(match)):
                matchDict[match[k][0]] = j2

          for j in range(1, natoms+1, 1):
            f_type.write("%5d%11s%5d%40s\n" %(j, matchDict[j][1], matchDict[j][0], matchDict[j][2]))

            # specific index is given to H and ordinary sp3 C
            if(matchDict[j][1] in self.H_comb):
              index = 0
              H_type = 'H'
              self.mol_gen_types[i].append(H_type)
            elif(matchDict[j][1] in self.C_sp3_H_comb):
              index = -1
              C_type = 'C3*'
              self.mol_gen_types[i].append(C_type)
            else:
              self.mol_gen_types[i].append(matchDict[j][1])

            conclusion.write("%5s%5d%5d%11s%5d\n" %(i[:3], count, j, matchDict[j][1], matchDict[j][0]))
            count += 1
      # Read txyz
      else:
        for mol in readfile("xyz", i + '.xyz'):
          print(i)
          [atoms, coord, order, types, connections] = readTXYZ(i + '.txyz')
          matchDict = {}
          self.mol_gen_types[i] = []
          natoms = len(mol.atoms)

          for j1, j2 in self.gen_types.items():
            smarts_0 = j2[3]
            match = smarts_0.findall(mol)
            if(match):
              for k in range(len(match)):
                matchDict[match[k][0]] = j2

          for j in range(1, natoms+1, 1):
            f_type.write("%5d%11s%5d%40s\n" %(j, matchDict[j][1], matchDict[j][0], matchDict[j][2]))

            if(matchDict[j][1] in self.H_comb):
              index = 0
              H_type = 'H'
              self.mol_gen_types[i].append(H_type)
            elif(matchDict[j][1] in self.C_sp3_H_comb):
              index = 1000
              C_type = 'C3*'
              self.mol_gen_types[i].append(C_type)
            else:
              self.mol_gen_types[i].append(matchDict[j][1])

            conclusion.write("%5s%5s%5d%11s%5d\n" %(i[:3], types[j-1], j, matchDict[j][1], matchDict[j][0]))

        f_type.close()
    conclusion.close()
    os.chdir(self.path)

  # This function is used to create equivalent sets for general types
  # E.g. HC, HO, HS, HN can be combined as H; C33, C32, C31, C30 can be combined as C3*.
  def types_equi(self, H_comb, C_sp3_comb):
    # Here list the original classifications.
    # All the H are combined
    H_set = ['HC', 'HO', 'HS', 'HN', 'HP','HCH3','HC2','HCH2','HCH','H+']
    # Different hybridizations for each element have combined without the types in those hydrogenated molecules
    # E.g. C in CH4 is named as C34, but it is not inside C_sp3 set. Likewise for NH3, H2S and so on.
    C_sp3_set = []
    C_sp3_H_set = ['C33', 'C32', 'C31', 'C30'] # The subset of C_sp3_set
    C_sp2_set = [] # aliphatic
    C_sp2_aro_set = [] # aromatic
    C_sp_set = []
    N_sp3_set = []
    N_sp2_set = []
    N_sp2_aro_set = []
    N_sp_set = []
    O_sp3_set = []
    O_sp2_set = []
    O_sp2_aro_set = []
    S_sp3_set = []
    S_sp2_set = []
    P_set = []
    X_set = [] # halogen

    combine = []
    num = 0
    if(H_comb):
      combine.append(H_set)
      num += 1
    if(C_sp3_comb):
      combine.append(C_sp3_H_set)
      num += 2

    combine.append(num)
    return combine
