import numpy as np
import os,sys
from coords_from_fchk import info_fchk

class bond_angle:
  def __init__(self, input_folder=None, fchk=None, zmat=False, verbose=False):
    self.path = os.getcwd()
    self.log = input_folder.split('/')[-1] + '.log'
    if(os.path.exists(os.path.join(input_folder, self.log))):
      self.input = os.path.join(self.path, input_folder)
    elif(os.path.exists(os.path.join(input_folder, "lig.log"))):
      self.log = 'lig.log'
      self.input = os.path.join(self.path, input_folder)
    else:
      raise FileError("Files do not exist")
      sys.exit(0)
    self.zmat = zmat
    self.bond = []
    self.bond_value = []
    self.angle = []
    self.angle_value = []
    self.fchk = fchk
    self.verbose = verbose

  def read_bond_angle(self):
    os.chdir(self.input)

    f_log = open(self.log)
    lines = f_log.readlines()
    f_log.close()

    i = 0
    j = 0
    for line in lines:
      if(j != 1 and len(line) > 80 and line[:81] == ' ! Name  Definition              Value          Derivative Info.                !'):
        i = 1
        j = 1
        continue

      if(i == 1):
        terms = line.split()
        if(len(terms) > 6 and terms[1][0] == 'R' and terms[-1] == '!'):
          self.bond.append([int(terms[2].split(',')[0][2:]), int(terms[2].split(',')[1][:-1])])
          bond_AB = np.linalg.norm(self.fchk.coords[self.bond[-1][1]-1]-self.fchk.coords[self.bond[-1][0]-1])
          self.bond_value.append(bond_AB)
          if(self.verbose):
            print("equilibrium bond length")
            print("%5d%5d%11.4f" %(self.bond[-1][0], self.bond[-1][1], bond_AB))
        elif(len(terms) > 6 and terms[1][0] == 'A' and terms[-1] == '!'):
          if(len(terms[2].split(',')) == 3):
            self.angle.append([int(terms[2].split(',')[0][2:]), int(terms[2].split(',')[1]), int(terms[2].split(',')[2][:-1])])
          else:
            self.angle.append([int(terms[2].split(',')[0][2:]), int(terms[2].split(',')[1]), int(terms[2].split(',')[2])])
          bond_AB = self.fchk.coords[self.angle[-1][0]-1] - self.fchk.coords[self.angle[-1][1]-1] # AB length
          bond_CB = self.fchk.coords[self.angle[-1][2]-1] - self.fchk.coords[self.angle[-1][1]-1] # CB length
          angle_ABC = np.degrees(np.arccos(bond_AB.T.dot(bond_CB)/(np.linalg.norm(bond_AB)*np.linalg.norm(bond_CB)))) # ABC angle value (degree)
          self.angle_value.append(angle_ABC)
          if(self.verbose):
            print("equilibrium angle value")
            print("%5d%5d%5d%11.3f" %(self.angle[-1][0], self.angle[-1][1], self.angle[-1][2], angle_ABC))
        elif(len(terms) == 1):
          continue
        else:
          i = 0
    os.chdir(self.path)

