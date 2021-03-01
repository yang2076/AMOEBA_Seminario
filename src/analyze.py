import os, sys
import numpy as np
from molecule import molecules
from force_constant import force_constant 

# This class is used to average all the bond length, angles and corresponding force constant 
# with the same types.

class analyze:

  def __init__(self, mol_class = None, b_a_list = [], b_f_con_list = [], a_f_con_list = []):
    self.mol_class = mol_class
    self.b_a_list = b_a_list
    self.b_f_con_list = b_f_con_list
    self.a_f_con_list = a_f_con_list
    self.bond_f_con = {}
    self.angle_f_con = {}
    self.b_f_con_ave = {}
    self.a_f_con_ave = {}
    self.b_f_con_std = {}
    self.a_f_con_std = {}
    self.bond_val = {}
    self.angle_val = {}
    self.b_val_ave = {}
    self.a_val_ave = {}
    self.b_val_std = {}
    self.a_val_std = {}

  # force constant average
  def f_con_average(self):
    i0 = 0
    for i in self.mol_class.mol_name.values():
      b_a = self.b_a_list[i0]
      print(len(b_a.bond))
      # bond list
      for j, k in enumerate(b_a.bond):
        b_1 = k[0]
        b_2 = k[1]
        b_1 = self.mol_class.mol_gen_types[i][b_1-1]
        b_2 = self.mol_class.mol_gen_types[i][b_2-1]
        if((b_1,b_2) in self.bond_f_con.keys()):
          self.bond_f_con[(b_1,b_2)].append(self.b_f_con_list[i0][j])
        elif((b_2,b_1) in self.bond_f_con.keys()):
          self.bond_f_con[(b_2,b_1)].append(self.b_f_con_list[i0][j])
        else:
          self.bond_f_con[(b_1,b_2)] = [self.b_f_con_list[i0][j]]

      # angle list
      for j, k in enumerate(b_a.angle):
        a_1 = k[0]
        a_2 = k[1]
        a_3 = k[2]
        a_1 = self.mol_class.mol_gen_types[i][a_1-1]
        a_2 = self.mol_class.mol_gen_types[i][a_2-1]
        a_3 = self.mol_class.mol_gen_types[i][a_3-1]
        if((a_1,a_2,a_3) in self.angle_f_con.keys()):
          self.angle_f_con[(a_1,a_2,a_3)].append(self.a_f_con_list[i0][j])
        elif((a_3,a_2,a_1) in self.angle_f_con.keys()):
          self.angle_f_con[(a_3,a_2,a_1)].append(self.a_f_con_list[i0][j])
        else:
          self.angle_f_con[(a_1,a_2,a_3)] = [self.a_f_con_list[i0][j]]
      i0 += 1
    # average and std
    count = 0
    for j, k in self.bond_f_con.items():
      k = np.array(k)
      ave = np.average(k)
      std = np.std(k)
      print("bond  %10s%10s%11.4f PRM_%-5d %11.4f" %(j[0], j[1], ave, count, std))
      self.b_f_con_ave[j] = ave
      self.b_f_con_std[j] = std
      count += 1
    for j, k in self.angle_f_con.items():
      k = np.array(k)
      ave = np.average(k)
      std = np.std(k)
      print("angle %10s%10s%10s%11.4f PRM_%-5d %11.4f" %(j[0], j[1], j[2], ave, count, std))
      self.a_f_con_ave[j] = ave
      self.a_f_con_std[j] = std
      count += 1

  # equilibrium bond length and angle average
  def val_average(self):
    i0 = 0
    for i in self.mol_class.mol_name.values():
      b_a = self.b_a_list[i0]
      # bond list
      for j, k in enumerate(b_a.bond):
        b_1 = k[0]
        b_2 = k[1]
        b_1 = self.mol_class.mol_gen_types[i][b_1-1]
        b_2 = self.mol_class.mol_gen_types[i][b_2-1]
        if((b_1,b_2) in self.bond_val.keys()):
          self.bond_val[(b_1,b_2)].append(b_a.bond_value[j])
        elif((b_2,b_1) in self.bond_val.keys()):
          self.bond_val[(b_2,b_1)].append(b_a.bond_value[j])
        else:
          self.bond_val[(b_1,b_2)] = [b_a.bond_value[j]]
      # angle list
      for j, k in enumerate(b_a.angle):
        a_1 = k[0]
        a_2 = k[1]
        a_3 = k[2]
        a_1 = self.mol_class.mol_gen_types[i][a_1-1]
        a_2 = self.mol_class.mol_gen_types[i][a_2-1]
        a_3 = self.mol_class.mol_gen_types[i][a_3-1]
        if((a_1,a_2,a_3) in self.angle_val.keys()):
          self.angle_val[(a_1,a_2,a_3)].append(b_a.angle_value[j])
        elif((a_3,a_2,a_1) in self.angle_val.keys()):
          self.angle_val[(a_3,a_2,a_1)].append(b_a.angle_value[j])
        else:
          self.angle_val[(a_1,a_2,a_3)] = [b_a.angle_value[j]]
      i0 += 1
    # average and std
    for j, k in self.bond_val.items():
      k = np.array(k)
      ave = np.average(k)
      std = np.std(k)
      print("%10s%10s%11.4f%11.4f" %(j[0], j[1], ave, std))
      self.b_val_ave[j] = ave
      self.b_val_std[j] = std
    for j, k in self.angle_val.items():
      k = np.array(k)
      ave = np.average(k)
      std = np.std(k)
      print("%10s%10s%10s%11.4f%11.4f" %(j[0], j[1], j[2], ave, std))
      self.a_val_ave[j] = ave
      self.a_val_std[j] = std

