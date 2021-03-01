import os, sys
from molecule import molecules
from coords_from_fchk import info_fchk
from bond_angle_lists import bond_angle
from force_constant import force_constant 
from analyze import analyze

mol_list = sys.argv[1]
mol = molecules(mol_list, given_txyz=True,H_comb=False,C_sp3_H_comb=False)
mol.typing_xyz()
b_a_list = []
b_f_con_list = []
a_f_con_list = []
for name in mol.mol_name.values():
  fchk_1 = info_fchk('../dir/' + name, name + '.fchk')
  fchk_1.info_anal()
  bond_angle_1 = bond_angle('../dir/' + name, fchk_1)
  bond_angle_1.read_bond_angle()
  b_a_list.append(bond_angle_1)
  force_constant_1 = force_constant('../dir/' + name, fchk_1, bond_angle_1, scale=0.943)
  force_constant_1.eigen()
  force_constant_1.bond_projection()
  force_constant_1.angle_projection()
  b_f_con_list.append(force_constant_1.k_b)
  a_f_con_list.append(force_constant_1.k_a)

ana = analyze(mol, b_a_list, b_f_con_list, a_f_con_list)
ana.f_con_average()
ana.val_average()
