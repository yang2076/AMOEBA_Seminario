import numpy as np
import os,sys
from coords_from_fchk import info_fchk
from bond_angle_lists import bond_angle

class force_constant:
  # The main factory to produce the final parameters and average length and angle degree.
  # The details of the calculating process could be found in the ref.
  def __init__(self, input_folder=None, fchk=None, b_a=None, scale=0.957):
    self.path = os.getcwd()
    self.input = os.path.join(self.path, input_folder)
    self.fchk = fchk
    self.b_a = b_a
    self.vec = []
    self.unit_vec = []
    self.scale = scale

  def eigen(self):
    N = self.fchk.tot_num
    self.eigenvectors = np.zeros((N,N,3,3))
    self.eigenvalues = np.zeros((N,N,3))
    for i in range(N):
      for j in range(N):
        v1, v2 = np.linalg.eig(self.fchk.hessian[(i*3):((i+1)*3), (j*3):((j+1)*3)])
        self.eigenvalues[i, j, :] = v1
        self.eigenvectors[i, j, :, :] = v2

  def bond_projection(self):
    # vectors for bond
    self.bond_num = len(self.b_a.bond)
    self.k_b = np.zeros(self.bond_num) # array for all force constants
    for i in range(self.bond_num):
      bond_AB = self.b_a.bond[i]

      # coords index from 0
      vec_AB = self.fchk.coords[bond_AB[1]-1] - self.fchk.coords[bond_AB[0]-1]
      eigenvalues_AB = self.eigenvalues[bond_AB[0]-1, bond_AB[1]-1, :]
      eigenvectors_AB = self.eigenvectors[bond_AB[0]-1, bond_AB[1]-1, 0:3, 0:3]
      eigenvalues_BA = self.eigenvalues[bond_AB[1]-1, bond_AB[0]-1, :]
      eigenvectors_BA = self.eigenvectors[bond_AB[1]-1, bond_AB[0]-1, 0:3, 0:3]
      unit_vec_AB = vec_AB / self.b_a.bond_value[i]
      unit_vec_BA = - vec_AB / self.b_a.bond_value[i]

      for j in range(3):
        self.k_b[i] = self.k_b[i] - 0.5 * eigenvalues_AB[j] * np.abs(unit_vec_AB.dot(eigenvectors_AB[:,j]))
        self.k_b[i] = self.k_b[i] - 0.5 * eigenvalues_BA[j] * np.abs(unit_vec_BA.dot(eigenvectors_BA[:,j]))
      self.k_b[i] *= 0.5 * (self.scale**2)

      #print(bond_AB[0], bond_AB[1], self.k_b[i])
      self.vec.append(vec_AB)
      self.unit_vec.append(unit_vec_AB)

    self.vec = np.array(self.vec)
    self.unit_vec = np.array(self.unit_vec)

  def angle_projection(self):
    # vector for angle
    self.angle_num = len(self.b_a.angle)
    self.k_a = np.zeros((self.angle_num,))
    self.p_dict = {}
    self.p_list_A = []
    self.p_list_C = []
    self.bond_AB = []
    self.bond_CB = []

    for i in range(self.angle_num):
      angle_ABC = self.b_a.angle[i]
      tmp_1 = "%d-%d" %(angle_ABC[1], angle_ABC[0]) # atom B - atom A
      tmp_2 = "%d-%d" %(angle_ABC[1], angle_ABC[2]) # atom B - atom C
      self.p_dict[tmp_1] = []
      self.p_dict[tmp_2] = []

    for i in range(self.angle_num):
      angle_ABC = self.b_a.angle[i]
      vec_AB = self.fchk.coords[angle_ABC[0]-1] - self.fchk.coords[angle_ABC[1]-1]
      vec_CB = self.fchk.coords[angle_ABC[2]-1] - self.fchk.coords[angle_ABC[1]-1]
      unit_AB = vec_AB / np.linalg.norm(vec_AB)
      unit_CB = vec_CB / np.linalg.norm(vec_CB)

      unit_N = np.cross(unit_CB, unit_AB)
      unit_N = unit_N / np.linalg.norm(unit_N)

      unit_PA = np.cross(unit_AB, unit_N)
      unit_PC = np.cross(unit_CB, unit_N)
      unit_PA = unit_PA / np.linalg.norm(unit_PA)
      unit_PC = unit_PC / np.linalg.norm(unit_PC)

      self.p_list_A.append(unit_PA)
      self.p_list_C.append(unit_PC)

      tmp_1 = "%d-%d" %(angle_ABC[1], angle_ABC[0]) # atom B - atom A
      tmp_2 = "%d-%d" %(angle_ABC[1], angle_ABC[2]) # atom B - atom C
      self.p_dict[tmp_1].append(unit_PA)
      self.p_dict[tmp_2].append(unit_PC)
      self.bond_AB.append(np.linalg.norm(vec_AB))
      self.bond_CB.append(np.linalg.norm(vec_CB))

    # modified force constant for angles
    for i in range(self.angle_num):
      angle_ABC = self.b_a.angle[i]
      sum_1 = np.zeros((3,))
      sum_2 = np.zeros((3,))
      eigenvalues_AB = self.eigenvalues[angle_ABC[0]-1, angle_ABC[1]-1, :]
      eigenvalues_CB = self.eigenvalues[angle_ABC[2]-1, angle_ABC[1]-1, :]
      eigenvectors_AB = self.eigenvectors[angle_ABC[0]-1, angle_ABC[1]-1, 0:3, 0:3]
      eigenvectors_CB = self.eigenvectors[angle_ABC[2]-1, angle_ABC[1]-1, 0:3, 0:3]
      unit_PA = self.p_list_A[i]
      unit_PC = self.p_list_C[i]
      for j in range(3):
        eig_AB = eigenvectors_AB[:,j]
        sum_1[j] = eigenvalues_AB[j] * abs(unit_PA.dot(eig_AB.T))
        eig_CB = eigenvectors_CB[:,j]
        sum_2[j] = eigenvalues_CB[j] * abs(unit_PC.dot(eig_CB.T))
      sum_1 = np.sum(sum_1)
      sum_2 = np.sum(sum_2)

      scale_A = 0.
      scale_C = 0.
      # scaling factor
      tmp_1 = "%d-%d" %(angle_ABC[1], angle_ABC[0]) # atom B - atom A
      if(len(self.p_dict[tmp_1]) == 1):
        scale_A = 1.
      else:
        for j in self.p_dict[tmp_1]:
          scale_A += np.dot(unit_PA, j.T) ** 2
        scale_A =  1 + ((scale_A - 1) / (len(self.p_dict[tmp_1]) - 1))
      tmp_2 = "%d-%d" %(angle_ABC[1], angle_ABC[2]) # atom B - atom C
      if(len(self.p_dict[tmp_2]) == 1):
        scale_C = 1.
      else:
        for j in self.p_dict[tmp_2]:
          scale_C += np.dot(unit_PC, j.T) ** 2
        scale_C =  1 + ((scale_C - 1) / (len(self.p_dict[tmp_2]) - 1))

      sum_1 = sum_1/scale_A
      sum_2 = sum_2/scale_C
      k_theta = 1. / ((self.bond_AB[i] ** 2) * sum_1) + 1. / ((self.bond_CB[i] ** 2) * sum_2)
      k_theta = 1./k_theta
      k_theta *= -0.5 * (self.scale**2)
      self.k_a[i] = k_theta

