import os,sys
import subprocess
import numpy as np

class info_fchk:
  # Extract the necessary info. from .fchk file
  # (1). atomic num (or index)
  # (2). atomic coordinates
  # (3). Hessian matrix
  def __init__(self, input_folder=None, fchk_file=None):
    self.path = os.getcwd()
    self.input = os.path.join(self.path, input_folder)
    self.fchk = os.path.join(self.input, fchk_file)
    if(not os.path.exists(os.path.join(input_folder, fchk_file))):
      chk = os.path.join(self.input, fchk_file.split('.')[0] + '.chk')
      os.system("formchk %s %s" %(chk, self.fchk))
    self.atm_num = []
    self.coords = []
    self.hessian = []

  def info_anal(self):
    os.chdir(self.input)

    # Get the atom number
    atm = os.popen("grep 'Atomic numbers' %s" %(self.fchk)).read()
    self.tot_num = int(atm.split()[4])
    list_atm = os.popen("grep -A%d 'Atomic numbers' %s" %(np.ceil(self.tot_num/6.), self.fchk)).read()
    self.atm_num = [int(i) for i in list_atm.split()[5:]]
    
    # Get the coordinates
    crd = os.popen("grep 'Current cartesian coordinates' %s" %(self.fchk)).read()
    self.crd_num = int(crd.split()[5])
    list_crd = os.popen("grep -A%d 'Current cartesian coordinates' %s" %(np.ceil(self.crd_num/5.), self.fchk)).read()
    self.coords = [float(i) for i in list_crd.split()[6:]]
    self.coords = np.array(self.coords) * 0.529  # Bohr to Ang
    self.coords = self.coords.reshape(self.tot_num, 3)

    # Get the hessian
    hes = os.popen("grep 'Cartesian Force Constants' %s" %(self.fchk)).read()
    print(self.input)
    self.hes_num = int(hes.split()[5])
    list_hes = os.popen("grep -A%d 'Cartesian Force Constants' %s" %(np.ceil(self.hes_num/5.), self.fchk)).read()
    unprocessed_hessian = [float(i) for i in list_hes.split()[6:]]
    unprocessed_hessian = np.array(unprocessed_hessian) * 627.509391 / (np.square(0.529)) # Hartree/bohr to kcal/mol/ang
    self.hessian = np.zeros((3*self.tot_num, 3*self.tot_num))
    for i in range(3*self.tot_num):
      for j in range(i+1):
        self.hessian[i][j] = unprocessed_hessian[int(i*(i+1)/2+j)]
        self.hessian[j][i] = unprocessed_hessian[int(i*(i+1)/2+j)]

    os.chdir(self.path)

  def general_typing(self):
    os.chdir(self.input)
    self.types = []

    f_types = open("%s.type" %(self.input.split('/')[-1]))
    for line in f_types:
      terms = line.split()
      self.types.append(terms[2])
    f_types.close()

    os.chdir(self.path)
