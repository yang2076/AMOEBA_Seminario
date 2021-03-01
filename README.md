# Modified Seminario on the style of AMOEBA(+)
A program for computing initial stretching and bending parameters (good guess) based on the style of AMOEBA+

Please cite the following reference:
1. Allen, A. E. A.; Payne, M. C.; Cole, D. J., Harmonic Force Constants for Molecular Mechanics Force Fields via Hessian Matrix Projection. Journal of Chemical Theory and Computation 2018, 14 (1), 274-281.
2.

# The purpose of this program
It will automatically derive the bond and angle parameters from QM Hessian, which shows sufficient efficiency and accuracy
to give us reliable initial guess for valence parametrization.
The modified Seminario is dependent on the classical simple harmonic model, which indicates that we can not directly use
the parameters because AMOEBA+ is using MM3-style valence terms (up to quartic for stretching and sextic for bending)

# The workflow of the whole program

(1). Some preparation needs to be done prior to using this program:
[1]. Use Gaussian to calculate all the molecules in your database in folder ./dir with each molecule using one single
subfolder. Please turn the opt=Tight and freq on. After calculation, please use formchk to generate .fchk file

*** Make sure you have set the correct scaling factor for frequency calculated by different QM method and bsis set
For MP2/6-31g*, the scaling factor is 0.943. (You can find the appropriate scaling factor in: https://cccbdb.nist.gov/vibscalejust.asp)
change scaling factor in dowork.py if necessary

[2]. Make a list of all the molecules in ./dir/filelist_all.txt
[3]. You need to produce correct tinker .txyz file for each molecule in subfolder with atomic index in order you made in file list
(e.g. in 1st molecule txyz, you already set indices 1~12, then in 2nd txyz, the indices should begin from 13)
(The type file will be automatically printed out by operating dowork.py)

(2). The main workflow:
dowork.py -- molecule.py : read all the molecules in filelist_all.txt and the types in amoebaplus_type 
          -- info_fchk.py : read necessary info. (atomic num, coordinates and hessian) from .fchk
          -- bond_angle_lists.py : read all the atomic indices forming the bonds and angles for each molecule, calculate the optimized bond length and angle degree
          -- force_constant.py : calculate the force constants
          -- analyze.py : collect all the results and do average

(3). Output:
In subfolder of dir:
.type : the results of type match
.xyz : the classical xyz file

In src:
results.log:
the results of the force constants and the average reference length and angle.
format:

[1]. Force constant for the bond:
bond   atom_1  atom_2  k_b  RMSE  (unit: kcal/(mol\*A))
[2]. Force constant for the angle:
angle   atom_1  atom_2  atom_3  k_a  RMSE  (unit: kcal/(mol\*rad))
[3]. Force constant for the bond:
bond   atom_1  atom_2  length  RMSE  (unit: A)
[4]. Force constant for the angle:
angle   atom_1  atom_2  atom_3  angle  RMSE  (unit: deg)
type.txt:
the corresponding relation between atomic class
format:
index_of_molecule   atom_index  orderly_local_index   type_shortname   atom_class
(atom_index : determined by tinker after combining those structrally equivalent atoms)
(orderly_local_index : the atomic index in only one single molecule)
e.g.
  000    1    1         OW   13
  000    2    2         HW   14
  000    2    3         HW   14
  001    3    1        C34   36
  001    4    2         HC   20
  001    4    3         HC   20
  001    4    4         HC   20
  001    4    5         HC   20
  002    5    1       N3H3  139
  002    6    2         HN   22
  002    6    3         HN   22
  002    6    4         HN   22
  003    7    1         N3  137
  003    8    2         HN   22
  003    8    3         HN   22
  003    9    4        C33   35
  003   10    5         HC   20
  003   10    6         HC   20
  003   10    7         HC   20
  004   11    1         OH  170
  004   12    2         HO   19
  004   13    3        C33   35
  004   14    4         HC   20
  004   14    5         HC   20
  004   14    6         HC   20
  005   15    1        C33   35
  005   16    2         HC   20
  005   16    3         HC   20
  005   16    4         HC   20
  005   17    5        CON   88
  005   18    6        O=C  178
  005   19    7       NC=O  141
  005   20    8         HN   22
  005   21    9        C33   35
  005   22   10         HC   20
  005   22   11         HC   20
  005   22   12         HC   20
