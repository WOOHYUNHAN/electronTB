import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import *
import numpy as np
import matplotlib.pyplot as plt
import cmath
from matplotlib.collections import LineCollection
import matplotlib.cm as cm
import sys
sys.path.insert(0, '/Users/woohyunhan/githup_reposi/electronTB')
from electronTB_HAN import *
#from numeric import *

'''
This is example for previous PRL paper (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.156401)
'''

lat=[[2.0, 0.0],[-1.0, np.sqrt(3)]]
# define coordinates of orbitals
orb = [[0.5, 0.0], [0.0, 0.5], [0.5,0.5]]
spinor = False
TMI = electronTB(2, 3, spinor)
TMI.set_geometry(lat, orb)
TMI.set_onsite([.0,.0,.0])
delta = 0#-0.01
t11 = -1.0
t21 = t11 / 5.0
fh = 1.0
frac = 0 # 0.05
t2c = 0 #0.01jTMI.set_hopping(0,1,[0.0,0.0],NN_hopping)
TMI.set_hopping(0, 1, [1, 0], t11)
TMI.set_hopping(0, 1, [0, -1], t11)

TMI.set_hopping(0, 2, [0, 0], t11)
TMI.set_hopping(0, 2, [0, -1], t11)

TMI.set_hopping(1, 2, [0, 0], t11)
TMI.set_hopping(1, 2, [-1, 0], t11)

TMI.set_hopping(0, 1, [0, 0], t21)
TMI.set_hopping(0, 1, [1, -1], t21)


TMI.set_hopping(0, 2, [1, 0], t21)
TMI.set_hopping(0, 2, [-1, -1],t21)

TMI.set_hopping(1, 2, [0, 1],t21)
TMI.set_hopping(1, 2, [-1, -1],t21)
#TMI.print_info()
################################################################################
q_path_line = [[0, 0], [0.5, 0.0], [1.0/3, 1.0/3], [0.0, 0]]
q_mesh = [100,100]
q_spacing_line = 20
TMI.get_electron_band(q_path_line, q_spacing_line)
TMI.draw_electron_band()
#eig_array = TMI.get_electron_eigval_mesh(q_mesh)
#print np.max(eig_array[2]), np.min(eig_array[2]), np.max(eig_array[2])- np.min(eig_array[2])
#print TMI.find_filling_specific_band(eig_array, 2, 1.95, 0.1)
#TMI.find_chemical_pot_for_specific_filling(eig_array, 2, 0.01, 0.666, 20000)

################################################################################

