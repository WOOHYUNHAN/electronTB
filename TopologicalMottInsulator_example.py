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

lat = [[1,   -np.sqrt(3)], [1,   np.sqrt(3)]]
orb = [[1.0/3, 2.0/3], [2.0/3, 1.0/3]]
spinor = False
TMI = electronTB(2, 2, spinor)
TMI.set_geometry(lat, orb)
TMI.set_onsite([.0,.0])
NN_hopping = -1.0
TMI.set_hopping(0,1,[0.0,0.0],NN_hopping)
TMI.set_hopping(0,1,[-1.0,0.0],NN_hopping)
TMI.set_hopping(0,1,[0.0,1.0],NN_hopping)
#TMI.print_info()
################################################################################
q_path_line = [[0, 0], [0.5, 0.0], [1.0/3, 1.0/3], [0.0, 0]]
q_spacing_line = 20
#TMI.get_electron_band(q_path_line, q_spacing_line)
#TMI.draw_electron_band()

################################################################################
q_mesh = [30,30]
max_steps = 60
threshold = 1e-4
filling_factor = 1.0/2 #half-filled case
temperature = 1e-9
#print TMI.spin

#print np.dot(TMI.latt_vec.T, TMI.recip_vec)

MFT_TMI = ManybodyInteraction_MFT(TMI)

MFT_TMI.set_order_parameter_type(0,0, [0,0], 'real')
MFT_TMI.set_order_parameter_type(1,1, [0,0], 'real')

MFT_TMI.set_order_parameter_type(1,1, [1,1], 'complex')
#MFT_TMI.set_order_parameter_type(1,1, [0,1], 'complex')
#MFT_TMI.set_order_parameter_type(1,1, [1,0], 'complex')

MFT_TMI.set_order_parameter_type(0,0, [-1,-1], 'complex')
#MFT_TMI.set_order_parameter_type(0,0, [0,-1], 'complex')
#MFT_TMI.set_order_parameter_type(0,0, [-1,0], 'complex')

V_1 = -1 * NN_hopping / 100000000000000000.0
V_2 = -1 * NN_hopping / 0.4

MFT_TMI.set_MB_interactions(1, 1, [0,0], V_1, 0)
MFT_TMI.set_MB_interactions(0, 0, [0,0], V_1, 1)

MFT_TMI.set_MB_interactions(1, 1, [-1,-1], V_2, 2)
MFT_TMI.set_MB_interactions(1, 1, [0,1], V_2, 2)
MFT_TMI.set_MB_interactions(1, 1, [1,0], V_2, 2)

MFT_TMI.set_MB_interactions(0, 0, [1,1], V_2, 3)
MFT_TMI.set_MB_interactions(0, 0, [0,-1], V_2, 3)
MFT_TMI.set_MB_interactions(0, 0, [-1,0], V_2, 3)

for i in range(51):
	MFT_TMI.sc_solver(q_mesh, max_steps, threshold, filling_factor, temperature)


TMI_new = electronTB(2, 2, spinor)
TMI_new.inherit_info(MFT_TMI)
TMI_new.get_electron_band(q_path_line, q_spacing_line)
TMI_new.draw_electron_band()
