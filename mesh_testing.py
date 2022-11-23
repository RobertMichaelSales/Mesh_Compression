""" Created: 23.11.2022  \\  Updated: 23.11.2022  \\   Author: Robert Sales """

#==============================================================================
# Import libraries

from vtk_to_npy import *

#==============================================================================

vtk_filename = "/home/rms221/Documents/Mesh_Compression/MeshComp_AuxFiles/test_vol.vts"
resolution = (133,49,49)
variable_name = "V"

volume = VTKtoNPY(vtk_filename, variable_name, resolution)

from geomdl import fitting
from geomdl.visualization import VisMPL as vis

points = np.reshape(volume[:, 0, 0,:2],(-1,2)).tolist()
curve = fitting.interpolate_curve(points=points,degree=3)

import numpy as np
import matplotlib.pyplot as plt
evalpts = np.array(curve.evalpts)
pts = np.array(points)
plt.plot(evalpts[:, 0], evalpts[:, 1])
plt.scatter(pts[:, 0], pts[:, 1], color="red")
plt.show()

#==============================================================================