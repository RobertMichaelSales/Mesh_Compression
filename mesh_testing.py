""" Created: 23.11.2022  \\  Updated: 23.11.2022  \\   Author: Robert Sales """

#==============================================================================
# Import libraries

import matplotlib.pyplot as plt
from vtk_to_npy import *

#==============================================================================

# vtk_filename = "/home/rms221/Documents/Mesh_Compression/MeshComp_AuxFiles/test_vol.vts"
# resolution = (133,49,49)
# variable_name = "V"

# %matplotlib auto

# volume = VTKtoNPY(vtk_filename, variable_name, resolution)

# PlotAll(volume)

# from geomdl import fitting
# from geomdl.visualization import VisMPL as vis

# points = np.reshape(volume[:, 0, 0,:2],(-1,2)).tolist()
# curve = fitting.interpolate_curve(points=points,degree=3)

# import numpy as np
# import matplotlib.pyplot as plt
# evalpts = np.array(curve.evalpts)
# pts = np.array(points)
# plt.plot(evalpts[:, 0], evalpts[:, 1])
# plt.scatter(pts[:, 0], pts[:, 1], color="red")
# plt.show()

#==============================================================================

def GenerateTrianglePoints(n_layers):
    
    import math
    import numpy as np
    
    # Set up an empty list for adding layers/points
    layers = []
    points = []
    
    # Generate points iteratively
       
    for i in range(n_layers):
        
        newpoints = []
        
        if (i==0):
            
            layers.append([(0,0)])

        else:
            
            for x_coord,y_coord in layers[-1]:
                
                new_x_coord = x_coord + 1.0
                new_y_coord = y_coord + 0.0
                
                newpoints.append((new_x_coord,new_y_coord))
            
            new_x_coord = x_coord + (1.0 * math.sin(math.pi/6))
            new_y_coord = y_coord + (1.0 * math.cos(math.pi/6))
            
            newpoints.append((new_x_coord,new_y_coord))
            layers.append(newpoints)
            
    # Turn the list of list of tuples to a list of tuples
    for layer in layers:
        for point in layer:
            points.append(point)

    return np.array(points)

#==============================================================================

from scipy.spatial import Delaunay
import alphashape

points = GenerateTrianglePoints(4)

points = points + (0.5*np.random.rand(*points.shape))

tri = Delaunay(points)
alp = alphashape.alphashape(points)

# Create a figure with 3-dimensional axes    
fig, ax = plt.subplots(1,1,constrained_layout=True,figsize=(12,12))

# Plot the edge points
ax.scatter(points[:,0],points[:,1],c="blue",s=100) 
ax.triplot(points[:,0],points[:,1],tri.simplices)

# Label each of the axes
ax.set_xlabel('x')
ax.set_ylabel('y')

# Display the figure 
plt.show()

#==============================================================================
