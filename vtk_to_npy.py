""" Created: 22.06.2022  \\  Updated: 22.11.2022  \\   Author: Robert Sales """

#==============================================================================

import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

import matplotlib.pyplot as plt

#==============================================================================
# Define a function to conver 

# Note: The reshaping order is non-standard as the VTK order is not contiguous
# Note: A transpose is required so index slicing corresponds to the right axes

def VTKtoNPY(vtk_filename,variable_name,resolution):

    # Create a vtk file reader object to read structured grids
    reader = vtk.vtkXMLStructuredGridReader()
    
    # Set the input filename for the VTK reader
    reader.SetFileName(vtk_filename)
    
    # Update the reader object to obtain values
    reader.Update()
        
    # Initialise a dictionary to store the scalars 
    variables = {}
    
    # Iterate through each point array
    for i in range(reader.GetNumberOfPointArrays()):
        
        # Get the variable name
        label = reader.GetOutput().GetPointData().GetArrayName(i)
        
        # Get the scalar values
        value = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(i))
        
        # Add the name and values to a dictionary
        variables[label] = value
        
    # Check the number of mesh points match the user-provided resolution
    if (reader.GetNumberOfPoints() != np.prod(resolution)): 
        print("ERROR: NUMBER OF MESH POINTS DOES NOT MATCH USER'S RESOLUTION.")
        return None
    
    if (variable_name not in list(variables.keys())):
        print("ERROR: VARIABLE NAME NOT IN LIST {}.".format(list(variables.keys())))    
        return None
    
    # Get the x,y,z coordinates of each grid point
    volume = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())
    
    # Get the scalar values at each grid point 
    values = np.expand_dims(variables[variable_name],axis=-1)
    
    # Cast the arrays to float32 precision
    volume = volume.astype(np.float32)
    values = values.astype(np.float32)
    
    # Concatenate volume and values arrays
    npy_data = np.concatenate((volume,values),axis=-1)
        
    # Reshape the array according to the user's resolution
    npy_data = np.reshape(npy_data,(resolution[2],resolution[0],resolution[1],-1),order="F")
    
    # Transpose from z,x,y -> x,y,z index ordering
    npy_data = np.transpose(npy_data,(1,2,0,3))
    
    return npy_data

#==============================================================================
# Define a function to find the external boundary of a structured mesh/grid

def FindBoundaries(volume):
    
    x_end1,x_end2 = volume[0,:,:,:],volume[-1,:,:,:]
    y_end1,y_end2 = volume[:,0,:,:],volume[:,-1,:,:]
    z_end1,z_end2 = volume[:,:,0,:],volume[:,:,-1,:]
    
    # Create a figure with 3-dimensional axes    
    fig = plt.figure(figsize=(28, 28))
    ax = fig.add_subplot(projection='3d')
    
    # Plot the x-end boundaries    
    ax.scatter(xs=x_end1[...,0], ys=x_end1[...,1], zs=x_end1[...,2], c=x_end1[...,3])
    ax.scatter(xs=x_end2[...,0], ys=x_end2[...,1], zs=x_end2[...,2], c=x_end2[...,3])
    
    # Plot the y-end boundaries    
    ax.scatter(xs=y_end1[...,0], ys=y_end1[...,1], zs=y_end1[...,2], c=y_end1[...,3])
    ax.scatter(xs=y_end2[...,0], ys=y_end2[...,1], zs=y_end2[...,2], c=y_end2[...,3])
    
    # Plot the y-end boundaries    
    ax.scatter(xs=z_end1[...,0], ys=z_end1[...,1], zs=z_end1[...,2], c=z_end1[...,3])
    ax.scatter(xs=z_end2[...,0], ys=z_end2[...,1], zs=z_end2[...,2], c=z_end2[...,3])
    
    # Label each of the axes
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    # Display the figure 
    plt.show()
    return None

#==============================================================================
# Define a function to stack the external boundary of a structured mesh/grid

def StackBoundaries(volume):
    
    x_end1 = np.reshape(volume[ 0,:,:,:],(-1,volume.shape[-1]))
    x_end2 = np.reshape(volume[-1,:,:,:],(-1,volume.shape[-1]))
    y_end1 = np.reshape(volume[:, 0,:,:],(-1,volume.shape[-1]))
    y_end2 = np.reshape(volume[:,-1,:,:],(-1,volume.shape[-1]))
    z_end1 = np.reshape(volume[:,:, 0,:],(-1,volume.shape[-1]))
    z_end2 = np.reshape(volume[:,:,-1,:],(-1,volume.shape[-1]))
    
   ####################################################################################################### YOU GOT HERE
   ####################################################################################################### YOU GOT HERE
   ####################################################################################################### YOU GOT HERE
   ####################################################################################################### YOU GOT HERE
   ####################################################################################################### YOU GOT HERE
    
#==============================================================================
# Define a function to plot the external boundary of a structured mesh/grid

def PlotBoundaries(volume):
    
    # %matplotlib inline
    # %matplotlib auto    

    # Create a figure with 3-dimensional axes    
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(projection='3d')
    
    # Plot the x-end boundaries
    ax.plot_surface(X=volume[ 0,:,:,0], Y=volume[ 0,:,:,1], Z=volume[ 0,:,:,2])
    ax.plot_surface(X=volume[-1,:,:,0], Y=volume[-1,:,:,1], Z=volume[-1,:,:,2])
    ax.scatter(xs=volume[ 0,:,:,0], ys=volume[ 0,:,:,1], zs=volume[ 0,:,:,2], c=volume[ 0,:,:,3])
    ax.scatter(xs=volume[-1,:,:,0], ys=volume[-1,:,:,1], zs=volume[-1,:,:,2], c=volume[-1,:,:,3])    

    # Plot the y-end boundaries
    ax.plot_surface(X=volume[:, 0,:,0], Y=volume[:, 0,:,1], Z=volume[:, 0,:,2])
    ax.plot_surface(X=volume[:,-1,:,0], Y=volume[:,-1,:,1], Z=volume[:,-1,:,2])
    ax.scatter(xs=volume[:, 0,:,0], ys=volume[:, 0,:,1], zs=volume[:, 0,:,2], c=volume[:, 0,:,3])
    ax.scatter(xs=volume[:,-1,:,0], ys=volume[:,-1,:,1], zs=volume[:,-1,:,2], c=volume[:,-1,:,3])

    # Plot the z-end boundaries
    ax.plot_surface(X=volume[:,:, 0,0], Y=volume[:,:, 0,1], Z=volume[:,:, 0,2])
    ax.plot_surface(X=volume[:,:,-1,0], Y=volume[:,:,-1,1], Z=volume[:,:,-1,2])
    ax.scatter(xs=volume[:,:, 0,0], ys=volume[:,:, 0,1], zs=volume[:,:, 0,2], c=volume[:,:, 0,3])
    ax.scatter(xs=volume[:,:,-1,0], ys=volume[:,:,-1,1], zs=volume[:,:,-1,2], c=volume[:,:,-1,3])

    # Label each of the axes
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    # Display the figure 
    plt.show()
    return None

#==============================================================================


vtk_filename = "/home/rms221/Documents/Mesh_Compression/MeshComp_AuxFiles/test_vol.vts"
resolution = (133,49,49)
variable_name = "V"

FindBoundaries(VTKtoNPY(vtk_filename, variable_name, resolution))

# PlotBoundary(VTKtoNPY(vtk_filename, variable_name, resolution))
