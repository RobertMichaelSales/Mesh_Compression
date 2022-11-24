""" Created: 22.06.2022  \\  Updated: 22.11.2022  \\   Author: Robert Sales """

#==============================================================================
# Import libraries

import vtk
import numpy as np
import matplotlib.pyplot as plt

from vtk.util.numpy_support import vtk_to_numpy

#==============================================================================

# For the purposes of this code a unit cube is defined with vertcies (x,y,z):
    
# A (0,0,0) 
# B (0,1,0)      
# C (1,0,0)      
# D (1,1,0)       
# E (0,0,1)
# F (0,1,1)
# G (1,0,1)
# H (1,1,1)

#==============================================================================
# Define a function to import data from a VTK file and output a NPY data array 

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
# Define a function to find all the coordinates of a structured 3D mesh

def ExtractAllCoords(volume,plot=True):
    
    all_coords = np.reshape(volume[:,:,:,:3],(-1,3))
    
    if plot:
        # Create a figure with 3-dimensional axes    
        fig = plt.figure(figsize=(28, 28))
        ax = fig.add_subplot(projection='3d')

        # Plot the face points        
        ax.scatter(xs=all_coords[...,0],ys=all_coords[...,1],zs=all_coords[...,2],c="blue",s=0.5)
            
        # Label each of the axes
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
    
        # Display the figure 
        plt.show()
    
    return all_coords

#==============================================================================
# Define a function to plot all the coordinates of a structured 3D mesh

def PlotAll(volume):
    
    all_coords = volume[:,:,:,:3]
    
    # Create a figure with 3-dimensional axes    
    fig = plt.figure(figsize=(28, 28))
    ax = fig.add_subplot(projection='3d')

    # Plot the face points        
    ax.scatter(xs=all_coords[...,0],ys=all_coords[...,1],zs=all_coords[...,2],c="blue",s=0.5)
        
    # Label each of the axes
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    # Display the figure 
    plt.show()
    
    return None
    
#==============================================================================
# Define a function to find the face coordinates of a structured 3D mesh

def ExtractFaceCoords(volume,plot=True):
    
    # Extract the x-end cap faces (ABFE and DCGH)
    face_x1 = np.reshape(volume[ 0,:,:,:3],(-1,3))
    face_x2 = np.reshape(volume[-1,:,:,:3],(-1,3))
    
    # Extract the y-end cap faces (ADHE and BCGF)
    face_y1 = np.reshape(volume[:, 0,:,:3],(-1,3))
    face_y2 = np.reshape(volume[:,-1,:,:3],(-1,3))
    
    # Extract the z-end cap faces (ABCD and EFGH)
    face_z1 = np.reshape(volume[:,:, 0,:3],(-1,3))
    face_z2 = np.reshape(volume[:,:,-1,:3],(-1,3))
    
    # Add flattened face points to a list 
    face_coords = np.vstack((face_x1,face_x2,face_y1,face_y2,face_z1,face_z2))
    
    # Remove all duplicate coordinate points
    face_coords = np.unique([tuple(row) for row in face_coords],axis=0)
    
    if plot:
        # Create a figure with 3-dimensional axes    
        fig = plt.figure(figsize=(28, 28))
        ax = fig.add_subplot(projection='3d')

        # Plot the face points        
        ax.scatter(xs=face_coords[...,0],ys=face_coords[...,1],zs=face_coords[...,2],c="blue",s=0.5)
            
        # Label each of the axes
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
    
        # Display the figure 
        plt.show()
    
    return face_coords

#==============================================================================
# Define a function to plot the face surface points of a structured 3D mesh

def PlotFaces(volume,plot=True):
    
    # Extract the x-end cap faces (ABFE and DCGH)
    face_x1 = volume[ 0,:,:,:]
    face_x2 = volume[-1,:,:,:]
    
    # Extract the y-end cap faces (ADHE and BCGF)
    face_y1 = volume[:, 0,:,:]
    face_y2 = volume[:,-1,:,:]
    
    # Extract the z-end cap faces (ABCD and EFGH)
    face_z1 = volume[:,:, 0,:]
    face_z2 = volume[:,:,-1,:]
    
    if plot:
        # Create a figure with 3-dimensional axes    
        fig = plt.figure(figsize=(28, 28))
        ax = fig.add_subplot(projection='3d')
        
        for face in [face_x1,face_x2,face_y1,face_y2,face_z1,face_z2]:
        
            # Plot the end boundaries    
            ax.scatter(xs=face[...,0],ys=face[...,1],zs=face[...,2],c="blue",s=0.5)
            
        # Label each of the axes
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
    
        # Display the figure 
        plt.show()
        
    return None

#==============================================================================
# Define a function to find the edge coordinates of a structured 3D mesh

def ExtractEdgeCoords(volume,plot=True):
        
    # Extract the x-aligned edges (AD, BC, EH and FG)
    edge_x1 = np.reshape(volume[:, 0, 0,:3],(-1,3))
    edge_x2 = np.reshape(volume[:,-1, 0,:3],(-1,3))
    edge_x3 = np.reshape(volume[:, 0,-1,:3],(-1,3))
    edge_x4 = np.reshape(volume[:,-1,-1,:3],(-1,3))
    edge_x = np.vstack((edge_x1,edge_x2,edge_x3,edge_x4))
    
    # Extract the y-aligned edges (AB, DC, EF and HG)
    edge_y1 = np.reshape(volume[ 0,:, 0,:3],(-1,3))
    edge_y2 = np.reshape(volume[-1,:, 0,:3],(-1,3))
    edge_y3 = np.reshape(volume[ 0,:,-1,:3],(-1,3))
    edge_y4 = np.reshape(volume[-1,:,-1,:3],(-1,3))
    edge_y = np.vstack((edge_y1,edge_y2,edge_y3,edge_y4))
    
    # Extract the z-aligned edges (AE, DH, BF and CG)
    edge_z1 = np.reshape(volume[ 0, 0,:,:3],(-1,3))
    edge_z2 = np.reshape(volume[-1, 0,:,:3],(-1,3))
    edge_z3 = np.reshape(volume[ 0,-1,:,:3],(-1,3))
    edge_z4 = np.reshape(volume[-1,-1,:,:3],(-1,3))
    edge_z = np.vstack((edge_z1,edge_z2,edge_z3,edge_z4))

    # Add flattened face points to a list 
    edge_coords = np.vstack((edge_x,edge_y,edge_z))
    
    # Remove all duplicate coordinate points
    edge_coords = np.unique([tuple(row) for row in edge_coords],axis=0)
    
    if plot:
        # Create a figure with 3-dimensional axes    
        fig = plt.figure(figsize=(28, 28))
        ax = fig.add_subplot(projection='3d')
        
        # Plot the edge points
        ax.scatter(xs=edge_coords[...,0],ys=edge_coords[...,1],zs=edge_coords[...,2],c="blue",s=0.5) 
        
        # Label each of the axes
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
    
        # Display the figure 
        plt.show()
    
    return edge_coords

#==============================================================================
# Define a function to plot the edges boundaries of a structured 3D mesh

def PlotEdges(volume,plot=True):
    
    # Extract the x-aligned edges (AD, BC, EH and FG)
    edge_x1 = volume[:, 0, 0,:]
    edge_x2 = volume[:,-1, 0,:]
    edge_x3 = volume[:, 0,-1,:]
    edge_x4 = volume[:,-1,-1,:]
    
    # Extract the y-aligned edges (AB, DC, EF and HG)
    edge_y1 = volume[ 0,:, 0,:]
    edge_y2 = volume[-1,:, 0,:]
    edge_y3 = volume[ 0,:,-1,:]
    edge_y4 = volume[-1,:,-1,:]
    
    # Extract the z-aligned edges (AE, DH, BF and CG)
    edge_z1 = volume[ 0, 0,:,:]
    edge_z2 = volume[-1, 0,:,:]
    edge_z3 = volume[ 0,-1,:,:]
    edge_z4 = volume[-1,-1,:,:]
    
    if plot:
        # Create a figure with 3-dimensional axes    
        fig = plt.figure(figsize=(28, 28))
        ax = fig.add_subplot(projection='3d')
        
        # Plot the x-aligned edges  
        for edge in [edge_x1,edge_x2,edge_x3,edge_x4]:
            ax.scatter(xs=edge[...,0],ys=edge[...,1],zs=edge[...,2],c="blue",s=0.5) 
        
        # Plot the y-aligned edges  
        for edge in [edge_y1,edge_y2,edge_y3,edge_y4]:
            ax.scatter(xs=edge[...,0],ys=edge[...,1],zs=edge[...,2],c="blue",s=0.5) 
        
        # Plot the z-aligned edges  
        for edge in [edge_z1,edge_z2,edge_z3,edge_z4]:
            ax.scatter(xs=edge[...,0],ys=edge[...,1],zs=edge[...,2],c="blue",s=0.5) 
        
        # Label each of the axes
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
    
        # Display the figure 
        plt.show()
    return None

#==============================================================================
