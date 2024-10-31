import numpy as np
import xarray as xr
from datetime import datetime
import sys
import copy
import esmpy as E
"""
--------------------------------------------
Strutcture of sst_HadISST 360x180 mesh file
---------------------------------------------
Dimensions:
origGridRank: 2 nodeCount: 65160  coordDim: 2  elementCount: 64800  maxNodePElement: 4

Data variables:
origGridDims (origGridRank) int32
nodeCoords (nodeCount, coordDim) float64
elementConn (elementCount, maxNodePElement) float64
numElementConn (elementCount) int32
centerCoords (elementCount, coordDim) float64
elementArea (elementCount) float64
elementMask (elementCount) int32
"""
def nxny_to_mesh( nx=None, ny=None, mesh_file_out=None, **kwargs ):
    ########################################################################
    # This code assumes a glolbal nx (lon, periodic) by ny (lat) grid
    ########################################################################
    # Define the latitude and longitude arrays
    latitudes = np.linspace(-90, 90, ny + 1)  # ny+1 grid points in latitude
    longitudes = np.linspace(0, 360, nx, endpoint=False)  # nx grid points in longitude, periodic in lon
    delta_lon = (360.) / nx  # Longitude cell width in radians
    center_lons = longitudes + 0.5*delta_lon 
    center_lats = 0.5*( latitudes[1:]+latitudes[0:-1] )

    #return latitudes,longitudes,center_lats,center_lons
    
    # Initialize the unique_nodeCoords array to hold (nx)*(ny+1) unique node coordinates
    unique_nodeCoords = np.zeros(  ((ny + 1) * nx, 2))
    
    # Set up a simple loop to fill in the coordinates
    index = 0
    for i in range(ny + 1):  # Loop over latitudes
        for j in range(nx):  # Loop over longitudes (nx points)
            unique_nodeCoords[index, 0] = longitudes[j]  # Longitude (x-axis)
            unique_nodeCoords[index, 1] = latitudes[i]   # Latitude (y-axis)
            index += 1
    
    # Now, unique_nodeCoords should have the shape (1038240, 2)
    print( unique_nodeCoords[0:2,:] )
    
    
    # Initialize the elementConn array: (nx * ny, 4) where 4 represents the 4 corners of each cell
    elementConn = np.zeros((nx * ny, 4), dtype=int)
    
    # Loop through each grid cell
    index = 0
    # Might need this for frotran indexing
    add_one_for_fortran = 1
    for i in range(ny):  # Loop over latitudes (no need for the last row of latitudes since it's a boundary)
        for j in range(nx):  # Loop over longitudes
            # Current grid cell corners
            SW = i * nx + j  + add_one_for_fortran # South-West corner
            SE = i * nx + (j + 1) % nx   + add_one_for_fortran # South-East (wrap-around if at the last column)
            NW = (i + 1) * nx + j   + add_one_for_fortran # North-West
            NE = (i + 1) * nx + (j + 1) % nx   + add_one_for_fortran # North-East (wrap-around)
    
            # Store the corners for the current element (grid cell)
            elementConn[index, :] = [SW, SE, NE, NW]
            index += 1
    
    # Now elementConn should have the shape (nx * ny, 4) and contain the correct connectivity


    # At this pint everything is there for basic mesh file
    # Initialize the ESMF mesh object

    
    # Define nodeCount, nodeOwners, and nodeMask
    nodeCount = unique_nodeCoords.shape[0]
    nodeOwners = np.zeros(nodeCount, dtype=int)  # All nodes owned by processor 0
    nodeMask = np.zeros(nodeCount, dtype=int)    # No masking
    nodeIds = np.arange(nodeCount)
    
    # Add the nodes (unique corner coordinates)
    # mesh.add_nodes(nodeCount, nodeIds, unique_nodeCoords, nodeOwners   )  #, nodeMask)
    
    # Calculate the number of elements (grid cells)
    num_elements = nx * ny  # Total number of quadrilateral cells in the grid
    
    # Define element_types and element_mask for the elements
    element_types = np.full(num_elements, E.MeshElemType.QUAD)  # All quadrilateral elements
    element_mask = np.zeros(num_elements, dtype=int)  + 1             # No masking
    element_ids = np.arange(num_elements)  
    num_element_conn = np.zeros(num_elements, dtype=int) + 4       # Always 4 nodes per element

    
    # Convert latitudes and longitudes to radians
    lat_radians = np.radians(latitudes)
    delta_lambda = (2 * np.pi) / nx  # Longitude cell width in radians
    
    # Calculate the area of each grid cell
    elementArea = np.zeros(nx * ny)
    
    index = 0
    for i in range(ny):
        # Latitude boundaries of the cell (radians)
        phi_south = lat_radians[i]
        phi_north = lat_radians[i + 1]
        
        # Area of the grid cell
        area = delta_lambda * (np.sin(phi_north) - np.sin(phi_south))
        
        # Fill the elementArea array for each longitude cell
        for j in range(nx):
            elementArea[index] = area
            index += 1
    # Now elementArea contains the area of each element in steradians

    # Calculate the center coord of each grid cell
    elementCoords = np.zeros( (nx * ny, 2) )
    
    index = 0
    for i in range(ny):
        # Latitude boundaries of the cell (radians)
        for j in range(nx):
            elementCoords[index, 0] = center_lons[j]  # Longitude (x-axis)
            elementCoords[index, 1] = center_lats[i]   # Latitude (y-axis)
            index += 1
    # Now elementCoords contains the coord each element in degrees


    """
    --------------------------------------------
    From ESMF documentation
    ---------------------------------------------
    rc = _ESMF.ESMC_MeshAddElements(mesh.struct.ptr, lec,
                                    elementIds, elementTypes,
                                    elementConn, elementMask, elementArea,
                                    elementCoords)
    """

    # generate output dataset
    dso = xr.Dataset()    
    dso['origGridDims'] = xr.DataArray(np.array([nx, ny], dtype=np.int32), 
                                    dims=('origGridRank',)) 
    dso.origGridDims.encoding = {'dtype': np.int32}


    dso['nodeCoords'] = xr.DataArray( unique_nodeCoords ,
                                    dims=('nodeCount','coordDim',) , 
                                    attrs={'units': 'degrees'} ) 
    dso.nodeCoords.encoding = {'dtype': np.float64}

    dso['elementConn'] = xr.DataArray( elementConn ,
                                    dims=('elementCount','maxNodePElement',)  , 
                                    attrs={'long_name': 'Node indices for element connectivity (Fortran)'}  ) 
    #dso.elementConn.encoding = {'dtype': np.float64}

    dso['numElementConn'] = xr.DataArray( num_element_conn.astype(np.int32) ,
                                    dims=('elementCount',)) 
    #dso.numElementConn.encoding = {'dtype': np.int32}

    dso['centerCoords'] = xr.DataArray( elementCoords ,
                                    dims=('elementCount','coordDim',) , 
                                    attrs={'units': 'degrees'}  ) 
    #dso.numElementConn.encoding = {'dtype': np.int32}

    dso['elementArea'] = xr.DataArray( elementArea ,
                                    dims=('elementCount',)) 
    #dso.numElementConn.encoding = {'dtype': np.int32}

    dso['elementMask'] = xr.DataArray( element_mask.astype(np.int32) ,
                                    dims=('elementCount',)) 
    #dso.numElementConn.encoding = {'dtype': np.int32}


    return dso
    
def file_to_mesh( mesh_file_in=None, **kwargs ):
    """
    Contents of mesh netcdf file:
    Data variables:
        origGridDims    (origGridRank) int32 8B ...
        nodeCoords      (nodeCount, coordDim) float64 1MB ...
        elementConn     (elementCount, maxNodePElement) int64 2MB ...
        numElementConn  (elementCount) int32 259kB ...
        centerCoords    (elementCount, coordDim) float64 1MB ...
        elementArea     (elementCount) float64 518kB ...
        elementMask     (elementCount) int32 259kB ...
    """

    Mx = xr.open_dataset( mesh_file_in )
    #print(Mx)
    nodeCoords = Mx.nodeCoords.values
    elementConn = Mx.elementConn.values
    numElementConn = Mx.numElementConn.values
    centerCoords = Mx.centerCoords.values
    elementArea = Mx.elementArea.values
    elementMask = Mx.elementMask.values
    

    # Define nodeCount, nodeOwners, and nodeMask
    nodeCount = nodeCoords.shape[0]
    nodeOwners = np.zeros(nodeCount, dtype=int)  # All nodes owned by processor 0
    nodeMask = np.zeros(nodeCount, dtype=int)    # No masking
    nodeIds = np.arange(nodeCount)

    # Calculate the number of elements (grid cells)
    num_elements = elementConn.shape[0]

    print( f"num nodes={nodeCount} , num_elements={num_elements}" )

    # Define element_types and element_mask for the elements
    element_types = np.full(num_elements, E.MeshElemType.QUAD)  # All quadrilateral elements
    element_mask = np.zeros(num_elements, dtype=int)  + 1             # No masking
    element_ids = np.arange(num_elements)  
    
    mesh = E.Mesh(parametric_dim=2, spatial_dim=2)
    
    
    #add_nodes(node_count, node_ids, node_coords, node_owners, node_mask=None)
    #Add nodes to a Mesh, this must be done before adding elements.
    # Add the nodes (unique corner coordinates)
    mesh.add_nodes(nodeCount, nodeIds, nodeCoords, nodeOwners   )  #, nodeMask)

    # Reset Python indexing
    add_one_for_fortran = 1
    elementConn = elementConn - add_one_for_fortran
    ## add_elements(element_count, element_ids, element_types, element_conn, element_mask=None, element_area=None, element_coords=None)    
    # Add elements (cells) to the mesh using the connectivity information, element types, and mask
    mesh.add_elements(num_elements, element_ids, element_types, elementConn, 
                     element_mask=element_mask , element_area=elementArea , element_coords=centerCoords ) # , element_num_nodes=numElementConn) 

    if ("debug_output" in kwargs):
        return mesh, elementConn,centerCoords,nodeCoords,elementArea,elementMask
    else:
        return mesh


    