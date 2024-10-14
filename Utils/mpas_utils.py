import numpy as np

def uv_cell_to_edge(uZonal, uMerid, nlev, lonEdge, latEdge, lonCell, latCell, edgeNormalVectors, cellsOnEdge):
    nCells = len(lonCell)
    nEdges = len(lonEdge)

    # Initialize east and north vectors
    east = np.zeros((3, nCells))
    north = np.zeros((3, nCells))

    # Calculate east and north vectors for each cell
    for iCell in range(nCells):
        # East vector
        east[0, iCell] = -np.sin(lonCell[iCell])
        east[1, iCell] = np.cos(lonCell[iCell])
        east[2, iCell] = 0.0

        # Normalize east vector
        east[:, iCell] /= np.linalg.norm(east[:, iCell])

        # North vector
        north[0, iCell] = -np.cos(lonCell[iCell]) * np.sin(latCell[iCell])
        north[1, iCell] = -np.sin(lonCell[iCell]) * np.sin(latCell[iCell])
        north[2, iCell] = np.cos(latCell[iCell])

        # Normalize north vector
        north[:, iCell] /= np.linalg.norm(north[:, iCell])

    # Initialize uNormal array
    uNormal = np.zeros((nlev, nEdges), dtype=uZonal.dtype)

    # Calculate uNormal for each edge
    for iEdge in range(nEdges):
        if iEdge % 10000 == 0:
            print(f"... UV_CELL_TO_EDGE: {100. * iEdge / (nEdges - 1):.2f}%")

        cell1 = cellsOnEdge[iEdge, 0] - 1
        cell2 = cellsOnEdge[iEdge, 1] - 1

        # Calculate uNormal for edge
        uNormal[:, iEdge] = (uZonal[:, cell1] * 0.5 * (edgeNormalVectors[iEdge, 0] * east[0, cell1] +
                                                       edgeNormalVectors[iEdge, 1] * east[1, cell1] +
                                                       edgeNormalVectors[iEdge, 2] * east[2, cell1]) +
                             uMerid[:, cell1] * 0.5 * (edgeNormalVectors[iEdge, 0] * north[0, cell1] +
                                                       edgeNormalVectors[iEdge, 1] * north[1, cell1] +
                                                       edgeNormalVectors[iEdge, 2] * north[2, cell1]) +
                             uZonal[:, cell2] * 0.5 * (edgeNormalVectors[iEdge, 0] * east[0, cell2] +
                                                       edgeNormalVectors[iEdge, 1] * east[1, cell2] +
                                                       edgeNormalVectors[iEdge, 2] * east[2, cell2]) +
                             uMerid[:, cell2] * 0.5 * (edgeNormalVectors[iEdge, 0] * north[0, cell2] +
                                                       edgeNormalVectors[iEdge, 1] * north[1, cell2] +
                                                       edgeNormalVectors[iEdge, 2] * north[2, cell2]))

    return uNormal
