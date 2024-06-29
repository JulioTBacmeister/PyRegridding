from mpi4py import MPI
import numpy as np

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Define the grid size and number of processes (assume size divides evenly for simplicity)
nx, ny = 100, 100  # Global grid size
px = int(np.sqrt(size))  # Processes in x direction
py = int(np.sqrt(size))  # Processes in y direction

local_nx = nx // px  # Local grid size for each process in x direction
local_ny = ny // py  # Local grid size for each process in y direction

# Create local grid with halo regions
local_grid = np.zeros((local_nx + 2, local_ny + 2))

# Initialize the local grid with some values
for i in range(1, local_nx + 1):
    for j in range(1, local_ny + 1):
        local_grid[i, j] = rank + 1  # Simple initialization for demonstration

# Exchange halo regions
def exchange_halos(comm, local_grid, local_nx, local_ny, rank, px, py):
    up = rank - px if rank >= px else MPI.PROC_NULL
    down = rank + px if rank < size - px else MPI.PROC_NULL
    left = rank - 1 if rank % px != 0 else MPI.PROC_NULL
    right = rank + 1 if (rank + 1) % px != 0 else MPI.PROC_NULL

    # Send up, receive down
    comm.Sendrecv(local_grid[1, 1:-1], dest=up, recvbuf=local_grid[-1, 1:-1], source=down)
    # Send down, receive up
    comm.Sendrecv(local_grid[-2, 1:-1], dest=down, recvbuf=local_grid[0, 1:-1], source=up)
    # Send left, receive right
    comm.Sendrecv(local_grid[1:-1, 1], dest=left, recvbuf=local_grid[1:-1, -1], source=right)
    # Send right, receive left
    comm.Sendrecv(local_grid[1:-1, -2], dest=right, recvbuf=local_grid[1:-1, 0], source=left)

# Perform halo exchange
exchange_halos(comm, local_grid, local_nx, local_ny, rank, px, py)

# Perform a simple finite difference calculation
finite_diff = np.zeros((local_nx, local_ny))
for i in range(1, local_nx + 1):
    for j in range(1, local_ny + 1):
        finite_diff[i-1, j-1] = (
            local_grid[i+1, j] + local_grid[i-1, j] +
            local_grid[i, j+1] + local_grid[i, j-1] -
            4 * local_grid[i, j]
        )

# Gather results for verification (optional, for demonstration purposes)
gathered = None
if rank == 0:
    gathered = np.empty((nx, ny))
comm.Gather(finite_diff, gathered, root=0)

if rank == 0:
    print("Gathered results:")
    print(gathered)

if __name__ == "__main__":
    print(f"Rank {rank} - local_grid with halos:\n", local_grid)
    print(f"Rank {rank} - finite_diff:\n", finite_diff)
