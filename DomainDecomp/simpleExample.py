from mpi4py import MPI
import numpy as np

def initialize_grid(nx, ny):
    """
    Initialize a 2D grid with random values
    """
    return np.random.random((nx, ny))

def process_subdomain(subdomain):
    """
    Perform some computation on the subdomain
    For simplicity, we'll just sum the values
    """
    return np.sum(subdomain)

def main():
    # Initialize MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Define the global grid size
    nx, ny = 100, 100

    # Determine the local grid size for each process
    local_nx = nx // size
    local_ny = ny

    # Initialize the grid on the root process
    if rank == 0:
        grid = initialize_grid(nx, ny)
    else:
        grid = None

    # Create a buffer for the local subdomain
    local_subdomain = np.empty((local_nx, local_ny), dtype='d')

    # Scatter the global grid to all processes
    comm.Scatter(grid, local_subdomain, root=0)

    # Process the local subdomain
    local_result = process_subdomain(local_subdomain)

    # Gather the results from all processes
    results = comm.gather(local_result, root=0)

    # Combine the results on the root process
    if rank == 0:
        total_result = np.sum(results)
        print(f"Total result: {total_result}")

if __name__ == '__main__':
    main()
