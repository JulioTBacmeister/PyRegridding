#!/usr/bin/env python
# Import packages 
from mpi4py import MPI
import numpy as np
import os

def main():
    # Initialize MPI

    print( "This is stoopyExample" )
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    print("Rank, size",rank,size)
    
    #nworkers = cpu_count()
    nworkers = len(os.sched_getaffinity(0))
    print("nworkers ", nworkers)
 
if __name__ == '__main__':
    main()
