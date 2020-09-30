from __future__ import print_function
"""
This files contains a setup utility to manage the parallelization with different
modules.
"""

import numpy as np 
import time

# Supports both pypar and mpi4py
__PYPAR__ = False 
__MPI4PY__ = False
try: 
    import mpi4py, mpi4py.MPI
    __MPI4PY__ = True
except:
    try:
        import pypar
        __PYPAR__ = True  
    except:
        pass

AllParallel = [__PYPAR__, __MPI4PY__]


def is_parallel():
    """
    Returns True if the MPI parallelization is active,
    False otherwise
    """
    if True in AllParallel:
        return True 
    return False

def am_i_the_master():
    if __PYPAR__:
        if pypar.rank() == 0:
            return True
        return False  
    elif __MPI4PY__:
        comm = mpi4py.MPI.COMM_WORLD
        if comm.rank == 0:
            return True 
        else:
            return False
    else:
        return True 

def pprint(*argv):
    """
    PARALLEL PRINTING
    =================

    This will print on stdout only once in parallel execution of the code
    """
    #print("pypar:", __PYPAR__)
    #print("mpi4py:", __MPI4PY__)
    if am_i_the_master():
        print(*argv)

