from __future__ import print_function
"""
This files contains a setup utility to manage the parallelization with different
modules.
"""

import numpy as np
import time
import os
import sys

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

def _force_stdout_blocking():
    """
    Make sure the standard output is in blocking mode.

    MPI launchers (mpirun, srun, ...) frequently attach the standard output
    of the processes to a pipe that is opened in *non-blocking* mode. When a
    large message fills the pipe buffer, a write on a non-blocking descriptor
    raises ``BlockingIOError`` ([Errno 11]) instead of waiting for the buffer
    to drain. This used to crash a parallel minimization while printing the
    table of imaginary frequencies (see issue #196).

    Restoring the blocking mode makes the write wait for the reader to consume
    the buffer, which is the expected behaviour for log output.
    """
    if not hasattr(os, "set_blocking"):
        return
    try:
        fd = sys.stdout.fileno()
    except (AttributeError, OSError, ValueError):
        # stdout has been replaced by an object without a real descriptor
        return
    try:
        os.set_blocking(fd, True)
    except (OSError, ValueError):
        pass


def pprint(*argv, **kwargs):
    """
    PARALLEL PRINTING
    =================

    This will print on stdout only once in parallel execution of the code.

    It is robust against a standard output opened in non-blocking mode, which
    is a frequent source of ``BlockingIOError`` when running under MPI (see
    issue #196): a log line must never abort a running calculation.
    """
    #print("pypar:", __PYPAR__)
    #print("mpi4py:", __MPI4PY__)
    if not am_i_the_master():
        return

    _force_stdout_blocking()
    try:
        print(*argv, **kwargs)
    except BlockingIOError:
        # The blocking mode could not be enforced (e.g. os.set_blocking is
        # not available or failed). Never let a print crash the calculation.
        pass
