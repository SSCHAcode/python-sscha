Introduction
============

What is python-sscha?
---------------------

python-sscha is both a python library and a stand-alone program to simulate quantum and thermal fluctuations in solid systems.




Why do I need python-sscha?
---------------------------


If you are simulating transport or thermal properties of materials, phase diagrams, or phonon-related properties of materials, you need python-sscha.
It is a package that enables you to include the effect of both thermal and quantum phonon fluctuations into your *ab initio* simulations.

The method used by this package is the  Stochastic self-consistent Harmonic Approximation (SSCHA). The SSCHA is a full-quantum method that optimizes the nuclear wave-function (or density matrix at finite temperature) to minimize the free energy.
In this way, you can simulate highly anharmonic systems, like those close to a second-order phase transition (as charge density waves and thermoelectric materials). 
Despite the full quantum and thermal nature of the algorithm, the overall computational cost is comparable to standard classical molecular dynamics. Since the algorithm correctly exploits the symmetries of the crystal, it is also much cheaper. 

python-sscha comes both as a python library that can be run inside your workflows and as stand-alone software, initialized by input scripts with the same syntax as Quantum ESPRESSO.

You can couple it with any *ab initio* engine for force and energy calculations. It can interact through the Atomic Simulation Environment (ASE) and has an implemented interface for automatic submission of jobs in a remote cluster.

Moreover, it is easy to use, with short input files highly human-readable.
What are you waiting for? Download and install python-sscha, and start enjoying the Tutorials!


