Introduction
============

What is python-sscha?
---------------------

python-sscha is a both a python library and a stand-alone program to simulate quantum and thermal fluctuations in solid systems.




Why do I need python-sscha?
---------------------------


If you are simulating transport or thermal properties of materials, phase diagrams, or phonon-related properties of materials, than you need python-sscha.
This is a package that enables you to include the effect of both thermal and quantum phonon fluctuations into your *ab initio* simulations.

The method used by this package is the  Stochastic self-consistent Harmonic Approximation (SSCHA). This is a full-quantum method that optimizes the nuclear wave-function (or density matrix at finite temperature) to minimize the free energy.
In this way you can simulatie highly anharmonic systems, like those close to a second order phase transition (as charge density waves and thermoelectric materials). 
Despite the full quantum and thermal nature of the algorithm, the overall computational cost is the comparable to standard classical molecular dynamics. Since the algorithm correctly exploits the symmetries of the crystal, in highly symmetric cases it is also much cheaper. 

python-sscha comes both as a python library that can be runned inside your own workflows and as a stand-alone software, initialized by input scripts with the same syntax as Quantum ESPRESSO.

It can be coupled with any *ab initio* engine for force and energy calculations. It can interact through the Atomic Simulation Environment (ASE), and has an implemented interface for automatic submission of jobs in a remote cluster.

Moreover, it is quite easy to use, a standard input is highly human readable and less than 10 lines!
So, what are you waiting? Download and install python-sscha, and start enjoing the Tutorials!


