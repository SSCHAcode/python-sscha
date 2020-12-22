from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt

import ase, ase.units, ase.atoms
import ase.calculators.calculator as calc

import cellconstructor as CC
import cellconstructor.Structure


class ToyModel2D(calc.Calculator):
    """
    TOY MODEL 2D
    ============

    This is a 2D toy model of the shape:
    
    .. math::
    
        a \cdot |R|^4 + b |R|^3 - c |R|^2 + d x + e z^2

    where :math:`|R| = \sqrt{x^2 + y^2}`. 
    The structure on input should have 2 atoms. 
    The x, y and z variables are intended as the distance between the two atoms, 
    while there is no force on the center of mass (assuring the presence of acoustic modes). 
    The z direction is perfectly harmonic, and decoupled by the other modes.
    """

    def __init__(self, *args, **kwargs):
        calc.Calculator.__init__(self, *args, **kwargs)
        self.implemented_properties = ["energy", "forces"]

        # conversion units
        convert = 10

        # Define the parameters
        self.a = 3 * convert**4
        self.b = 0.5 * convert ** 3
        self.c = -3 * convert ** 2
        self.d = 1 / 3. * convert
        self.e = 1 * convert ** 2


    def calculate(self, atoms = None, *args, **kwargs):
        calc.Calculator.calculate(self, atoms, *args, **kwargs)

        # Get the vector of distance
        pos = atoms.get_positions()

        nat, dim = pos.shape

        if nat < 2:
            raise ValueError("Error, this calculator requires at least two atoms")
        
        r_vect = pos[1,:] - pos[0,:]

        r_xy = np.sqrt(np.sum(r_vect[:2]**2))
        x = r_vect[0]
        y = r_vect[1]
        z = r_vect[2]
        
        energy = self.a * r_xy**4 + self.b * r_xy**3 + self.c * r_xy**2 + self.d * x + self.e * z**2

        force = np.zeros(3, dtype = np.double)
        
        r_xy_proj = r_vect.copy() / r_xy
        r_xy_proj[2] = 0

        # Get the derivative
        force[:] = (4 * self.a * r_xy**3 + 3 * self.b * r_xy**2 + 2 * self.c * r_xy) * r_xy_proj
        force[0] += self.d
        force[2] += 2 * self.e * z

        force *= -1

        total_force = pos.copy()
        total_force[0, :] = - force
        total_force[1, :] = force

        self.results = {"energy" : energy, "forces" : total_force}
        

        
        
def test_toy_model():
    calc = ToyModel2D()

    # Define a simple system
    atoms = ase.atoms.Atoms("H2", [(0,0,0), (0,0,0)])
    atoms.calc = calc

    x = []
    delta = 0.0001
    energy = []
    force = []
    for i in range(0, 1000):
        # Update the positions
        pos = atoms.get_positions().copy()
        pos[1,0] += delta
        atoms.set_positions(pos)

        # Get energy and force
        energy.append(atoms.get_total_energy())
        force.append(atoms.get_forces()[1,0])
        x.append(pos[1,0])
        
    plt.plot(x, energy, label = "energy")
    plt.legend()
    plt.tight_layout()
    
    plt.figure()
    plt.plot(x, force, label = "model force")
    plt.plot(x, -np.gradient(energy) / delta, label = "numerical force")
    plt.legend()
    plt.tight_layout()
    plt.show()



if __name__ == "__main__":
    
    test_toy_model()
