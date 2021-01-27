"""
Here we set up a calculator for the CP2K water Toy Model
"""

from ase.calculators.calculator import FileIOCalculator, all_changes
from ase.io import write
from ase import Atoms
from ase.units import Rydberg, Bohr, GPa
import numpy as np
import os

class CP2K_water_calculator(FileIOCalculator):

    # The default command to run CP2K
    command = "cp2k -i PREFIX.in -o PREFIX.out"
    implemented_properties = ['energy', 'forces', 'stress']


    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='cp2k_toy_water', atoms=None, **kwargs):

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

    def write_input(self, atoms, properties = None, system_charges = None):
        """
        """

        # Initialize the directories and setup correctly the working path
        FileIOCalculator.write_input(self, atoms, properties, system_charges)

        # Write the input
        input_dist = []

        # Prepare the toy model info for CP2K
        input_dist.append(__TOY_MODEL_INFO__)
        
        # Insert the cell
        cell = atoms.get_cell()
        STR1 = "A [angstrom]  %.16f  %.16f  %.16f\n" % (cell[0, 0], cell[0,1], cell[0,2])
        STR2 = "B [angstrom]  %.16f  %.16f  %.16f\n" % (cell[1, 0], cell[1,1], cell[1,2])
        STR3 = "C [angstrom]  %.16f  %.16f  %.16f\n" % (cell[2, 0], cell[2,1], cell[2,2])

        input_dist.append(STR1)
        input_dist.append(STR2)
        input_dist.append(STR3)

        input_dist.append(__INPUT_SECOND_PART__)

        # Add the coordinate file
        write("%s.xyz" % self.label, atoms)

        STR_ATM = "     COORD_FILE_NAME   %s.xyz\n" % self.label
        input_dist.append(STR_ATM)

        PRINT_INFO = """&END
 &END SUBSYS
  &PRINT
  &FORCES ON
  ADD_LAST NUMERIC
  COMMON_ITERATION_LEVELS 1
  FILENAME forces
  &END FORCES
  &STRESS_TENSOR ON
   ADD_LAST NUMERIC
   COMMON_ITERATION_LEVELS 1
   FILENAME stress
  &END STRESS_TENSOR
 &END PRINT
 """
        input_dist.append(PRINT_INFO)

        input_dist.append(__INPUT_END__)

        # Write the input in the PREFIX.in file specified in the command section
        f = open("%s.in" % self.label, "w")
        f.writelines(input_dist)
        f.close()



    def read(self, label):
        # Set the label
        FileIOCalculator.read(self, label)

        # We can read the output file from cp2k
        outfile = file("%s.out" % label, "r")
        outlines = [line.strip() for line in  outfile.readlines() ]
        outfile.close()

        print "READING!"

        unit_cell = np.zeros((3,3), dtype = np.float64)
        reading_atoms = False
        nat = 0
        atm_coords = []
        atm_type = []
        __start__ = 0
        for i, line in  enumerate(outlines):
            data_line = line.strip()

            if len(data_line) == 0:
                continue

            if data_line[0] == "CELL|" and data_line[1] == "Vector" :
                index = {"a" : 0, "b" : 1, "c" : 2}
                v_i = index[data_line[2]]
                unit_cell[v_i, 0] = float( data_line[4])
                unit_cell[v_i, 1] = float( data_line[5])
                unit_cell[v_i, 2] = float( data_line[6])
            
            if "ATOMIC COORDINATES" in line:
                reading_atoms = True
            
            if reading_atoms:
                if len(data_line) != 8:
                    reading_atoms = False
                    continue

                if data_line[0] == "Atom":
                    continue
                
                atm_type_conv = {"H": "H", "D": "H", "O":"O"} # Avoid deuterium
                nat += 1
                atm_coords.append( [float(data_line[3]), float(data_line[4]), float(data_line[5])])
                atm_type.append( atm_type_conv[data_line[2]])

        # Generate the ASE ATOM structure
        self.atoms = Atoms(atm_type, atm_coords, cell=unit_cell)
        self.read_results()

        print "Test results:"
        print "energy:", self.results["energy"]
        print "forces:", self.results["forces"]
        print "stress:", self.results["stress"]


    def read_results(self):
        """
        Read from the output forces and stress file the results
        """

        outfile = file("%s.out" % self.label, "r")
        outlines = [line.strip().split()[-1] for line in  outfile.readlines() if "ENERGY|" in line]
        outfile.close()

        energy = float(outlines[-1]) * 2 * Rydberg

        # Now read the forces and the stress
        fforc = file("TM-forces-1.xyz")
        f_lines = [l.strip() for l in fforc.readlines()]
        fforc.close()

        forces = []
        for line in f_lines:
            data = line.split()
            if len(data) == 6:
                fx = [float(data[3]), float(data[4]), float(data[5])]
                fx = [x * 2 * Rydberg / Bohr for x in fx]
                forces.append(fx)
        
        forces = np.array(forces, dtype = np.float64)

        # Read the stress
        fstress = file("TM-stress-1.stress_tensor")
        s_lines = [l.strip() for l in fstress.readlines()]
        fstress.close()

        stress = np.zeros((3,3), dtype = np.float64)

        i = 0
        for line in s_lines:
            data = line.split()
            if len(data) != 4:
                continue
            if data[0] != "X" :
                continue

            stress[i,0] = float(data[1])
            stress[i,1] = float(data[2])
            stress[i, 2] = float(data[3])
            i += 1
        voigth_stress = [stress[0,0], stress[1,1], stress[2,2], stress[1,2], stress[0,2], stress[0,1]]
        voigth_stress = np.array(voigth_stress, dtype = np.float64) * GPa
        self.results = {"energy" : energy, "forces" : forces, "stress" : voigth_stress}
        #print "RESULTS!!! ", self.results


    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        # Remove, if any the previous calculation file
        if os.path.exists("%s/%s.out" % (self.directory, self.label)):
            os.remove("%s/%s.out" % (self.directory, self.label))
        if os.path.exists("%s/TM-forces-1.xyz" % self.directory):
            os.remove("%s/TM-forces-1.xyz" % self.directory)
        if os.path.exists("%s/TM-stress-1.stress_tensor" % self.directory):
            os.remove("%s/TM-stress-1.stress_tensor" % self.directory)
        
        FileIOCalculator.calculate(self, atoms, properties, system_changes)





                

                
                









    


__TOY_MODEL_INFO__ = """
&FORCE_EVAL
 STRESS_TENSOR ANALYTICAL
 METHOD FIST
 &MM
   &FORCEFIELD
     &SPLINE
       EMAX_SPLINE 1.e8
     &END
     &CHARGE
       ATOM O
       CHARGE -0.8476
     &END CHARGE
     &CHARGE
       ATOM H
       CHARGE 0.4238
     &END CHARGE
     &CHARGE
       ATOM D
       CHARGE 0.4238
     &END CHARGE
     &BOND
   ATOMS O  H
       K [angstrom^-2kcalmol]  1214.4  -4166.0  7410.3
       R0  [angstrom]            0.9419
       KIND   QUARTIC
     &END BOND
     &BOND
   ATOMS O  D
       K [angstrom^-2kcalmol]  1214.4  -4166.0  7410.3
       R0  [angstrom]            0.9419
       KIND   QUARTIC
     &END BOND
    &BEND
     ATOMS  H  O  H
     K [rad^-2kcalmol] 75.90
     THETA0 [deg] 112.0
     KIND HARMONIC
    &END BEND
    &BEND
     ATOMS  H  O  D
     K [rad^-2kcalmol] 75.90
     THETA0 [deg] 112.0
     KIND HARMONIC
    &END BEND
    &NONBONDED
       &LENNARD-JONES
         ATOMS O O
     EPSILON [kcalmol] 0.1554
     SIGMA [angstrom] 3.1655
     RCUT [angstrom] 9.85
       &END LENNARD-JONES
   &LENNARD-JONES
         atoms H O
     EPSILON [hartree] 0.d0
     SIGMA [bohr] 0.d0
     RCUT [angstrom] 9.85
       &END LENNARD-JONES
    &LENNARD-JONES
         atoms H H
     EPSILON [hartree] 0.d0
     SIGMA [bohr] 0.d0
     RCUT [angstrom] 9.85
!      RMIN [angstrom] 1.0
       &END LENNARD-JONES
   &LENNARD-JONES
         atoms D O
     EPSILON [hartree] 0.d0
     SIGMA [bohr] 0.d0
     RCUT [angstrom] 9.85
       &END LENNARD-JONES
    &LENNARD-JONES
         atoms D H
     EPSILON [hartree] 0.d0
     SIGMA [bohr] 0.d0
     RCUT [angstrom] 9.85
!      RMIN [angstrom] 1.0
       &END LENNARD-JONES
    &LENNARD-JONES
         atoms D D
     EPSILON [hartree] 0.d0
     SIGMA [bohr] 0.d0
     RCUT [angstrom] 9.85
!      RMIN [angstrom] 1.0
       &END LENNARD-JONES
    &END NONBONDED
   &END FORCEFIELD
   &POISSON
     &EWALD
       EWALD_TYPE EWALD
       ALPHA .40
       EWALD_ACCURACY 1.0E-2
       GMAX 11
     &END EWALD
   &END POISSON
 &END MM
 &SUBSYS
   &CELL
        """

__INPUT_SECOND_PART__ = """
   &END CELL
   &KIND D
     ELEMENT H
!      ATOMIC_MASS 1.0995
   &END KIND
   &TOPOLOGY
     CONNECTIVITY MOL_SET
     &MOL_SET
!        &MOLECULE
!          NMOL 1
!          CONN_FILE_NAME topo_hdo.psf
!        &END
       &MOLECULE
         NMOL 4
         CONN_FILE_NAME topology_fist_WAT.psf
       &END
     &END
     COORDINATE XYZ
"""

__INPUT_END__ = """&END FORCE_EVAL
&GLOBAL
 PROJECT TM
! RUN_TYPE GEO_OPT
 RUN_TYPE ENERGY_FORCE
 PRINT_LEVEL HIGH
&END GLOBAL
&MOTION
  &GEO_OPT
   OPTIMIZER CG
  &END GEO_OPT
   &PRINT
      &TRAJECTORY ON
         FORMAT XMOL
         FILENAME ./traj.xyz
         &EACH
           MD 1
!           JUST_ENERGY 1
         &END EACH
      &END TRAJECTORY
   &END PRINT
&END MOTION
"""