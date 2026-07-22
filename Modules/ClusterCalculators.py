# -*- coding: utf-8 -*-
"""Calculators for the sscha cluster backends (Axis B of the cluster design).

Two families:

* FileIOCalculator: calculators driven through input/output files, consumed
  by the scheduler-based clusters (sscha.Cluster.Cluster,
  sscha.LocalCluster.LocalCluster). One subclass per simulation code
  (EspressoCalculator today; VaspCalculator, AbinitCalculator,
  Cp2kCalculator, GaussianCalculator, ... in the future).

* DirectCalculator: calculators computed in-process, consumed by
  sscha.BaseCluster.DirectCluster.

For any code that has an ASE interface, ASEFileCalculator (file-based) and
ASEDirectCalculator (in-process) already provide full support with no
dedicated subclass: e.g. ASEFileCalculator(ase.calculators.vasp.Vasp(...)).
"""
import copy
import os
import pickle
import sys

import cellconstructor as CC
import cellconstructor.calculators


# ---------------------------------------------------------------------- #
#  Generic contract                                                       #
# ---------------------------------------------------------------------- #

class ClusterCalculator(object):
    """Documentation-level base class for cluster calculators.

    Any calculator consumed by a sscha cluster must provide:

        copy() -> an independent instance (one per worker thread)

    and produce result dicts with the contract documented in
    BaseCluster.compute_jobarray ("energy" [eV], "forces" [eV/Ang],
    "stress" [ASE Voigt, -eV/Ang^3], "structure" [CC.Structure], extras).
    """


# ---------------------------------------------------------------------- #
#  File-based family (scheduler clusters)                                 #
# ---------------------------------------------------------------------- #

class FileIOCalculator(CC.calculators.FileIOCalculator, ClusterCalculator):
    """Mid-level interface for file-based calculators.

    Class attributes consumed by the scheduler machinery:

        input_extension / output_extension : str
            Extensions of the per-label input/output files.
        extra_input_files : list of str
            Additional files (relative to the local workdir) to ship to the
            cluster together with the inputs (e.g. "calculator.pkl").
    """

    input_extension = ".pwi"
    output_extension = ".pwo"
    extra_input_files = []

    def get_execution_command(self):
        """The command template run on the cluster, with the PREFIX
        placeholder where the calculation label must be inserted, e.g.:

            "pw.x -npool 4 -i PREFIX.pwi > PREFIX.pwo"

        Assign it to `cluster.binary` (or use it in the namelist).
        """
        raise NotImplementedError


class EspressoCalculator(CC.calculators.Espresso, FileIOCalculator):
    """Quantum ESPRESSO (pw.x) calculator for the sscha clusters.

    Thin formalization of cellconstructor.calculators.Espresso, which
    already satisfies the interface. Using this class (instead of the CC
    one) only makes the extensions and the execution command explicit.
    """

    input_extension = ".pwi"
    output_extension = ".pwo"
    extra_input_files = []

    def get_execution_command(self, n_pool=1):
        return "pw.x -npool {} -i PREFIX.pwi > PREFIX.pwo".format(n_pool)

    def copy(self):
        """Return an identical instance of THIS class (per-thread copies)."""
        return EspressoCalculator(self.input_data, self.pseudopotentials,
                                  self.masses, self.command, self.kpts, self.koffset)


class ASEFileCalculator(FileIOCalculator):
    """
    RUN ANY ASE CALCULATOR THROUGH A SCHEDULER-BASED SSCHA CLUSTER
    ==============================================================

    Bridge between a plain ASE calculator (EMT, LJ, VASP, CP2K, ...) and the
    file-based `Cluster` submission:

    - the ASE calculator is pickled ONCE into the working directory
      (`calculator.pkl`, listed in `extra_input_files` so it is shipped to
      the cluster together with the inputs);
    - each input file is the ASE Atoms serialized with ase.io.jsonio;
    - the execution command runs `python -m sscha.ASEClusterRunner`, which
      loads both, computes energy/forces/stress and writes a JSON result
      file (see Modules/ASEClusterRunner.py).

    Usage:
        calc = ASEFileCalculator(EMT())
        cluster.binary = calc.command   # the runner command with PREFIX placeholder
        cluster.mpi_cmd = ""            # the runner is a serial python process
        cluster.compute_ensemble(ensemble, calc)

    Requirements on the (remote) cluster: python + ase (+ the actual
    calculator dependencies, and a python able to unpickle the calculator —
    same ASE version recommended). For sscha.LocalCluster this is
    automatically satisfied by the local environment.
    """

    input_extension = "_input.json"
    output_extension = ".json"
    calc_pickle_name = "calculator.pkl"
    extra_input_files = [calc_pickle_name]

    def __init__(self, ase_calc, python_exe=None):
        """
        Parameters
        ----------
            ase_calc : ase.calculators.calculator.Calculator
                The ASE calculator to execute on the cluster. Must be picklable.
            python_exe : str, optional
                The python interpreter used ON THE CLUSTER to run the runner.
                Defaults to the local interpreter (correct for LocalCluster).
        """
        super().__init__()

        # Fail early with a clear error if the calculator cannot be pickled
        pickle.dumps(ase_calc)

        self.ase_calc = ase_calc
        self.python_exe = python_exe or sys.executable
        self.command = self.get_execution_command()

    def get_execution_command(self):
        return ("{exe} -m sscha.ASEClusterRunner "
                "--calc {pkl} "
                "--input PREFIX{in_ext} "
                "--output PREFIX{out_ext}"
                ).format(exe=self.python_exe, pkl=self.calc_pickle_name,
                         in_ext=self.input_extension, out_ext=self.output_extension)

    def copy(self):
        """Return an identical instance, without inheriting calculation info."""
        return ASEFileCalculator(self.ase_calc, self.python_exe)

    def set_directory(self, directory):
        """Set the working directory and pickle the calculator there (once)."""
        CC.calculators.FileIOCalculator.set_directory(self, directory)
        pkl_path = os.path.join(directory, self.calc_pickle_name)
        if not os.path.exists(pkl_path):
            with open(pkl_path, "wb") as handle:
                pickle.dump(self.ase_calc, handle)

    def write_input(self, structure):
        """Serialize the structure as JSON in {directory}/{label}{input_extension}."""
        import ase.io.jsonio
        CC.calculators.FileIOCalculator.write_input(self, structure)

        atoms = structure.get_ase_atoms()
        fname = os.path.join(self.directory, self.label + self.input_extension)
        with open(fname, "w") as handle:
            handle.write(ase.io.jsonio.encode(atoms))

    def read_results(self):
        """Read {directory}/{label}{output_extension} into .results/.structure."""
        import ase.io.jsonio

        fname = os.path.join(self.directory, self.label + self.output_extension)
        with open(fname, "r") as handle:
            payload = ase.io.jsonio.decode(handle.read())

        self.results = {
            "energy": payload["energy"],   # eV
            "forces": payload["forces"],   # eV / Angstrom
            "stress": payload["stress"],   # ASE Voigt (xx,yy,zz,yz,xz,xy), -eV/Angstrom^3
        }
        self.structure = CC.Structure.Structure()
        self.structure.generate_from_ase_atoms(payload["atoms"])


# ---------------------------------------------------------------------- #
#  Direct (in-process) family (DirectCluster)                             #
# ---------------------------------------------------------------------- #

class DirectCalculator(ClusterCalculator):
    """Mid-level interface for in-process calculators (DirectCluster)."""

    def compute(self, structure):
        """Compute one structure and return the raw results dict.

        Parameters
        ----------
            structure : cellconstructor.Structure.Structure

        Returns
        -------
            dict with "energy" [eV], "forces" [eV/Ang], optionally
            "stress" [ASE Voigt, -eV/Ang^3] and "structure" [CC.Structure].
        """
        raise NotImplementedError


class ASEDirectCalculator(DirectCalculator):
    """Wrap ANY ASE calculator for in-process use with DirectCluster.

    No files are written: the ASE calculator is called directly on the
    ASE image of each structure.
    """

    def __init__(self, ase_calc):
        self.ase_calc = ase_calc

    def copy(self):
        """Independent per-thread copy (ASE calculators are stateful)."""
        return ASEDirectCalculator(copy.deepcopy(self.ase_calc))

    def compute(self, structure):
        atoms = structure.get_ase_atoms()
        atoms.calc = self.ase_calc
        results = {
            "energy": float(atoms.get_potential_energy()),  # eV
            "forces": atoms.get_forces(),                   # eV / Ang
            "structure": structure.copy(),                  # trivial consistency check
        }
        try:
            results["stress"] = atoms.get_stress()          # ASE Voigt, -eV/Ang^3
        except Exception:
            pass  # calculator without stress: call compute_ensemble with get_stress=False
        return results
