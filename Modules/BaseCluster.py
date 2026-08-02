# -*- coding: utf-8 -*-
"""Generic interface for computing SSCHA ensembles 'on a cluster'.

A 'cluster' here is any execution backend capable of computing energies,
forces and (optionally) stresses for a list of structures: a remote HPC
with a job scheduler (sscha.Cluster.Cluster), the same machinery mocked
locally (sscha.LocalCluster.LocalCluster), or a direct in-process
calculator (DirectCluster).

The ensemble driver (compute_ensemble_batch), including the on-the-fly
(FLARE) learning loop, is implemented ONCE in this module. Backends only
implement `compute_jobarray`, which obtains the raw results for a list of
ensemble indices.
"""
from __future__ import print_function

import sys
import os
import threading
import difflib

import numpy as np

# SETUP THE CODATA 2006, To match the QE definition of Rydberg (as in Cluster.py)
try:
    import ase.units
    units = ase.units.create_units("2006")
except Exception:
    units = {"Ry": 13.605698066, "Bohr": 1 / 1.889725989}


class BaseCluster(object):
    """Abstract execution backend for ensemble calculations.

    Subclasses must implement :func:`compute_jobarray` and declare the
    calculator interface they consume via `_required_calculator_interface`.
    """

    # Attributes a calculator must expose to be used with this cluster.
    _required_calculator_interface = ("copy",)

    def __init__(self, batch_size=1000, job_number=1, max_recalc=10):
        """
        Parameters
        ----------
            batch_size : int
                Maximum number of job arrays computed in parallel per cycle.
                With on-the-fly learning active, one cycle is the learning
                batch: the GP is updated/retrained (and the remaining
                structures re-predicted) after each cycle, so the effective
                OTF batch size is `batch_size * job_number` structures.
            job_number : int
                Number of structures grouped in a single job array.
            max_recalc : int
                Maximum number of resubmission cycles for failed jobs.
        """
        self.batch_size = batch_size
        self.job_number = job_number
        self.max_recalc = max_recalc
        self.lock = None  # threading.Lock(), created per compute_ensemble_batch

        # NOTE: subclasses set their own attributes and then call
        # self._lock_attributes() at the end of their __init__.
        # Setting attributes BEFORE super().__init__() also works (they are
        # picked up into __total_attributes__), as OpticalQECluster does.

    # ------------------ attribute locking machinery ------------------ #
    # (moved verbatim from sscha.Cluster.Cluster)

    def _lock_attributes(self):
        """Forbid setting attributes not defined up to this point."""
        self.__total_attributes__ = [item for item in self.__dict__.keys()]
        self.fixed_attributes = True  # This must be the last attribute to be setted

    def __setattr__(self, name, value):
        if "fixed_attributes" in self.__dict__:
            if name in self.__total_attributes__:
                super(BaseCluster, self).__setattr__(name, value)
            elif self.fixed_attributes:
                similar_objects = str(difflib.get_close_matches(name, self.__total_attributes__))
                ERROR_MSG = """
        Error, the attribute '{}' is not a member of '{}'.
        Suggested similar attributes: {} ?
        """.format(name, type(self).__name__, similar_objects)
                raise AttributeError(ERROR_MSG)

            if name.endswith("_name"):
                key = "use_{}".format(name.split("_")[0])
                self.__dict__[key] = True
        else:
            super(BaseCluster, self).__setattr__(name, value)

    def __getstate__(self):
        """Return the picklable state (the thread lock cannot be pickled)."""
        state = self.__dict__.copy()
        state["lock"] = None
        return state

    def __setstate__(self, state):
        state["lock"] = None
        self.__dict__.update(state)

    # ------------------ template-method hooks ------------------ #

    def _check_calculator_interface(self, calc):
        """Fail early with a clear error if calc cannot run on this cluster."""
        missing = [m for m in self._required_calculator_interface
                   if not hasattr(calc, m)]
        if missing:
            raise TypeError(
                "Error, the calculator {} cannot be used with {}.\n"
                "Missing methods/attributes: {}".format(
                    type(calc).__name__, type(self).__name__, missing))

    def _pre_compute_hook(self, ensemble, calc):
        """Called once before the first submission cycle (default: nothing)."""

    def compute_jobarray(self, ensemble, calc, jobs_id):
        """Compute one job array and return the raw results.

        MUST be implemented by subclasses.

        Parameters
        ----------
            ensemble : sscha.Ensemble.Ensemble
                The ensemble being computed.
            calc : a calculator of the interface required by this cluster
                A private copy for this job array (created with calc.copy()).
            jobs_id : list of int
                The ensemble indices to compute.

        Returns
        -------
            list, aligned with jobs_id, of raw result dicts or None (failure).
            Each dict must provide "energy" [eV], "forces" [eV/Angstrom],
            optionally "stress" (ASE Voigt xx,yy,zz,yz,xz,xy, -eV/Angstrom^3)
            and "structure" (CC.Structure, for the consistency check);
            any extra key is stored in ensemble.all_properties.
        """
        raise NotImplementedError("compute_jobarray must be implemented by subclasses")

    # ------------------ the generic ensemble driver ------------------ #

    def compute_ensemble(self, ensemble, calc, get_stress=True, timeout=None):
        """Run the whole ensemble on this cluster (see compute_ensemble_batch)."""
        self.compute_ensemble_batch(ensemble, calc, get_stress, timeout)

    def compute_ensemble_batch(self, ensemble, calc, get_stress=True, timeout=None):
        """
        RUN THE ENSEMBLE WITH BATCH SUBMISSION (generic driver)
        =======================================================

        If the ensemble has an active on-the-fly ML model
        (``ensemble.gp_model is not None``, set via ``ensemble.set_otf``),
        each cycle of ab-initio computations is followed by a GP
        update/retrain, and the remaining structures are re-checked against
        the model: those predicted with acceptable uncertainty are filled
        with the ML prediction and never computed. One learning cycle
        computes up to ``batch_size * job_number`` structures.
        """
        self._check_calculator_interface(calc)

        # Track the remaining configurations
        success = [False] * ensemble.N

        # Setup if the ensemble has the stress
        ensemble.has_stress = get_stress

        self._pre_compute_hook(ensemble, calc)

        # Get the expected number of batch
        num_batch_offset = int(ensemble.N / self.batch_size)

        # ==================== OTF SETUP (FLARE) ====================
        use_otf = getattr(ensemble, "gp_model", None) is not None
        dft_counts = 0
        if use_otf:
            number_of_atoms = ensemble.structures[0].get_ase_atoms().get_global_number_of_atoms()
            ensemble._otf_setup_defaults(number_of_atoms)
            remaining = list(range(ensemble.N))
            ensemble._otf_predict(ensemble.structures, remaining)
            for i in range(ensemble.N):
                if ensemble.force_computed[i]:
                    success[i] = True

        # Run until some work has not finished
        recalc = 0
        self.lock = threading.Lock()
        while np.sum(np.array(success, dtype=int) - 1) != 0:
            threads = []
            cycle_results = {}

            print("[CYCLE] SUCCESS: ", success)
            print("[CYCLE] STOPPING CONDITION:", np.sum(np.array(success, dtype=int) - 1))

            # Get the remaining jobs
            false_mask = np.array(success) == False
            false_id = np.arange(ensemble.N)[false_mask]

            count = 0
            # Submit in parallel
            jobs = [false_id[i:i + self.job_number] for i in range(0, len(false_id), self.job_number)]
            # Create a local copy of the calculator for each thread, to avoid conflicting modifications
            calculators = [calc.copy() for i in range(0, len(jobs))]

            for k_th, job in enumerate(jobs):
                # Submit only the batch size
                if count >= self.batch_size:
                    break
                t = threading.Thread(target=self._compute_jobarray_thread,
                                     args=(ensemble, calculators[k_th], job,
                                           get_stress, cycle_results, success))
                t.start()
                threads.append(t)
                count += 1

            # Wait until all the job have finished
            for t in threads:
                t.join(timeout)

            # ============ OTF UPDATE / TRAIN / PREDICT ============
            if use_otf and cycle_results:
                # Main thread only, deterministic order (G7)
                dft_counts += len(cycle_results)
                for num in sorted(cycle_results):
                    res = cycle_results[num]
                    dft_stress = np.array(res["stress"], dtype=float) \
                        if (get_stress and "stress" in res) else None
                    ensemble._otf_update_from_structure(
                        ensemble.structures[num],
                        dft_energy=res["energy"],     # eV
                        dft_frcs=res["forces"],       # eV/Ang
                        dft_stress=dft_stress,        # ASE Voigt, -eV/Ang^3
                    )
                ensemble._otf_maybe_train_and_write(dft_counts)
                # Re-check the remaining structures against the model
                remaining = [int(i) for i in np.arange(ensemble.N)[np.array(success) == False]]
                ensemble._otf_predict(ensemble.structures, remaining)
                for i in range(ensemble.N):
                    if ensemble.force_computed[i]:
                        success[i] = True

            print("[CYCLE] [END] SUCCESS: ", success)
            print("[CYCLE] [END] STOPPING CONDITION:", np.sum(np.array(success, dtype=int) - 1))

            recalc += 1
            if recalc > num_batch_offset + self.max_recalc:
                print("Expected batch ordinary resubmissions:", num_batch_offset)
                raise ValueError("Error, resubmissions exceeded the maximum number of %d" % self.max_recalc)

        if use_otf:
            ensemble._clean_runs(dft_counts)

        print("CALCULATION ENDED: all properties: {}".format(ensemble.all_properties))

    def _compute_jobarray_thread(self, ensemble, calc, jobs_id,
                                 get_stress, cycle_results, success):
        """Thread worker: obtain raw results, then ingest them (locked)."""
        raw_results = self.compute_jobarray(ensemble, calc, jobs_id)

        # Thread safe operation
        self.lock.acquire()
        try:
            print("[THREAD {}] submitted calculations: {}".format(
                threading.get_native_id(), list(jobs_id)))
            for pos, res in enumerate(raw_results):
                num = int(jobs_id[pos])
                print("[THREAD {}] ADDING RESULT {} = {}".format(
                    threading.get_native_id(), num, res))
                ok = self._ingest_result(ensemble, res, num, get_stress)
                success[num] = ok
                if ok:
                    cycle_results[num] = res   # keep raw eV results for the OTF update (G6)
        finally:
            self.lock.release()

    def _ingest_result(self, ensemble, res, num, get_stress):
        """Validate one raw result and write it into the ensemble arrays.

        Returns True if the result was complete and stored, False otherwise
        (the job will be resubmitted, up to max_recalc).
        """
        if res is None:
            return False

        # Check if the run was good
        check_e = "energy" in res
        check_f = "forces" in res
        check_s = "stress" in res

        # Check the structure
        if "structure" in res:
            error_struct = np.linalg.norm(ensemble.structures[num].coords.ravel()
                                          - res["structure"].coords.ravel())
            if error_struct > 1e-2:
                print("ERROR IDENTIFYING STRUCTURE!")
                MSG = """
                    Error in thread {}.
                    Displacement between the expected structure {}
                    and the one readed from the calculator
                    is of {} A.
                """.format(threading.get_native_id(), num, error_struct)
                print(MSG)
                ensemble.structures[num].save_scf(
                    't_{}_error_struct_generated_{}.scf'.format(threading.get_native_id(), num))
                res["structure"].save_scf(
                    't_{}_error_struct_readed_{}.scf'.format(threading.get_native_id(), num))
                return False
        else:
            print("[WARNING] no check on the structure.")

        is_success = check_e and check_f
        if get_stress:
            is_success = is_success and check_s

        if not is_success:
            return False

        res_only_extra = {x: res[x] for x in res if x not in ["energy", "forces", "stress", "structure"]}
        ensemble.all_properties[num].update(res_only_extra)
        ensemble.energies[num] = res["energy"] / units["Ry"]
        ensemble.forces[num, :, :] = res["forces"] / units["Ry"]
        ensemble.force_computed[num] = True

        if get_stress:
            stress = np.zeros((3, 3), dtype=np.float64)
            stress[0, 0] = res["stress"][0]
            stress[1, 1] = res["stress"][1]
            stress[2, 2] = res["stress"][2]
            stress[1, 2] = res["stress"][3]
            stress[2, 1] = res["stress"][3]
            stress[0, 2] = res["stress"][4]
            stress[2, 0] = res["stress"][4]
            stress[0, 1] = res["stress"][5]
            stress[1, 0] = res["stress"][5]
            # Remember, ase has a very strange definition of the stress
            ensemble.stresses[num, :, :] = -stress * units["Bohr"]**3 / units["Ry"]
            ensemble.stress_computed[num] = True
        return True


class DirectCluster(BaseCluster):
    """
    DIRECT (IN-PROCESS) CLUSTER
    ===========================

    A 'cluster' that computes the ensemble directly in the current process,
    without writing any file and without any job scheduler. It consumes
    `DirectCalculator` objects (e.g. ASEDirectCalculator wrapping any ASE
    calculator).

    The learning-cycle size is still `batch_size * job_number`; threads are
    used across job arrays exactly like in the scheduler-based clusters
    (each thread owns a private calculator copy). Note that pure-Python ASE
    calculators are bound by the GIL; the parallelism is still useful for
    calculators that release it (NumPy-heavy or subprocess-based codes), and
    the interface stays identical to the other clusters.
    """

    _required_calculator_interface = ("copy", "compute")

    def __init__(self, batch_size=1000, job_number=1, max_recalc=10):
        super().__init__(batch_size=batch_size, job_number=job_number,
                         max_recalc=max_recalc)
        self._lock_attributes()

    def compute_jobarray(self, ensemble, calc, jobs_id):
        """Compute each structure in-process via calc.compute(structure)."""
        results = []
        for num in jobs_id:
            try:
                results.append(calc.compute(ensemble.structures[int(num)]))
            except Exception as exc:
                sys.stderr.write("JOB {} resulted in error:\n{}\n".format(num, exc))
                sys.stderr.flush()
                results.append(None)
        return results
