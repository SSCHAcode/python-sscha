# -*- coding: utf-8 -*-
import sys, os
import threading
import numpy as np

from ase.units import Rydberg, Bohr
import ase, ase.io

"""
This is an untility script that is able to manage the submission into
a cluster of an ensemble
"""


class Cluster:
    
    hostname=""
    pwd=""
    sshcmd="ssh"
    workdir=r""
    submit_command="sbatch --wait"
    submit_name="SBATCH"
    terminal="#!/bin/bash"
    v_nodes="-N"
    use_nodes = True
    v_cpu="-n"
    use_cpu = True
    v_account="-A"
    use_account = True
    v_time="--time="
    use_time = True
    v_memory="--mem="
    use_memory = False
    v_partition="--partition="
    use_partition= False
    load_modules=""
    n_nodes = 1
    n_cpu = 1
    n_pool = 1
    max_recalc = 10
    
    
    def __init__(self, hostname, pwd=None, extra_options="", workdir = "",
                 account_name = "", partition_name = "", binary="pw.x -npool $NPOOL -i PREFIX.pwi > PREFIX.pwo",
                 mpi_cmd=r"srun --mpi=pmi2 -n $NPROC"):
        """
        SETUP THE CLUSTER
        =================
        It is strongly suggested to use public/private keys.
        However, sometimes for some reasons it is not possible to do so.
        Then you must install sshpass, and provide the password 
        for establishing the connection.
        
        Parameters
        ----------
            hostname:
                The name of the host toward you want to connect.
                E.G. 
                pippo@login.cineca.marconi.it
            pwd:
                The password for the connection. 
                If you use a private key, do not pass this argument.
                In the other case you must have sshpass package to allow the
                password communication.
            extra_options:
                The extra options to be passed to the ssh command.
                For example, connect to port 24 use 'extra_options="-p 24"'
            workdir:
                The workinig directory in the cluster. This is the directory
                in which the job will be runned.
            account_name:
                The name of the account for the job submission
            partition_name:
                The partition name for the job submission
        """
        
        self.hostname = hostname
        if not pwd is None:
            self.pwd = pwd
            res = os.system("sshpass")
            if res != 0:
                raise ValueError("Error, sshpass command not found, required to connect through password to server.")
            
            self.sshcmd = "sshpass -p '" + pwd + "' ssh " + extra_options
            self.scpcmd = "sshpass -p '" + pwd + "' scp " + extra_options.replace("-p", "-P")
        else:
            self.sshcmd = "ssh " + extra_options
            self.scpcmd = "scp " + extra_options.replace("-p", "-P")

    
        self.account_name = account_name
        self.partition_name = partition_name
        if partition_name:
            self.use_partition = True
            
        self.binary = binary
        self.mpi_cmd = mpi_cmd
            
            
    def run_atoms(self, ase_calc, ase_atoms, label="ESP", 
                  in_extension = ".pwi", out_extension=".pwo",
                  n_nodes = 1, n_cpu=1, npool=1, ram=10000, time="00:02:00"):
        """
        RUN ATOMS ON THE CLUSTER
        ========================
        
        This function runs the given atoms in the cluster, using the ase_calculator.
        Note: the ase_calc must be a FileIOCalculator. 
        For now it works with quantum espresso.
        """
        
        # Setup the calculator with the atomic structure
        ase_atoms.set_calculator(ase_calc)
        
        # Write the input file
        ase.io.write("%s%s" % (label, in_extension),ase_atoms,**ase_calc.parameters)
        
        # prepare the submission script
        submission = self.terminal + "\n"
        
        # Add the submission options
        if self.use_nodes:
            submission += "#%s %s%d\n" % (self.submit_name, self.v_nodes, n_nodes)
        if self.use_cpu:
            submission += "#%s %s%d\n" % (self.submit_name, self.v_cpu, n_cpu)
        if self.use_time:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_time, time)
        if self.use_account:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_account, self.account_name)
        if self.use_memory:
            submission += "#%s %s%d\n" % (self.submit_name, self.v_memory, ram)
        if self.use_partition:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_partition, self.partition_name)
        
        # Add the loading of the modules
        submission += self.load_modules + "\n"
        
        # Go to the working directory
        submission += "cd " + self.workdir + "\n"
        
        # Get the real calculation command
        mpicmd = self.mpi_cmd.replace("$NPROC", str(n_cpu))
        binary = self.binary.replace("$NPOOL", str(npool)).replace("PREFIX", label)
        
        submission += mpicmd + " " + binary + "\n"
        
        
        # Copy the input files into the target directory
        f = file("%s.sh" % label, "w")
        f.write(submission)
        f.close()
        cmd = self.scpcmd + " %s.sh %s:%s" % (label, self.hostname, self.workdir)
        cp_res = os.system(cmd)
        if cp_res != 0:
            print "Error while executing:", cmd
            print "Return code:", cp_res
            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
            
        cmd = self.scpcmd + " %s%s %s:%s" % (label, in_extension, self.hostname, self.workdir)
        cp_res = os.system(cmd)
        if cp_res != 0:
            print "Error while executing:", cmd
            print "Return code:", cp_res
            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
        
        # Run the simulation
        cmd = self.sshcmd + " %s '%s %s/%s.sh'" % (self.hostname, self.submit_command, self.workdir, label)
        cp_res = os.system(cmd)
        if cp_res != 0:
            print "Error while executing:", cmd
            print "Return code:", cp_res
            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
            
        # Get the response
        cmd = self.scpcmd + " %s:%s/%s%s ./" % (self.hostname, self.workdir, label, out_extension) 
        cp_res = os.system(cmd)
        if cp_res != 0:
            print "Error while executing:", cmd
            print "Return code:", cp_res
            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
            
        # Get the results
        ase_calc.set_label(label)
        ase_calc.read_results()
        
        return ase_calc.results
    
    def compute_ensemble(self, ensemble, ase_calc, get_stress = True, timeout=None):
        """
        RUN THE WHOLE ENSEMBLE ON THE CLUSTER
        =====================================
        
        Parameters
        ----------
            ensemble :
                The ensemble to be runned.
        """
        # Track the remaining configurations
        success = [False] * ensemble.N
        
        # Prepare the function for the simultaneous submission
        def compute_single(num, calc):
            atm = ensemble.structures[num].get_ase_atoms()
            res = self.run_atoms(calc, atm, "ESP_" + str(num), 
                                 n_nodes = self.n_nodes, 
                                 n_cpu=self.n_cpu,
                                 npool = self.n_pool)
            if res:
                ensemble.energies[num] = res["energy"] / Rydberg
                ensemble.forces[num, :, :] = res["forces"] / Rydberg
                if get_stress:
                    stress = np.zeros((3,3), dtype = np.float64)
                    stress[0,0] = res["stress"][0]
                    stress[1,1] = res["stress"][1]
                    stress[2,2] = res["stress"][2]
                    stress[1,2] = res["stress"][3]
                    stress[2,1] = res["stress"][3]
                    stress[0,2] = res["stress"][4]
                    stress[2,0] = res["stress"][4]
                    stress[0,1] = res["stress"][5]
                    stress[1,0] = res["stress"][5]
                    ensemble.stresses[num, :, :] = stress * Bohr**3 / Rydberg
                success[num] = True
        
        # Run until some work has not finished
        recalc = 0
        while np.sum(np.array(success, dtype = int) - 1) != 0:
            threads = []
            
            # Get the remaining jobs
            false_mask = np.array(success) == False
            false_id = np.arange(ensemble.N)[false_mask]
            
            # Submit in parallel
            for i in false_id:
                t = threading.Thread(target = compute_single, args=(i, ase_calc, ))
                t.start()
                threads.append(t)
            
            # Wait until all the job have finished
            for t in threads:
                t.join(timeout)
            
            recalc += 1
            if recalc > self.max_recalc:
                raise ValueError("Error, recalculations exceeded the maximum number of %d" % self.max_recalc)
                break
            
            
        