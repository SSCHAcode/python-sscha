# -*- coding: utf-8 -*-
import sys, os
import subprocess
import threading

import numpy as np

from ase.units import Rydberg, Bohr
import ase, ase.io

import cellconstructor as CC
import cellconstructor.Methods

# SETUP THE CODATA 2006, To match the QE definition of Rydberg
units = ase.units.create_units("2006")

"""
This is an untility script that is able to manage the submission into
a cluster of an ensemble
"""
__CLUSTER_NAMELIST = "cluster"
__CLUSTER_HOST = "hostname"
__CLUSTER_PWD = "pwd"


class Cluster:
    
    # The host name for the connection
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
    
    # This is a set of lines to be runned before the calculation
    # It can be used to load the needed modules in the cluster
    load_modules=""
    
    # Setup the number of cpu, nodes and pools for the calculation
    n_nodes = 1
    n_cpu = 1
    n_pool = 1
    
    # This is the default label, change it if you are going to submit
    # two different calculations in the same working directory
    label = "ESP_"
    
    # This is the maximum number of resubmissions after the expected one
    # from the batch size. It can be use to resubmit the failed jobs.
    max_recalc = 10
    
    # This is the time string. Faster job will be the less in the queue,
    time="00:02:00" # 2minutes
    
    # The ram required for the calculation
    ram="10000Mb" # 10Gb
    
    # Still unused
    prefix_name = "prefix" # Variable in the calculator for differentiating the calculations
    
    # This directory is used to work with clusters.
    local_workdir = "cluster_work/"
    
    # The batch size is the maximum number of job to be submitted together.
    # The new jobs will be submitted only after a batch is compleated
    # Useful if the cluster has a limit for the maximum number of jobs allowed.
    batch_size = 1000
    
    
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
            res = os.system("sshpass > /dev/null")
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
        
        
            
    def CheckCommunication(self):
        """
        CHECK IF THE SERVER IS AVAILABLE
        ================================
        
        This function return true if the server respond correctly, 
        false otherwise.
        """
        
        cmd = self.sshcmd + " %s 'echo ciao'" % self.hostname
        print cmd
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
        output, err = p.communicate()
        output = output.strip()
        status = p.wait()
        if not err is None:
            sys.stderr.write(err)
        if status != 0:
            sys.stderr.write("Error, cmd: " + cmd + "\n")
            sys.stderr.write("Exit status:" + str(status))
            return False
    
        if output != "ciao":
            sys.stderr.write("Code exited correctly but echo did not answer properly\n")
            sys.stderr.write("cmd: " + cmd + "\n")
            sys.stderr.write("out: " + output + "\n")
            sys.stderr.write("expected " + "ciao" + "\n")
            return False
        
        return True
            
    def run_atoms(self, ase_calc, ase_atoms, label="ESP", 
                  in_extension = ".pwi", out_extension=".pwo",
                  n_nodes = 1, n_cpu=1, npool=1):
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
        ase.io.write("%s/%s%s" % (self.local_workdir, label, in_extension),ase_atoms,**ase_calc.parameters)
        
        # prepare the submission script
        submission = self.terminal + "\n"
        
        # Add the submission options
        if self.use_nodes:
            submission += "#%s %s%d\n" % (self.submit_name, self.v_nodes, n_nodes)
        if self.use_cpu:
            submission += "#%s %s%d\n" % (self.submit_name, self.v_cpu, n_cpu)
        if self.use_time:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_time, self.time)
        if self.use_account:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_account, self.account_name)
        if self.use_memory:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_memory, self.ram)
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
        
        # First of all clean eventually input/output file of this very same calculation
        cmd = self.sshcmd + " %s 'rm -f %s/%s%s %s/%s%s'" % (self.hostname, 
                                                             self.workdir, label, in_extension,
                                                             self.workdir, label, out_extension)
        cp_res = os.system(cmd)
        if cp_res != 0:
            print "Error while executing:", cmd
            print "Return code:", cp_res
            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
        
        # Copy the input files into the target directory
        f = file("%s/%s.sh" % (self.local_workdir, label), "w")
        f.write(submission)
        f.close()
        cmd = self.scpcmd + " %s/%s.sh %s:%s" % (self.local_workdir, label, self.hostname, self.workdir)
        cp_res = os.system(cmd)
        if cp_res != 0:
            print "Error while executing:", cmd
            print "Return code:", cp_res
            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
            
        cmd = self.scpcmd + " %s/%s%s %s:%s" % (self.local_workdir, label, in_extension, self.hostname, self.workdir)
        cp_res = os.system(cmd)
        if cp_res != 0:
            print "Error while executing:", cmd
            print "Return code:", cp_res
            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
            return
        
        # Run the simulation
        cmd = self.sshcmd + " %s '%s %s/%s.sh'" % (self.hostname, self.submit_command, self.workdir, label)
        cp_res = os.system(cmd)
        if cp_res != 0:
            print "Error while executing:", cmd
            print "Return code:", cp_res
            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
            
        # Get the response
        cmd = self.scpcmd + " %s:%s/%s%s %s/" % (self.hostname, self.workdir, label, out_extension, 
                                                 self.local_workdir) 
        cp_res = os.system(cmd)
        if cp_res != 0:
            print "Error while executing:", cmd
            print "Return code:", cp_res
            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
            return
            
        # Get the results
        ase_calc.set_label("%s/%s" % (self.local_workdir, label))
        ase_calc.read_results()
        
        return ase_calc.results
    
    def setup_from_namelist(self, namelist):
        """
        SETUP THE CLUSTER WITH AN INPUTFILE
        ===================================
        
        This method setup the cluster using a custom input file.
        The inputfile must have the same shape of QuantumESPRESSO ones.
        The information about the cluster must be located in a namespace called
        as __CLUSTER_NAMELIST
        
        Parameters
        ----------
            namelist: 
                The parsed namelist dictionary.
        """
        
        
        # Check if the cluster namelist is present
        if not __CLUSTER_NAMELIST in namelist.keys():
            raise ValueError("Error, the parsed namelist must contain %s" % __CLUSTER_NAMELIST)
            
        c_info = namelist[__CLUSTER_NAMELIST]
        keys = c_info.keys()
        # Setup all the info
        if __CLUSTER_HOST in keys:
            self.hostname = c_info[__CLUSTER_HOST]
        
        if __CLUSTER_PWD in keys:
            self.pwd = c_info[__CLUSTER_PWD]
            # Check if sshpass is present
            res = os.system("sshpass > /dev/null") 
            if res != 0:
                raise ValueError("Error, sshpass command not found, required to connect through password to server.")
            
            self.sshcmd = "sshpass -p '" + pwd + "' ssh " + extra_options
            self.scpcmd = "sshpass -p '" + pwd + "' scp " + extra_options.replace("-p", "-P")
        
    
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
        
        # Setup if the ensemble has the stress
        ensemble.has_stress = get_stress
        
        # Check if the working directory exists
        if not os.path.isdir(self.local_workdir):
            os.makedirs(self.local_workdir)
            
        # Prepare the function for the simultaneous submission
        def compute_single(num, calc):
            atm = ensemble.structures[num].get_ase_atoms()
            res = self.run_atoms(calc, atm, self.label + str(num), 
                                 n_nodes = self.n_nodes, 
                                 n_cpu=self.n_cpu,
                                 npool = self.n_pool)
            if res:
                ensemble.energies[num] = res["energy"] / units["Ry"]
                ensemble.forces[num, :, :] = res["forces"] / units["Ry"]
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
                    # Remember, ase has a very strange definition of the stress
                    ensemble.stresses[num, :, :] = -stress * units["Bohr"]**3 / units["Ry"]
                success[num] = True
                
        # Get the expected number of batch
        num_batch_offset = int(ensemble.N / self.batch_size)
        
        # Run until some work has not finished
        recalc = 0
        while np.sum(np.array(success, dtype = int) - 1) != 0:
            threads = []
            
            # Get the remaining jobs
            false_mask = np.array(success) == False
            false_id = np.arange(ensemble.N)[false_mask]
            
            count = 0
            # Submit in parallel
            for i in false_id:
                # Submit only the batch size
                if count > self.batch_size:
                    break
                t = threading.Thread(target = compute_single, args=(i, ase_calc, ))
                t.start()
                threads.append(t)
                count += 1
                
            
            # Wait until all the job have finished
            for t in threads:
                t.join(timeout)
            
            recalc += 1
            if recalc > num_batch_offset + self.max_recalc:
                print "Expected batch ordinary resubmissions:", num_batch_offset
                raise ValueError("Error, resubmissions exceeded the maximum number of %d" % self.max_recalc)
                break
            
            
        