# -*- coding: utf-8 -*-
from __future__ import print_function
import sys, os
import subprocess
import threading
import copy

__DIFFLIB__ = False
try:
    __DIFFLIB__ = True
    import difflib
except:
    pass

import numpy as np

from ase.units import Rydberg, Bohr
import ase, ase.io

import cellconstructor as CC
import cellconstructor.Methods

# SETUP THE CODATA 2006, To match the QE definition of Rydberg
try:
    units = ase.units.create_units("2006")
except:
    units = {"Ry": 13.605698066, "Bohr": 1/1.889725989}

"""
This is an untility script that is able to manage the submission into
a cluster of an ensemble
"""
__CLUSTER_NAMELIST__ = "cluster"
__CLUSTER_TEMPLATE__ = "template"
__TEMPLATE_ENV__ = "SSCHA_CLUSTERS_DIR"
__CLUSTER_HOST__ = "hostname"
__CLUSTER_PWD__ = "pwd"
__CLUSTER_ACCOUNT__ = "account"
__CLUSTER_BINARY__ = "binary_path"
__CLUSTER_MPICMD__ = "mpicmd"
__CLUSTER_ATTEMPTS__ = "reconnect_attempts"
__CLUSTER_PORT__ = "port"


__CLUSTER_TERMINAL__ = "shell"
__CLUSTER_SUBCMD__ = "submit_cmd"
__CLUSTER_SUBNAME__ = "queue_directive"
__CLUSTER_VNODES__ = "v_nodes"
__CLUSTER_NNODES__ = "n_nodes"
__CLUSTER_UNODES__ = "use_nodes"
__CLUSTER_VCPU__ = "v_cpu" 
__CLUSTER_NCPU__ = "n_cpu"
__CLUSTER_UCPU__ = "use_cpu"
__CLUSTER_VTIME__ = "v_time"
__CLUSTER_NTIME__ = "n_time"
__CLUSTER_NPOOLS__ = "n_pools"
__CLUSTER_UTIME__ = "use_time"
__CLUSTER_VMEM__ = "v_memory"
__CLUSTER_NMEM__ = "max_ram"
__CLUSTER_UMEM__ = "use_memory"
__CLUSTER_VPART__ = "v_partition"
__CLUSTER_NPART__ = "partition_name"
__CLUSTER_UPART__ = "use_partition"
__CLUSTER_INITSCRIPT__ = "init_script"
__CLUSTER_MAXRECALC__ = "max_recalc"
__CLUSTER_BATCHSIZE__ = "batch_size"
__CLUSTER_LOCALWD__ = "local_workdir"
__CLUSTER_VACCOUNT__ = "v_account"
__CLUSTER_UACCOUNT__ = "use_account"
__CLUSTER_SSHCMD__ = "sshcmd"
__CLUSTER_SCPCMD__ = "scpcmd"
__CLUSTER_TIMEOUT__ = "timeout"
__CLUSTER_JOBNUMBER__ = "job_numbers"
__CLUSTER_NPARALLEL__ = "n_together"


__CLUSTER_WORKDIR__ = "workdir"


# List of the requested keywords
__CLUSTER_RKW__ = [__CLUSTER_WORKDIR__, __CLUSTER_HOST__,
                   __CLUSTER_ACCOUNT__, __CLUSTER_BINARY__]

# List all the possible keys
__CLUSTER_KEYS__ = [__CLUSTER_NAMELIST__, __CLUSTER_TEMPLATE__, __CLUSTER_HOST__, __CLUSTER_PWD__,
                    __CLUSTER_ACCOUNT__, __CLUSTER_BINARY__, __CLUSTER_MPICMD__, __CLUSTER_TERMINAL__,
                    __CLUSTER_SUBCMD__, __CLUSTER_SUBNAME__, __CLUSTER_VNODES__, __CLUSTER_NNODES__,
                    __CLUSTER_UNODES__, __CLUSTER_VCPU__, __CLUSTER_NCPU__, __CLUSTER_UCPU__,
                    __CLUSTER_VTIME__, __CLUSTER_NTIME__, __CLUSTER_UTIME__, __CLUSTER_VMEM__,
                    __CLUSTER_NMEM__, __CLUSTER_UMEM__, __CLUSTER_VPART__, __CLUSTER_NPART__,
                    __CLUSTER_UPART__, __CLUSTER_INITSCRIPT__, __CLUSTER_MAXRECALC__, __CLUSTER_BATCHSIZE__,
                    __CLUSTER_LOCALWD__, __CLUSTER_VACCOUNT__, __CLUSTER_UACCOUNT__, __CLUSTER_SSHCMD__,
                    __CLUSTER_SCPCMD__, __CLUSTER_WORKDIR__, __CLUSTER_TIMEOUT__, 
                    __CLUSTER_JOBNUMBER__, __CLUSTER_NPARALLEL__, __CLUSTER_NPOOLS__,
                    __CLUSTER_ATTEMPTS__, __CLUSTER_PORT__]


class Cluster(object):
    
    def __init__(self, hostname=None, pwd=None, extra_options="", workdir = "",
                 account_name = "", partition_name = "", binary="pw.x -npool NPOOL -i PREFIX.pwi > PREFIX.pwo",
                 mpi_cmd=r"srun --mpi=pmi2 -n NPROC"):
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
        self.pwd = pwd

        self.sshcmd = "ssh"
        self.sshcmd = "ssh " + extra_options
        self.scpcmd = "scp " + extra_options.replace("-p", "-P")

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

        self.workdir=r""
        self.submit_command="sbatch --wait"
        self.submit_name="SBATCH"
        self.terminal="#!/bin/bash"
        self.v_nodes="-N "
        self.use_nodes = True
        self.v_cpu="-n "
        self.use_cpu = True
        self.v_account="-A "
        self.use_account = True
        self.v_time="--time="
        self.use_time = True
        self.v_memory="--mem="
        self.use_memory = False
        self.v_partition="--partition="
        self.use_partition= False
        self.timeout = 1000
        self.use_timeout = False

        # This is the number of configurations to be computed for each jub submitted
        # This times the self.batch_size is the total amount of configurations submitted toghether
        self.job_number = 1
        self.n_together_def = 1
        self.use_multiple_submission = False


        # This is a set of lines to be runned before the calculation
        # It can be used to load the needed modules in the cluster
        self.load_modules=""
        
        # Setup the number of cpu, nodes and pools for the calculation
        self.n_nodes = 1
        self.n_cpu = 1
        self.n_pool = 1
        
        # This is the default label, change it if you are going to submit
        # two different calculations in the same working directory
        self.label = "ESP_"
        
        # This is the maximum number of resubmissions after the expected one
        # from the batch size. It can be use to resubmit the failed jobs.
        self.max_recalc = 10
        self.connection_attempts = 1
        
        # This is the time string. Faster job will be the less in the queue,
        self.time="00:02:00" # 2minutes
        
        # The ram required for the calculation
        self.ram="10000Mb" # 10Gb
        
        # The default partition in which to submit calculations
        self.partition_name = ""
        
        # Still unused
        self.prefix_name = "prefix" # Variable in the calculator for differentiating the calculations
        
        # This directory is used to work with clusters.
        self.local_workdir = "cluster_work/"
        
        # The batch size is the maximum number of job to be submitted together.
        # The new jobs will be submitted only after a batch is compleated
        # Useful if the cluster has a limit for the maximum number of jobs allowed.
        self.batch_size = 1000


        # Allow to setup additional custom extra parameters
        self.custom_params = {}

        # Fix the attributes

        # Setup the attribute control
        self.__total_attributes__ = [item for item in self.__dict__.keys()]
        self.fixed_attributes = True # This must be the last attribute to be setted


    def __setattr__(self, name, value):
        """
        This method is used to set an attribute.
        It will raise an exception if the attribute does not exists (with a suggestion of similar entries)
        """

        
        if "fixed_attributes" in self.__dict__:
            if name in self.__total_attributes__:
                super(Cluster, self).__setattr__(name, value)
            elif self.fixed_attributes:
                similar_objects = str( difflib.get_close_matches(name, self.__total_attributes__))
                ERROR_MSG = """
        Error, the attribute '{}' is not a member of '{}'.
        Suggested similar attributes: {} ?
        """.format(name, type(self).__name__,  similar_objects)

                raise AttributeError(ERROR_MSG)
        else:
            super(Cluster, self).__setattr__(name, value)


        # Setting the account or partition name will automatically result in
        # activating the corresponding flags
        if name.endswith("_name"):
            key = "use_{}".format(name.split("_")[0])
            self.__dict__[key] = True
        
        
        
    def ExecuteCMD(self, cmd, raise_error = True, return_output = False):
        """
        EXECUTE THE CMD ON THE CLUSTER
        ==============================
        
        This subroutine execute the cmd in the cluster,
        with the specified number of attempts.
        
        Parameters
        ----------
            cmd: string
                The whole command, including the ssh/scp.
            raise_error : bool, optional
                If True (default) raises an error upon failure.
            return_output : bool, optional
                If True (default False) the output of the command is 
                returned as second value.
                
        Returns
        -------
            success : bool
                If True, the command has been executed with success,
                False otherwise
            output : string
                Returned only if return_output is True
        """
        
        success = False
        output = ""
        for i in range(self.connection_attempts):
            # Subprocess opening
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True)
            output, err = p.communicate()
            if not isinstance(output, str):
                try:
                    output = str(output.decode("utf-8"))
                except:
                    pass
            output = output.strip()
            status = p.wait()
            if not err is None:
                sys.stderr.write(err)
            
            if status != 0:
                sys.stderr.write("Error with cmd: "+ cmd + "\n")
                sys.stderr.write(output + "\n")
                sys.stderr.write("EXITSTATUS: " + str(status) + "; attempt = " + str(i+1) + "\n")
            else:
                success = True
                break
            
        if raise_error and not success:
            raise IOError("Error while communicating with the cluster. More than %d attempts failed." % (i+1))
            
        
        if return_output:
            return success, output
        return success
            

    def set_timeout(self, timeout):
        """
        Set a timeout time for each single calculation.
        This is very usefull as sometimes the calculations gets stucked after some times on clusters.

        Parameters
        ----------
            timeout: int
                The timeout in seconds after which a single calculation is killed.
        """

        self.use_timeout = True
        self.timeout = timeout
        
        
            
    def CheckCommunication(self):
        """
        CHECK IF THE SERVER IS AVAILABLE
        ================================
        
        This function return true if the server respond correctly, 
        false otherwise.
        """
        
        cmd = self.sshcmd + " %s 'echo ciao'" % self.hostname
        
        status, output = self.ExecuteCMD(cmd, return_output = True)
        
#        print cmd
#        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
#        output, err = p.communicate()
#        output = output.strip()
#        status = p.wait()
#        if not err is None:
#            sys.stderr.write(err)
#        if status != 0:
#            sys.stderr.write("Error, cmd: " + cmd + "\n")
#            sys.stderr.write("Exit status:" + str(status))
#            return False
    
        if output != "ciao":
            sys.stderr.write("Code exited correctly but echo did not answer properly\n")
            sys.stderr.write("cmd: " + cmd + "\n")
            sys.stderr.write("out: " + output + "\n")
            sys.stderr.write("expected " + "ciao" + "\n")
            return False
        
        return True
    
    def batch_submission(self, list_of_structures, calc, indices, 
                         in_extension, out_extension,
                         label = "ESP", n_togheder=1):
        """
        BATCH SUBMISSION
        ================
        
        This is a different kind of submission, it exploits xargs to perform
        a parallel submission of serveral structures in one single job.
        This is very good to avoid overloading the queue system manager, or
        when a limited number of jobs per user are allowed.
        
        NOTE: the number of structure in the list must be a divisor of the
        total number of processors.
        
        Parameters
        ----------
            list_of_structures : list
                List of the structures to be computed.
            calc : ase FileIOCalculator
                The FileIOCalculator to perform the minimization
            indices : list(int)
                The indices of the configurations, this avoids interferring with
                other jobs when multiple jobs are lunched togheder.
            in_extension : string
                Extension of the input filename
            out_extension : string
                Extension of the output filename.
            label : string, optional
                The root of the input file.
            n_togheder : int, optional (DO NOT USE)
                If present, the job will lunch a new job immediately after the other 
                is ended. This is usefull to further reduce the number of submitted 
                jobs.
        
        Results
        -------
            list_of_results.
                Returns a list of results dicts, one for each structure.
        """
        N_structs  = len(list_of_structures)
        
        if n_togheder != 1:
            raise NotImplementedError("Error, n_togheder != 1 does not work!")
            
        
        # Prepare the input atoms
        app_list = ""
        new_ncpu = self.n_cpu * n_togheder
        new_mpicmd  = self.mpi_cmd.replace("NPROC", str(self.n_cpu))
        results = [None] * N_structs
        submitted = []
        for i in range(N_structs):
            # Prepare a typical label
            lbl = label + "_" + str(indices[i])
            
            atm = list_of_structures[i].get_ase_atoms()
            atm.set_calculator(calc)
            ase.io.write("%s/%s%s"% (self.local_workdir, lbl, in_extension),
                         atm,**calc.parameters)
            
            
            # Add the file in the applist
            binary = self.binary.replace("NPOOL", str(self.n_pool)).replace("PREFIX", lbl)
            
                
            # First of all clean eventually input/output file of this very same calculation
            cmd = self.sshcmd + " %s 'rm -f %s/%s%s %s/%s%s'" % (self.hostname, 
                                                                 self.workdir, lbl, in_extension,
                                                                 self.workdir, lbl, out_extension)
            self.ExecuteCMD(cmd, False)
#            cp_res = os.system(cmd + " > /dev/null")
#            if cp_res != 0:
#                print "Error while executing:", cmd
#                print "Return code:", cp_res
#                sys.stderr.write(cmd + ": exit with code " + str(cp_res) + "\n")
#            
            # Copy the file into the cluster
            cmd = self.scpcmd + " %s/%s%s %s:%s/" % (self.local_workdir, lbl, 
                                                    in_extension, self.hostname, 
                                                    self.workdir)
            cp_res = self.ExecuteCMD(cmd, False)
            
            #cp_res = os.system(cmd + " > /dev/null")
            if not cp_res:
                print ("Error while executing:", cmd)
                print ("Return code:", cp_res)
                sys.stderr.write(cmd + ": exit with code " + str(cp_res) + "\n")
                continue
            
            tmt_str = ""
            if self.use_timeout:
                tmt_str = "timeout %d " % self.timeout
            app_list += "%s%s %s\n" % (tmt_str, new_mpicmd, binary)
            submitted.append(i)
            
        # Save the app list and copy it to the destination
        #app_list_name = "%s_app.list" % (label + "_" + str(indices[0]))
        #app_list_path = "%s/%s" % (self.local_workdir, app_list_name)
        #f = file(app_list_path, "w")
        #f.write(app_list)
        #f.close()
        
#        # Copy the app_list into the destination
#        cmd = self.scpcmd + " %s %s:%s" % (app_list_path, self.hostname, 
#                                           self.workdir)
#        cp_res = os.system(cmd)
#        if cp_res != 0:
#            print "Error while executing:", cmd
#            print "Return code:", cp_res
#            sys.stderr.write(cmd + ": exit with code " + str(cp_res) + "\n")
#            return results #[None] * N_structs
        
        
        # prepare the submission script
        submission = self.terminal + "\n"
        
        # Add the submission options
        if self.use_nodes:
            submission += "#%s %s%d\n" % (self.submit_name, self.v_nodes, self.n_nodes)
        if self.use_cpu:
            submission += "#%s %s%d\n" % (self.submit_name, self.v_cpu, new_ncpu)
        if self.use_time:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_time, self.time)
        if self.use_account:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_account, self.account_name)
        if self.use_memory:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_memory, self.ram)
        if self.use_partition:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_partition, self.partition_name)

        # Append the additional parameters
        for add_parameter in self.custom_params:
            submission += "#{} --{}={}\n".format(self.submit_name, add_parameter, self.custom_params[add_parameter])
        
        # Add the set -x option
        submission += "set -x\n"
        
        # Add the loading of the modules
        submission += self.load_modules + "\n"
        
        # Go to the working directory
        submission += "cd " + self.workdir + "\n"
        
        # Use the xargs trick
        #submission += "xargs -d " + r"'\n'" + " -L1 -P%d -a %s -- bash -c\n" % (n_togheder, 
        submission += app_list 
        
        # Copy the submission script
        sub_fpath = "%s/%s.sh" % (self.local_workdir, label + "_" + str(indices[0]))
        f = open(sub_fpath, "w")
        f.write(submission)
        f.close()
        cmd = self.scpcmd + " %s %s:%s" % (sub_fpath, self.hostname, self.workdir)
        cp_res = self.ExecuteCMD(cmd, False)
        #cp_res = os.system(cmd  + " > /dev/null")
        if not cp_res:
            print ("Error while executing:", cmd)
            print ("Return code:", cp_res)
            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
            return results#[None] * N_structs
        
        
        # Run the simulation
        cmd = self.sshcmd + " %s '%s %s/%s.sh'" % (self.hostname, self.submit_command, 
                                                   self.workdir, label+ "_" + str(indices[0]))
        cp_res = self.ExecuteCMD(cmd, False)
#        cp_res = os.system(cmd + " > /dev/null")
#        if cp_res != 0:
#            print "Error while executing:", cmd
#            print "Return code:", cp_res
#            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
#            
        
        # Collect the output back
        for i in submitted:
            # Prepare a typical label
            lbl = label + "_" + str(indices[i])
            out_filename = "%s/%s%s"% (self.workdir, lbl, out_extension)
            
            # Get the response
            cmd = self.scpcmd + " %s:%s %s/" % (self.hostname, out_filename,
                                                self.local_workdir)
            cp_res = self.ExecuteCMD(cmd, False)
            #cp_res = os.system(cmd + " > /dev/null")
            if not cp_res:
                #print "Error while executing:", cmd
                #print "Return code:", cp_res
                #sys.stderr.write(cmd + ": exit with code " + str(cp_res))
                continue
            
            # Get the results
            try:
                calc.set_label("%s/%s" % (self.local_workdir, lbl))
                calc.read_results()
                results[i] = copy.deepcopy(calc.results)
            except:
                pass
        
        return results
        
            
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


        # Append the additional parameters
        for add_parameter in self.custom_params:
            submission += "#{} --{}={}\n".format(self.submit_name, add_parameter, self.custom_params[add_parameter])
        
        # Add the loading of the modules
        submission += self.load_modules + "\n"
        
        # Go to the working directory
        submission += "cd " + self.workdir + "\n"
        
        # Get the real calculation command
        mpicmd = self.mpi_cmd.replace("NPROC", str(n_cpu))
        binary = self.binary.replace("NPOOL", str(npool)).replace("PREFIX", label)
        
        submission += mpicmd + " " + binary + "\n"
        
        # First of all clean eventually input/output file of this very same calculation
        cmd = self.sshcmd + " %s 'rm -f %s/%s%s %s/%s%s'" % (self.hostname, 
                                                             self.workdir, label, in_extension,
                                                             self.workdir, label, out_extension)
        self.ExecuteCMD(cmd, False)
#        cp_res = os.system(cmd)
#        if cp_res != 0:
#            print "Error while executing:", cmd
#            print "Return code:", cp_res
#            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
#        
        # Copy the input files into the target directory
        f = open("%s/%s.sh" % (self.local_workdir, label), "w")
        f.write(submission)
        f.close()
        cmd = self.scpcmd + " %s/%s.sh %s:%s" % (self.local_workdir, label, self.hostname, self.workdir)
        self.ExecuteCMD(cmd, False)
#        cp_res = os.system(cmd)
#        if cp_res != 0:
#            print "Error while executing:", cmd
#            print "Return code:", cp_res
#            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
#            
        cmd = self.scpcmd + " %s/%s%s %s:%s" % (self.local_workdir, label, in_extension, self.hostname, self.workdir)
        cp_res = self.ExecuteCMD(cmd, False)
        #cp_res = os.system(cmd)
        if not cp_res:
            #print "Error while executing:", cmd
            #print "Return code:", cp_res
            #sys.stderr.write(cmd + ": exit with code " + str(cp_res))
            return
        
        # Run the simulation
        cmd = self.sshcmd + " %s '%s %s/%s.sh'" % (self.hostname, self.submit_command, self.workdir, label)
        self.ExecuteCMD(cmd, False)
#        cp_res = os.system(cmd)
#        if cp_res != 0:
#            print "Error while executing:", cmd
#            print "Return code:", cp_res
#            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
            
        # Get the response
        cmd = self.scpcmd + " %s:%s/%s%s %s/" % (self.hostname, self.workdir, label, out_extension, 
                                                 self.local_workdir) 
        cp_res = self.ExecuteCMD(cmd, False)
        #cp_res = os.system(cmd)
        if not cp_res:
            print ("Error while executing:", cmd)
            print ("Return code:", cp_res)
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
        # Parse the namelist if needed
        if isinstance(namelist, str):
            namelist = CC.Methods.read_namelist(namelist)
        
        # Check if the cluster namelist is present
        if not __CLUSTER_NAMELIST__ in namelist.keys():
            raise ValueError("Error, the parsed namelist must contain %s" % __CLUSTER_NAMELIST__)
            
        c_info = namelist[__CLUSTER_NAMELIST__]
        keys = c_info.keys()
        
        # Check if there is an unknown key
        for k in keys:
            if not k in __CLUSTER_KEYS__:
                print ("Error with the key:", k)
                s = "Did you mean something like:" + str( difflib.get_close_matches(k, __CLUSTER_KEYS__))
                print (s)
                raise IOError("Error in cluster namespace: key '" + k +"' not recognized.\n" + s)
        
        # First of all, check if a template is present
        if __CLUSTER_TEMPLATE__ in keys:
            # Check if the environment variable has been defined
            fname = c_info[__CLUSTER_TEMPLATE__]
            if os.path.exists(fname):
                print ("Reading cluster info from:", os.path.abspath(fname))
                self.setup_from_namelist(CC.Methods.read_namelist(fname))
            elif __TEMPLATE_ENV__ in os.environ.keys():
                if os.path.exists(os.environ[__TEMPLATE_ENV__]):
                    newfname = os.environ[__TEMPLATE_ENV__] + "/" + fname
                    print ("Reading cluster info from:", os.path.abspath(newfname))
                    if not os.path.exists(newfname):
                        print ("Error, Environ variable", __TEMPLATE_ENV__, "exists, but no file", fname, "found.")
                        raise IOError("Error while reading the cluster template.")
                    
                    self.setup_from_namelist(CC.Methods.read_namelist(newfname))
                else:
                    sys.stderr.write("Error, Environ variable" + __TEMPLATE_ENV__ + "exists, but not the directory it is pointing to.\n")
                    raise IOError("Error while reading the cluster template.")
            else:
                sys.stderr.write("Error, no file " + fname + "found.\n")
                raise IOError("Error while reading the cluster template.")
        else:
            # Check if the required keywords are present
            for req_key in __CLUSTER_RKW__:
                if not req_key in keys:
                    raise IOError("Error, the cluster configuration namelist requires the keyword: '" + req_key + "'")
                    

        
        # Setup all the info
        if __CLUSTER_HOST__ in keys:
            self.hostname = c_info[__CLUSTER_HOST__]
            #print "HOST:", c_info[__CLUSTER_HOST__]
            
        if __CLUSTER_SSHCMD__ in keys:
            self.sshcmd = c_info[__CLUSTER_SSHCMD__]
        if __CLUSTER_SCPCMD__ in keys:
            self.scpcmd = c_info[__CLUSTER_SCPCMD__]
        
        if __CLUSTER_PWD__ in keys:
            self.pwd = c_info[__CLUSTER_PWD__]
            # Check if sshpass is present
            res = os.system("sshpass > /dev/null") 
            if res != 0:
                raise ValueError("Error, sshpass command not found, required to connect through password to server.")
            
            self.sshcmd = "sshpass -p '" + self.pwd + "' " + self.sshcmd
            self.scpcmd = "sshpass -p '" + self.pwd + "' " + self.scpcmd

        # Add the port
        if __CLUSTER_PORT__ in keys:
            # Check if the password has been setup
            self.sshcmd += " -p {:.0f}".format(c_info[__CLUSTER_PORT__])
            self.scpcmd += " -P {:.0f}".format(c_info[__CLUSTER_PORT__])
            
        if __CLUSTER_ACCOUNT__ in keys:
            self.account_name = c_info[__CLUSTER_ACCOUNT__]
            
        if __CLUSTER_MPICMD__ in keys:
            self.mpi_cmd = c_info[__CLUSTER_MPICMD__]
            
        if __CLUSTER_TIMEOUT__ in keys:
            self.use_timeout = True
            self.timeout = int(c_info[__CLUSTER_TIMEOUT__])
        if __CLUSTER_JOBNUMBER__ in keys:
            self.use_multiple_submission = True
            self.job_number = int(c_info[__CLUSTER_JOBNUMBER__])
            if self.job_number < 1:
                print ("Error, the number of job per batch must be >= 1")
                raise ValueError("Error in the %s input variable." % __CLUSTER_JOBNUMBER__)
        
        if __CLUSTER_NPARALLEL__ in keys:
            self.n_together_def = int(c_info[__CLUSTER_NPARALLEL__])
            
            if self.n_together_def < 1:
                print ("Error, the number of parallel jobs must be >= 1")
                raise ValueError("Error in the %s input variable." % __CLUSTER_NPARALLEL__)
            if self.n_together_def > self.job_number:
                print ("Error, the number of parallel runs must be <= than the number of runs per job")
                raise ValueError("Error, check the cluster keys %s and %s" % (__CLUSTER_NPARALLEL__, __CLUSTER_JOBNUMBER__))
            
            
        if __CLUSTER_TERMINAL__ in keys:
            self.terminal = "#!" + c_info[__CLUSTER_TERMINAL__]
        if __CLUSTER_SUBCMD__ in keys:
            self.submit_command = c_info[__CLUSTER_SUBCMD__]
        if __CLUSTER_SUBNAME__ in keys:
            self.submit_name = c_info[__CLUSTER_SUBNAME__]
        if __CLUSTER_VNODES__ in keys:
            self.v_nodes = c_info[__CLUSTER_VNODES__]
        if __CLUSTER_NNODES__ in keys:
            self.n_nodes = int(c_info[__CLUSTER_NNODES__])
        if __CLUSTER_UNODES__ in keys:
            self.use_nodes = c_info[__CLUSTER_UNODES__]
        if __CLUSTER_VCPU__ in keys:
            self.v_cpu = c_info[__CLUSTER_VCPU__]
        if __CLUSTER_NCPU__ in keys:
            self.n_cpu = int(c_info[__CLUSTER_NCPU__])
        if __CLUSTER_UCPU__ in keys:
            self.use_cpu = c_info[__CLUSTER_UCPU__]
        if __CLUSTER_VTIME__ in keys:
            self.v_time = c_info[__CLUSTER_VTIME__]
        if __CLUSTER_NTIME__ in keys:
            self.time = c_info[__CLUSTER_NTIME__]
        if __CLUSTER_UTIME__ in keys:
            self.use_time = c_info[__CLUSTER_UTIME__]
        if __CLUSTER_VMEM__ in keys:
            self.v_memory = c_info[__CLUSTER_VMEM__]
        if __CLUSTER_NMEM__ in keys:
            self.ram = c_info[__CLUSTER_NMEM__]
            self.use_memory = True
        #if __CLUSTER_UMEM__ in keys:
        #    self.use_memory = c_info[__CLUSTER_UMEM__]
        if __CLUSTER_VPART__ in keys:
            self.v_partition = c_info[__CLUSTER_VPART__]
        if __CLUSTER_NPART__ in keys:
            self.partition_name = c_info[__CLUSTER_NPART__]
            self.use_partition = True
        #if __CLUSTER_UPART__ in keys:
        #    self.use_partition = c_info[__CLUSTER_UPART__]
        if __CLUSTER_UACCOUNT__ in keys:
            self.use_account = c_info[__CLUSTER_UACCOUNT__]
        if __CLUSTER_VACCOUNT__ in keys:
            self.v_account = c_info[__CLUSTER_VACCOUNT__]
        if __CLUSTER_NPOOLS__ in keys:
            self.n_pool = int(c_info[__CLUSTER_NPOOLS__])
            
        if __CLUSTER_ATTEMPTS__ in keys:
            self.connection_attempts = int(c_info[__CLUSTER_ATTEMPTS__])
            
        if __CLUSTER_INITSCRIPT__ in keys:
            # Load the init script.
            # First lets parse the local environmental variables
            k_env = os.environ.keys()
            script_path = c_info[__CLUSTER_INITSCRIPT__]
            for key in k_env:
                script_path = script_path.replace("$" + key, os.environ[key])
            
            script_path = os.path.abspath(script_path)
            
            # Check if the file exists
            if not os.path.exists(script_path):
                raise IOError("Error, the provided script path %s does not exists." % script_path)
            
            # Read the script and store the contnet
            f = open(script_path, "r")
            self.load_modules = f.read()
            f.close()
        
        if __CLUSTER_MAXRECALC__ in keys:
            self.max_recalc = int(c_info[__CLUSTER_MAXRECALC__])
        if __CLUSTER_BATCHSIZE__ in keys:
            self.batch_size = int(c_info[__CLUSTER_BATCHSIZE__])
        if __CLUSTER_LOCALWD__ in keys:
            wdir = c_info[__CLUSTER_LOCALWD__]
            
            # Parse the local environmental variables
            for ekey in os.environ.keys():
                wdir = wdir.replace("$" + ekey, os.environ[ekey])
                
            self.local_workdir = wdir
            
        
        if __CLUSTER_BINARY__ in keys:
            print ("Binary before parsing:")
            print (c_info[__CLUSTER_BINARY__])
            self.binary = self.parse_string(c_info[__CLUSTER_BINARY__])
            print ("Cluster binary setted to:")
            print (self.binary)
            
        # If all the cluster has been setted, setup the working directory
        if __CLUSTER_WORKDIR__ in keys:
            self.workdir = c_info[__CLUSTER_WORKDIR__]
            self.setup_workdir()

            
    def setup_workdir(self, verbose = True):
        """
        SETUP THE WORKING DIRECTORY
        ===========================
        
        Parse the line contained in self.workdir in the claster to get working directory.
        It needs that the communication with the cluster has been correctly setted up.
        
        It will parse correctly environmental variables of the cluster.
        """
        workdir = self.parse_string(self.workdir)
        
        sshcmd = self.sshcmd + " %s 'mkdir -p %s'" % (self.hostname, 
                                                      workdir)
        
        self.ExecuteCMD(sshcmd)
#        
#        retval = os.system(sshcmd)
#        if retval != 0:
#            raise IOError("Error, while executing cmd: " + sshcmd)
#        
#        if verbose:
#            print "Cluster workdir setted to:"
#            print workdir
        
        self.workdir = workdir
    
    def parse_string(self, string):
        """
        PARSE STRING
        ============
        
        Parse the given string on the cluster. 
        It can be used to resolve environmental variables defined in the cluster.
        
        It will execute on the cluster the command:
            echo "string"
        and return the result of the cluster.
        
        Parameter
        ---------
            string :
                String to be parsed in the cluster.
        Result
        ------
            string :
                The same as input, but with the cluster environmental variables correctly
                parsed.
        """
        
        # Open a pipe with the server
        # Use single ' to avoid string parsing by the local terminal
        cmd = "%s %s 'echo \"%s\"'" % (self.sshcmd, self.hostname, string)
        #print cmd
        
        status, output = self.ExecuteCMD(cmd, return_output = True)
#        
#        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
#        output, err = p.communicate()
#        status = p.wait()
#        output = output.strip()
#
#        if not err is None: 
#            sys.stderr.write(err)
#        if status != 0:
#            print "Command:", cmd
#            print "Error, status:", status
#            raise ValueError("Error, while connecting with the server.")
        
        return str(output)
            
    def compute_ensemble_batch(self, ensemble, ase_calc, get_stress = True, timeout=None):
        """
        RUN THE ENSEMBLE WITH BATCH SUBMISSION
        ======================================
        """
        
        # Track the remaining configurations
        success = [False] * ensemble.N
        
        # Setup if the ensemble has the stress
        ensemble.has_stress = get_stress
        
        # Check if the working directory exists
        if not os.path.isdir(self.local_workdir):
            os.makedirs(self.local_workdir)
    
        
        # Get the expected number of batch
        num_batch_offset = int(ensemble.N / self.batch_size)
        
        def compute_single_jobarray(jobs_id, calc):
            structures = [ensemble.structures[i] for i in jobs_id]
            n_together = min(len(structures), self.n_together_def)
            results = self.batch_submission(structures, calc, jobs_id, ".pwi",
                                            ".pwo", "ESP", n_together)
            
            for i, res in enumerate(results):
                num = jobs_id[i]
                
                # Check if the run was good
                try:
                    resk = res.keys()
                except:
                    continue
                check_e = "energy" in resk
                check_f = "forces" in resk
                check_s = "stress" in resk
                
                is_success =  check_e and check_f
                if get_stress:
                    is_success = is_success and check_s
                
                if not is_success:
                    continue
                
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
                success[num] = is_success
        
        # Run until some work has not finished
        recalc = 0
        while np.sum(np.array(success, dtype = int) - 1) != 0:
            threads = []
            
            # Get the remaining jobs
            false_mask = np.array(success) == False
            false_id = np.arange(ensemble.N)[false_mask]
            
            count = 0
            # Submit in parallel
            jobs = [false_id[i : i + self.job_number] for i in range(0, len(false_id), self.job_number)]
            
            for job in jobs:
                # Submit only the batch size
                if count > self.batch_size:
                    break
                t = threading.Thread(target = compute_single_jobarray, args=(job, ase_calc, ))
                t.start()
                threads.append(t)
                count += 1
            
            # Wait until all the job have finished
            for t in threads:
                t.join(timeout)
            
            recalc += 1
            if recalc > num_batch_offset + self.max_recalc:
                print ("Expected batch ordinary resubmissions:", num_batch_offset)
                raise ValueError("Error, resubmissions exceeded the maximum number of %d" % self.max_recalc)
                break
            
            
    
    def compute_ensemble(self, ensemble, ase_calc, get_stress = True, timeout=None):
        """
        RUN THE WHOLE ENSEMBLE ON THE CLUSTER
        =====================================
        
        Parameters
        ----------
            ensemble :
                The ensemble to be runned.
        """
        
        # Check if the compute_ensemble batch must be done
        if self.job_number != 1:
            self.compute_ensemble_batch(ensemble, ase_calc, get_stress, timeout)
            return
        
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
                print ("Expected batch ordinary resubmissions:", num_batch_offset)
                raise ValueError("Error, resubmissions exceeded the maximum number of %d" % self.max_recalc)
                break
            