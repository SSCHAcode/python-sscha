# -*- coding: utf-8 -*-
from __future__ import print_function
import sys, os
import subprocess
import threading
import copy
import time, datetime
import tarfile

__DIFFLIB__ = False
try:
    __DIFFLIB__ = True
    import difflib
except:
    pass

import numpy as np
import time, datetime

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
__CLUSTER_VQOS__ = "v_qos"
__CLUSTER_NQOS__ = "qos_name"
__CLUSTER_UQOS__ = "use_qos"
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
                    __CLUSTER_UPART__, __CLUSTER_VQOS__, __CLUSTER_NQOS__, __CLUSTER_UQOS__,
                    __CLUSTER_INITSCRIPT__, __CLUSTER_MAXRECALC__, __CLUSTER_BATCHSIZE__,
                    __CLUSTER_LOCALWD__, __CLUSTER_VACCOUNT__, __CLUSTER_UACCOUNT__, __CLUSTER_SSHCMD__,
                    __CLUSTER_SCPCMD__, __CLUSTER_WORKDIR__, __CLUSTER_TIMEOUT__,
                    __CLUSTER_JOBNUMBER__, __CLUSTER_NPARALLEL__, __CLUSTER_NPOOLS__,
                    __CLUSTER_ATTEMPTS__, __CLUSTER_PORT__]


SPECIAL_SYMBOLS = ["$", ";", "|"]


def parse_symbols(string):
    r"""
    REPLACE SPECIAL SYMBOLS
    =======================

    In a string that must be used on a shell command, replace all the special symbols
    in a way that they are correctly passed through the SHELL.

    for example $USER  => \$USER

    Parameters
    ----------
        string : str
            The string to be parsed

    Results
    -------
        output : str
            The string after the symbols are replaced
    """
    new_str = string
    for symbol in SPECIAL_SYMBOLS:
        new_str = new_str.replace(symbol, "\\" + symbol)
    return new_str




class Cluster(object):

    def __init__(self, hostname=None, pwd=None, extra_options="", workdir = "",
                 account_name = "", partition_name = "", qos_name = "", binary="pw.x -npool NPOOL -i PREFIX.pwi > PREFIX.pwo",
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
            qos_name:
                The QOS name for the job submission
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
        self.qos_name = qos_name
        if qos_name:
            self.use_qos = True

        self.binary = binary
        self.mpi_cmd = mpi_cmd

        self.workdir=""
        self.submit_command="sbatch"
        self.submit_name="SBATCH"
        self.terminal="/bin/bash"
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
        self.v_qos = "--qos="
        self.use_qos = False
        self.timeout = 1000
        self.use_timeout = False

        # Check the status of the job every TOT seconds
        self.check_timeout = 300
        self.nonblocking_command = True # True if you use a different version of slurm that does not accept blocking commands

        # Enforce ssh to open a shell for each command in the cluster
        self.use_active_shell = False
        # If the following is true, then use the active shell always
        self.use_active_shell_for_parsing = False

        # This is the number of configurations to be computed for each jub submitted
        # This times the self.batch_size is the total amount of configurations submitted toghether
        self.job_number = 1
        self.n_together_def = 1
        self.use_multiple_submission = False

        # If true, add the set -x option at the beggining of the script
        # This options makes the system print on stdout all executed commands.
        # Very usefull to debug if something goes wrong.
        # On some clusters may cause malfunctions.
        self.add_set_minus_x = True


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

        # The default qos in which to submit calculations
        self.qos_name = ""

        # Still unused
        self.prefix_name = "prefix" # Variable in the calculator for differentiating the calculations

        # This directory is used to work with clusters.
        self.local_workdir = "cluster_work/"

        # The batch size is the maximum number of job to be submitted together.
        # The new jobs will be submitted only after a batch is compleated
        # Useful if the cluster has a limit for the maximum number of jobs allowed.
        self.batch_size = 1000

        # This could be a function that generates for each input file additional
        # text in the submission file.
        # Usefull if you want to copy things in a different directory to run the calculation on the cluster.
        self.additional_script_parameters = None



        # Allow to setup additional custom extra parameters
        self.custom_params = {}

        self.lock = None  # Use to lock the threads

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

            # Setting the account or partition name will automatically result in
            # activating the corresponding flags
            if name.endswith("_name"):
                key = "use_{}".format(name.split("_")[0])
                self.__dict__[key] = True
        else:
            super(Cluster, self).__setattr__(name, value)



    def copy_file(self, source, destination, server_source = False, server_dest = True, raise_error=False, **kwargs):
        """
        COPY A FILE
        ===========

        This function copies a file or directory from the source to the destination.
        The destination is on the cluster, the source is in the local machine.

        It uses scp to perform the copy.
        Alternative implementations can be performed by overloading this function.

        args and kwargs are passed to the ExecuteCMD function.

        The result is the output of the ExecuteCMD function.

        Parameters
        ----------
            source : string
                The source file to be copied
            destination : string
                The destination file
            server_source : bool
                If true, the source is on the server
            server_dest : bool
                If true, the destination is on the server
            raise_error : bool
                If True, raises an error upon failure
        """
        server_path = "%s:" % self.hostname
        source_path = f"{source}"
        dest_path = f"{destination}"

        if server_source:
            source_path = server_path + source_path 
        if server_dest:
            dest_path = server_path + dest_path 
        
        cmd = self.scpcmd + f" {source_path} {dest_path}"
        result = self.ExecuteCMD(cmd, raise_error = raise_error, **kwargs)
        return result


    def ExecuteCMD(self, cmd, raise_error = False, return_output = False, on_cluster = False, use_active_shell = False):
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
            on_cluster : bool
                If true, the command is executed directly on the cluster through ssh
            use_active_shell : bool
                If true, the command is executed in a new shell on the cluster
                This is usefull if the command is a script that must be executed
                in a new shell, or if the command requires .bashrc to be sourced.

        Returns
        -------
            success : bool
                If True, the command has been executed with success,
                False otherwise
            output : string
                Returned only if return_output is True
        """

        if on_cluster:
            cmd = self.sshcmd + " {} '{}'".format(self.hostname, cmd)
            if use_active_shell:
                cmd = "{ssh} {host} -t '{shell} --login -c \"{command}\"'".format(ssh = self.sshcmd, 
                         host = self.hostname, 
                         command = parse_symbols(cmd), 
                         shell = self.terminal)





        success = False
        output = ""
        for i in range(self.connection_attempts):
            # Subprocess opening
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True)
            output, err = p.communicate()
            if not isinstance(output, str):
                try:
                    output = str(output.decode("utf-8"))
                except Exception as e:
                    sys.stderr.write("Error in the following command:\n")
                    sys.stderr.write(e)
                    sys.stderr.write("\n")
                    sys.stderr.flush()

            output = output.strip()
            status = p.wait()


            print('THREAD {} EXECUTE COMMAND: {}'.format(threading.get_native_id(), cmd))
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

        #cmd = self.sshcmd + " %s 'echo ciao'" % self.hostname
        cmd = "echo ciao"
        status, output = self.ExecuteCMD(cmd, return_output = True, on_cluster = True)

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

    def create_submission_script(self, labels):
        """
        CREATE THE SUBMISSION SCRIPT
        ===========================================

        This is a function that is general and does not depend on the specific
        calculator. It is usefull to create the header of the submission script.

        Parameters
        ----------
            labels : list
                It is a list of the labels of the calculations to be done.

        Returns
        -------
            submission_header : string
                The text of the submission header.
        """

        # prepare the submission script
        submission = "#!" + self.terminal + "\n"



        # Add the submission options
        if self.use_nodes:
            submission += "#%s %s%d\n" % (self.submit_name, self.v_nodes, self.n_nodes)
        if self.use_cpu:
            submission += "#%s %s%d\n" % (self.submit_name, self.v_cpu, self.n_cpu)
        if self.use_time:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_time, self.time)
        if self.use_account:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_account, self.account_name)
        if self.use_memory:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_memory, self.ram)
        if self.use_partition:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_partition, self.partition_name)
        if self.use_qos:
            submission += "#%s %s%s\n" % (self.submit_name, self.v_qos, self.qos_name)

        # Append the additional parameters
        for add_parameter in self.custom_params:

            adder_string = "--{}".format(add_parameter)
            if add_parameter.startswith("-"):
                adder_string = add_parameter

            if self.custom_params[add_parameter] is None:
                submission += "#{} {}\n".format(self.submit_name, adder_string)
            else:
                submission += "#{} {}={}\n".format(self.submit_name, adder_string, self.custom_params[add_parameter])


        # Add the set -x option
        if self.add_set_minus_x:
            submission += "set -x\n"

        # Add the loading of the modules
        submission += self.load_modules + "\n"

        # Go to the working directory
        submission += "cd " + self.workdir + "\n"

        # If any, apply the extra text before and after the calculation
        other_input = ""
        other_output = ""
        if (self.additional_script_parameters is not None):
            other_input, other_output = self.additional_script_parameters(labels)

        submission += other_input

        # Use the xargs trick
        #submission += "xargs -d " + r"'\n'" + " -L1 -P%d -a %s -- bash -c\n" % (n_togheder,
        for i, lbl in enumerate(labels):
            submission += self.get_execution_command(lbl)


        submission += other_output

        return submission

    def get_execution_command(self, label):
        """
        GET THE EXECUTION COMMAND
        =========================

        Return the command used in the submission script to actually execute the calculation.

        Parameters
        ----------
            label : string
                The label of the calculation

        Returns
        -------
            commnad : string
                The command to be appended to the submission script
        """

        # Get the MPI command replacing NPROC
        new_mpicmd  = self.mpi_cmd.replace("NPROC", str(self.n_cpu))

        # Replace the NPOOL variable and the PREFIX in the binary
        binary = self.binary.replace("NPOOL", str(self.n_pool)).replace("PREFIX", label)


        tmt_str = ""
        if self.use_timeout:
            tmt_str = "timeout %d " % self.timeout
        return "%s%s %s\n" % (tmt_str, new_mpicmd, binary)

    def prepare_input_file(self, structures, calc, labels):
        """
        PREPARE THE INPUT FILE
        ======================

        This is specific for quantum espresso and must be inherit and replaced for
        other calculators.

        This crates the input files in the local working directory
        self.local_workdir and it returns the list of all the files generated.


        Parameters
        ----------
            structures : List of cellconstructor.Structure.Structure
                The atomic structures.
            calc : the ASE or CellConstructor calculator.
                In this case, it works with quantum espresso
            labels : List of strings
                The unique name of this calculation

        Returns
        -------
            List_of_input : list
                List of strings containing all the input files
            List_of_output : list
                List of strings containing the output files expected
                for the calculation
        """

        # Prepare the input file
        list_of_inputs = []
        list_of_outputs = []
        for i, (label, structure) in enumerate(zip(labels, structures)):
            # Avoid thread conflict
            self.lock.acquire()

            try:
                calc.set_directory(self.local_workdir)
                calc.set_label(label)
                calc.write_input(structure)

                print("[THREAD {}] LBL: {} | PREFIX: {}".format(threading.get_native_id(), label, calc.input_data["control"]["prefix"]))

                #ase.io.write("%s/%s.pwi"% (self.local_workdir, label),
                #                atm, **calc.parameters)

                input_file = '{}.pwi'.format(label)
                output_file = '{}.pwo'.format(label)

                list_of_inputs.append(input_file)
                list_of_outputs.append(output_file)
            except Exception as e:
                MSG = '''
Error while writing input file {}.

Error message:
'''.format(label)
                MSG += str(repr(e))
                print(MSG)


            # Release the lock on the threads
            self.lock.release()


        return list_of_inputs, list_of_outputs

    def clean_localworkdir(self):
        """
        CLEAN THE LOCAL WORKDIR FROM ALL INPUT/OUTPUT FILES IN TAR
        ===================================================
        """

        all_files = [x for x in os.listdir(self.local_workdir) if x.startswith('input') and x.endswith('.tar')]
        for f in all_files:
            os.remove(os.path.join(self.local_workdir, f))

    def copy_files(self, list_of_input, list_of_output, to_server):
        """
        COPY INPUT/OUTPUT FILES FROM/TO THE SERVER
        ==========================================

        This function copies the input files toward the HPC if to_server is True
        or retrive the output files from the HPC if to_server is False

        Parameters
        ----------
            list_of_input : list
                List of the path to the input files to be copied into the server.
            list_of_output : list
                List of the output files to be copied from the HPC server.
                It is needed also when submitting the input file from the server,
                as files that are named in the same way will be cleaned.
                (No risk to mistake them with the actual result of the calculation)
            to_server : bool
                If true, we copy the input files into the server from the local machine.
                If false, we copy the output files from the server to the local machine.
        """

        thread_id = threading.get_native_id()

        if to_server:
            # Compress the files
            tar_name = 'inputs_id{}.tar'.format(thread_id)
            tar_file = os.path.join(self.local_workdir, tar_name)
            cmd = 'tar -cf {} -C {}'.format(tar_file, self.local_workdir)

            # Remove the old tar file and the old input/output files\
            rm_cmd = ''
            #rm_cmd = 'rm -f {}; '.format(os.path.join(self.workdir, tar_name))

            for i, fname in enumerate(list_of_input):
                cmd += ' ' +  fname

                # Remove all the previous output files on cluster
                rm_cmd = 'rm -f {}; '.format(os.path.join(self.workdir, fname))



            os.system(cmd)

            if not os.path.exists(tar_file):
                raise IOError("""
Error, for some reason I'm unable to generate the tar.
       command = {}
""".format(cmd))

            for fname in list_of_output:
                rm_cmd += 'rm -f {}; '.format(os.path.join(self.workdir, fname))

                # Check if the output files are already present.
                # In this case, remove them to avoid conflict
                if os.path.exists(os.path.join(self.local_workdir, fname)):
                    os.remove(os.path.join(self.local_workdir, fname))


            # Clean eventually input/output file of this very same calculation
            self.ExecuteCMD(rm_cmd, False, on_cluster = True)
#            cp_res = os.system(cmd + " > /dev/null")
#            if cp_res != 0:
#                print "Error while executing:", cmd
#                print "Return code:", cp_res
#                sys.stderr.write(cmd + ": exit with code " + str(cp_res) + "\n")
#

            # Copy the file into the cluster
            cp_res = self.copy_file(tar_file, self.workdir, raise_error=False)
            #cmd = self.scpcmd + " %s %s:%s/" % (tar_file, self.hostname, self.workdir)
            #cp_res = self.ExecuteCMD(cmd, False)
            if not cp_res:
                print ("Error while executing:", cmd)
                print ("Return code:", cp_res)
                sys.stderr.write(cmd + ": exit with code " + str(cp_res) + "\n")
                print("""
Error while connecting to the cluster to copy the files:
{}
""".format(list_of_input))
                return cp_res

            # Unpack the input files and remove the archive
            decompress = 'cd {}; tar xf {};'.format(self.workdir, tar_name)
            #cmd = self.sshcmd + " %s '%s'" % (self.hostname, decompress)
            #cp_res = self.ExecuteCMD(cmd, False)
            cp_res = self.ExecuteCMD(decompress, False, on_cluster = True)
            if not cp_res:
                print ("Error while executing:", cmd)
                print ("Return code:", cp_res)
                sys.stderr.write(cmd + ": exit with code " + str(cp_res) + "\n")
            return cp_res
            #cp_res = os.system(cmd + " > /dev/null")
        else:
            # Compress all the output files at once
            tar_name = 'outputs_id{}.tar'.format(thread_id)
            tar_command = 'tar cf {} '.format(tar_name)
            for output in list_of_output:
                tar_command += ' ' + output

            compress_cmd = 'cd {}; {}'.format(self.workdir, tar_command)

            # cmd = self.sshcmd + " %s '%s'" % (self.hostname, compress_cmd)
            # cp_res = self.ExecuteCMD(cmd, False)
            cp_res = self.ExecuteCMD(compress_cmd, False, on_cluster = True)
            if not cp_res:
                print ("Error while compressing the outputs:", cmd, list_of_output, "\nReturn code:", cp_res)
                #return cp_res


            # Copy the tar from the server to the local and unpack
            # cmd = self.scpcmd + "%s:%s %s/" % (self.hostname, os.path.join(self.workdir, tar_name), self.local_workdir)
            # cp_res = self.ExecuteCMD(cmd, False)
            cp_res = self.copy_file(os.path.join(self.workdir, tar_name), self.local_workdir, raise_error=False, server_source=True, server_dest=False)
            if not cp_res:
                print ("Error while executing:", cmd)
                print ("Return code:", cp_res)
                sys.stderr.write(cmd + ": exit with code " + str(cp_res) + "\n")
                print("""
Error while connecting to the cluster to copy the files:
{}
""".format(list_of_input))
                return cp_res

            # Extract the tar file
            os.system('tar xf {} -C {}'.format(os.path.join(self.local_workdir, tar_name), self.local_workdir))
            #os.remove(os.path.join(self.local_workdir, tar_name))

            return cp_res




    def submit(self, script_location):
        """
        SUBMIT THE CALCULATION
        ======================

        Submit the calculation. Compose the command into a cmd variable, then submit it through:

        .. code ::

            return self.ExecuteCMD(cmd, True, return_output=True)



        Parameters
        ----------
            script_localtion : string
                Path to the submission script inside the cluster.

        Results
        -------
            success : bool
                Result of the execution of the submission command.
                It is what returned from self.ExecuteCMD(cmd, False)
        """

        cmd = f"{self.submit_command} {script_location}"
        #cmd = "{ssh} {host} '{submit_cmd} {script}'".format(ssh = self.sshcmd, host = self.hostname,
        #                 submit_cmd = self.submit_command, script = script_location)
        
        result, output = self.ExecuteCMD(cmd, False, return_output=True, on_cluster = True,
                                        use_active_shell = self.use_active_shell)


        # if self.use_active_shell:
        #     cmd = "{ssh} {host} -t '{shell} --login -c \"{submit_cmd} {script}\"'".format(ssh = self.sshcmd,
        #                  host = self.hostname,
        #                  submit_cmd = self.submit_command, script = script_location,
        #                  shell = self.terminal)



        #cmd = self.sshcmd + " %s '%s %s/%s.sh'" % (self.hostname, self.submit_command,
        #                                           self.workdir, label+ "_" + str(indices[0]))

        return result, output

    def get_output_path(self, label):
        """
        Given the label of the submission, retrive the path of all the output files of that calculation
        """

        out_filename = os.path.join(self.workdir, label + ".pwo")
        return [out_filename]

    def read_results(self, calc, label):
        """
        Return a dictionary of the computed property for the given calculation label
        """
        print("[READ RESULTS] THREAD ID {} ENTERED".format(threading.get_native_id()))
        #calc.set_label("%s/%s" % (self.local_workdir, label))

        # Perform a safe thread read of the results
        calc.set_directory(self.local_workdir)
        calc.set_label(label)
        calc.read_results()
        #calc.structure.save_scf("thread_{}_justreaded_{}.scf".format(threading.get_native_id(), label))
        results = copy.deepcopy(calc.results)

        if results is not None:
            results["structure"] = calc.structure.copy()  # Report also the structure to check consistency

        return results

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
        results = [None] * N_structs
        submitted = []
        submission_labels = []



        for i in range(N_structs):
            # Prepare a typical label
            lbl = label + "_" + str(indices[i])
            submission_labels.append(lbl)
            submitted.append(i)


        # Create the input files
        input_files, output_files =  self.prepare_input_file(list_of_structures, calc, submission_labels)

        # Create the submission script
        submission = self.create_submission_script(submission_labels)

        # Add the submission script to the input files
        sub_name =  label + "_" + str(indices[0]) + '.sh'
        sub_fpath = os.path.join(self.local_workdir, sub_name)
        f = open(sub_fpath, "w")
        f.write(submission)
        f.close()
        input_files.append(sub_name)

        if not self.copy_files(input_files, output_files, to_server = True):
            # Submission failed
            return submitted, indices, label


        # Run the simulation
        sub_script_loc = os.path.join(self.workdir, sub_name)
        cp_res, submission_output = self.submit(sub_script_loc)

        if self.nonblocking_command:
            job_id = self.get_job_id_from_submission_output(submission_output)

            now = datetime.datetime.now()
            sys.stderr.write("{}/{}/{} - {}:{}:{} | submitted job id {} ({})\n".format(now.year, now.month, now.day, now.hour, now.minute, now.second, job_id, submission_output))
            sys.stderr.flush()
            time.sleep(self.check_timeout)

            while not self.check_job_finished(job_id):
                time.sleep(self.check_timeout)

        # Collect back the output
        if not self.copy_files(input_files, output_files, to_server = False):
            # Submission failed
            print('[SUBMISSION {}] FAILED RETRIVING OUTPUT'.format(threading.get_native_id()))
        else:
            print('[SUBMISSION {}] GOT OUTPUT'.format(threading.get_native_id()))


        return submitted, indices, label

    def collect_results(self, calc, submitted, indices, label):
        """
        Collect back all the results from the submitted data.
        This operation needs to be performed in thread safe mode
        (all threads must be locked between this operation and when the results are actually assigned to the ensemble
        to avoid conflicts)
        """
        # Read the output files
        results = [None] * len(submitted)
        for i in submitted:
            # Prepare a typical label
            lbl = label + "_" + str(indices[i])

            # Get the results
            try:
                results[i] = self.read_results(calc, lbl)
            except FileNotFoundError:
                sys.stderr.write("JOB {} | {} resulted in error:\n".format(i, lbl))
                sys.stderr.write('File not found!\n')
                sys.stderr.flush()
            except Exception as e:
                sys.stderr.write("JOB {} | {} resulted in error:\n".format(i, lbl))
                print(e, file=sys.stderr)
                sys.stderr.flush()


        print('[SUBMISSION {}] GOT RESULTS:\n{}'.format(threading.get_native_id(), results))
        return results

    def get_job_id_from_submission_output(self, output):
        """
        GET THE JOB ID

        Retreive the job id from the output of the submission.
        This depends on the software employed. It works for slurm.

        Returns None if the output contains an error
        """
        print('INSIDE GET JOB ID WITH {}'.format(output))


        try:
            id = output.split()[-1]

            print('Understanding output "{}" => {}'.format(output, id))
            return id
        except:
            print("Error, expected a standard output, but the result of the submission was: {}".format(output))
            return None

    def check_job_finished(self, job_id, verbose = True):
        """
        Check if the job identified by the job_id is finished

        Parameters
        ----------
            job_id : string
                The string that identifies uniquely the job
        """

        status, output = self.ExecuteCMD("squeue -u $USER", False, return_output = True, on_cluster = True, )
        lines = output.split("\n")

        if len(lines):
            for l in lines:
                l = l.strip()
                if not l:
                    if verbose:
                        now = datetime.datetime.now()
                        sys.stderr.write("{}/{}/{} - {}:{}:{} | job {}: No response from the server \n".format(now.year, now.month, now.day, now.hour, now.minute, now.second, job_id))
                        sys.stderr.flush()
                    return False

                data = l.split()
                if len(data) == 0:
                    if verbose:
                        now = datetime.datetime.now()
                        sys.stderr.write("{}/{}/{} - {}:{}:{} | job {}: No response from the server \n".format(now.year, now.month, now.day, now.hour, now.minute, now.second, job_id))
                        sys.stderr.flush()
                    return False
                elif data[0] == job_id:
                    if verbose:
                        now = datetime.datetime.now()
                        sys.stderr.write("{}/{}/{} - {}:{}:{} | job {} still running\n".format(now.year, now.month, now.day, now.hour, now.minute, now.second, job_id))
                        sys.stderr.flush()
                    return False

            # If I'm here it means I did not find the job, but the command returned at least 1 line (so it was correctly executed).
            if verbose:
                now = datetime.datetime.now()
                sys.stderr.write("{}/{}/{} - {}:{}:{} | job {} finished\n".format(now.year, now.month, now.day, now.hour, now.minute, now.second, job_id))
                sys.stderr.flush()
            return True
        if verbose:
            now = datetime.datetime.now()
            sys.stderr.write("{}/{}/{} - {}:{}:{} | error while interrogating the cluster for job {}\n".format(now.year, now.month, now.day, now.hour, now.minute, now.second, job_id))
            sys.stderr.flush()
        return False


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

            adder_string = "--{}".format(add_parameter)
            if add_parameter.startswith("-"):
                adder_string = add_parameter

            if self.custom_params[add_parameter] is None:
                submission += "#{} {}\n".format(self.submit_name, adder_string)
            else:
                submission += "#{} {}={}\n".format(self.submit_name, adder_string, self.custom_params[add_parameter])

        # Add the loading of the modules
        submission += self.load_modules + "\n"

        # Go to the working directory
        submission += "cd " + self.workdir + "\n"

        # Get the real calculation command
        mpicmd = self.mpi_cmd.replace("NPROC", str(n_cpu))
        binary = self.binary.replace("NPOOL", str(npool)).replace("PREFIX", label)

        submission += mpicmd + " " + binary + "\n"

        # First of all clean eventually input/output file of this very same calculation
        cmd = "rm -f %s/%s%s %s/%s%s" % (self.workdir, label, in_extension,
                                         self.workdir, label, out_extension)
        self.ExecuteCMD(cmd, False, on_cluster = True)
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
        #cmd = self.scpcmd + " %s/%s.sh %s:%s" % (self.local_workdir, label, self.hostname, self.workdir)
        self.copy_file("%s/%s.sh" % (self.local_workdir, label), self.workdir, server_source=False, server_dest=True)
        #self.ExecuteCMD(cmd, False)
#        cp_res = os.system(cmd)
#        if cp_res != 0:
#            print "Error while executing:", cmd
#            print "Return code:", cp_res
#            sys.stderr.write(cmd + ": exit with code " + str(cp_res))
#
        cp_res = self.copy_file("%s/%s%s" % (self.local_workdir, label, in_extension), self.workdir, server_source=False, server_dest=True, raise_error=False)
        #cmd = self.scpcmd + " %s/%s%s %s:%s" % (self.local_workdir, label, in_extension, self.hostname, self.workdir)
        #cp_res = self.ExecuteCMD(cmd, False)
        #cp_res = os.system(cmd)
        if not cp_res:
            #print "Error while executing:", cmd
            #print "Return code:", cp_res
            #sys.stderr.write(cmd + ": exit with code " + str(cp_res))
            return

        # Run the simulation
        cmd = "%s %s/%s.sh" % (self.submit_command, self.workdir, label)
        #cmd = self.sshcmd + " %s '%s %s/%s.sh'" % (self.hostname, self.submit_command, self.workdir, label)
        self.ExecuteCMD(cmd, False, on_cluster = True)
#        cp_res = os.system(cmd)
#        if cp_res != 0:
#            print "Error while executing:", cmd
#            print "Return code:", cp_res
#            sys.stderr.write(cmd + ": exit with code " + str(cp_res))

        # Get the response
        #cmd = self.scpcmd + " %s:%s/%s%s %s/" % (self.hostname, self.workdir, label, out_extension,
        #self.local_workdir)
        #cp_res = self.ExecuteCMD(cmd, False)
        cp_res = self.copy_file("%s/%s%s" % (self.workdir, label, out_extension), self.local_workdir, server_source=True, server_dest=False)
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
        if __CLUSTER_VQOS__ in keys:
            self.v_qos = c_info[__CLUSTER_VQOS__]
        if __CLUSTER_NQOS__ in keys:
            self.qos_name = c_info[__CLUSTER_NQOS__]
            self.use_qos = True
        if __CLUSTER_UQOS__ in keys:
            self.use_qos = c_info[__CLUSTER_UQOS__]
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

        cmd = "mkdir -p %s" % workdir
        # sshcmd = self.sshcmd + " %s 'mkdir -p %s'" % (self.hostname,
        #                                               workdir)

        self.ExecuteCMD(cmd, raise_error= True, on_cluster=True)
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
        #cmd = "%s %s 'echo \"%s\"'" % (self.sshcmd, self.hostname, string)
        cmd = f"echo \"{string}\""

        status, output = self.ExecuteCMD(cmd, return_output = True, raise_error= True, use_active_shell = self.use_active_shell_for_parsing, on_cluster = True)

                #print cmd

        #print(cmd)

        #status, output = self.ExecuteCMD(cmd, return_output = True, raise_error= True)
        
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

    def compute_ensemble_batch(self, ensemble, cellconstructor_calc, get_stress = True, timeout=None):
        """
        RUN THE ENSEMBLE WITH BATCH SUBMISSION
        ======================================
        """

        # Track the remaining configurations
        success = [False] * ensemble.N

        # Setup if the ensemble has the stress
        ensemble.has_stress = get_stress
        #ensemble.all_properties = [None] * ensemble.N

        # Check if the working directory exists
        if not os.path.isdir(self.local_workdir):
            os.makedirs(self.local_workdir)


        # Get the expected number of batch
        num_batch_offset = int(ensemble.N / self.batch_size)

        def compute_single_jobarray(jobs_id, calc):
            structures = [ensemble.structures[i].copy() for i in jobs_id]
            n_together = min(len(structures), self.n_together_def)
            subs, indices, labels = self.batch_submission(structures, calc, jobs_id, ".pwi",
                                            ".pwo", "ESP", n_together)

            # Thread safe operation
            self.lock.acquire()
            print("[THREAD {}] submitted calculations: {}".format(threading.get_native_id(), indices))
            results = self.collect_results(calc, subs, indices, labels)

            for i, res in enumerate(results):
                print("[THREAD {}] ADDING RESULT {} = {}".format(threading.get_native_id(), jobs_id[i], res))
                num = jobs_id[i]

                if res is None:
                    continue

                # Check if the run was good
                check_e = "energy" in res
                check_f = "forces" in res
                check_s = "stress" in res

                # Check the structure
                if "structure" in res:
                    error_struct = np.linalg.norm(ensemble.structures[jobs_id[i]].coords.ravel() - res["structure"].coords.ravel())
                    if error_struct > 1e-2:
                        print("ERROR IDENTIFYING STRUCTURE!")
                        MSG = """
                            Error in thread {}.
                            Displacement between the expected structure {}
                            and the one readed from the calculator
                            is of {} A.
                        """.format(threading.get_native_id(), jobs_id[i], error_struct)
                        print(MSG)
                        ensemble.structures[jobs_id[i]].save_scf('t_{}_error_struct_generated_{}.scf'.format(threading.get_native_id(), jobs_id[i]))
                        structures[i].save_scf('t_{}_error_struct_cmp_local_{}.scf'.format(threading.get_native_id(), jobs_id[i]))
                        res["structure"].save_scf('t_{}_error_struct_readed_{}.scf'.format(threading.get_native_id(), jobs_id[i]))

                        continue
                else:
                    print("[WARNING] no check on the structure.")

                is_success =  check_e and check_f
                if get_stress:
                    is_success = is_success and check_s

                if not is_success:
                    continue

                res_only_extra = {x : res[x] for x in res if x not in ["energy", "forces", "stress", "structure"]}
                ensemble.all_properties[num].update(res_only_extra)
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

            self.lock.release()

        # Run until some work has not finished
        recalc = 0
        self.lock = threading.Lock()
        while np.sum(np.array(success, dtype = int) - 1) != 0:
            threads = []

            print("[CYCLE] SUCCESS: ", success)
            print("[CYCLE] STOPPING CONDITION:", np.sum(np.array(success, dtype = int) - 1))

            # Get the remaining jobs
            false_mask = np.array(success) == False
            false_id = np.arange(ensemble.N)[false_mask]

            count = 0
            # Submit in parallel
            jobs = [false_id[i : i + self.job_number] for i in range(0, len(false_id), self.job_number)]
            # Create a local copy of the calculator for each thread, to avoid conflicting modifications
            calculators = [cellconstructor_calc.copy() for i in range(0, len(jobs))]

            for k_th, job in enumerate(jobs):
                # Submit only the batch size
                if count >= self.batch_size:
                    break
                t = threading.Thread(target = compute_single_jobarray, args=(job, calculators[k_th], ))
                t.start()
                threads.append(t)
                count += 1

            # Wait until all the job have finished
            for t in threads:
                t.join(timeout)

            print("[CYCLE] [END] SUCCESS: ", success)
            print("[CYCLE] [END] STOPPING CONDITION:", np.sum(np.array(success, dtype = int) - 1))

            recalc += 1
            if recalc > num_batch_offset + self.max_recalc:
                print ("Expected batch ordinary resubmissions:", num_batch_offset)
                raise ValueError("Error, resubmissions exceeded the maximum number of %d" % self.max_recalc)
                break

        print("CALCULATION ENDED: all properties: {}".format(ensemble.all_properties))



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
        #if self.job_number != 1:
        self.compute_ensemble_batch(ensemble, ase_calc, get_stress, timeout)
        return

