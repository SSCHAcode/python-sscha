import cellconstructor as CC
import sscha, sscha.Cluster
import threading, copy

import numpy as np
import scipy, scipy.interpolate

import sys, os

MODULE_DESCRIPTION = '''
This module contains useful functions to do post-processing analysis of the SSCHA.

For example, it implements custom cluster calculators to compute the electron-phonon effect on
absorbption and the bandgap.

'''


class OpticalQECluster(sscha.Cluster.Cluster):
    '''
    This class does the same as the sscha.Cluster to submit the calculation of an ensemble
    but it also computes the spectral properties.
    '''

    def __init__(self, new_k_grid = None, random_offset = True, epsilon_data = None,
                epsilon_binary = 'epsilon.x -npool NPOOL -i PREFIX.pwi > PREFIX.pwo', 
                **kwargs):
        '''
        Initialize the cluster object.

        Parameters
        ----------

            new_k_grid : list
                The dimension of the new k mesh to run the nscf calculation
            random_offset : bool
                If True, a random offset is added to the calculation
            epsilon_data : dict
                The dictionary with the information of the namespace for the
                epsilon.x file
            epsilon_binary : string
                The path to the epsilon.x binary inside the cluster.
            **kwargs : 
                All other arguments to be passed to the cluster.

        '''
        self.new_k_grid = new_k_grid

        if random_offset:
            self.random_offset = np.random.uniform(size = 3)
        else:
            self.random_offset = np.zeros(3)

        self.kpts = None
        self.epsilon_binary = epsilon_binary


        self.read_sigma = False

        self.epsilon_data = {
            'inputpp' : {
                'calculation' : 'sigma',
                'prefix' : 'AHAH'
            },
            'energy_grid' : {
                'wmin' : 0,
                'wmax' : 40,
                'nw' : 10000,
                'temperature' : 300
            }   
        }
        if epsilon_data is not None:
            self.epsilon_data = epsilon_data

        super().__init__(**kwargs)


    def __setattr__(self, __name, __value):
        super().__setattr__(__name, __value)

        # Always regenerate the kpts
        if __name == 'new_k_grid' and not __value is None: 
            assert len(__value) == 3, 'Error, new_k_grid must be a tuple with 3 elements'
            self.generate_kpts()
    

    def generate_kpts(self):
        '''
        Generate the kpts for the non self-consistent calculation
        '''

        nkpts = int(np.prod(self.new_k_grid))

        self.kpts = np.zeros((nkpts, 3), dtype = np.double)

        for i in range(nkpts):
            x = float(i // (self.new_k_grid[1] * self.new_k_grid[2]))
            res = i % (self.new_k_grid[1] * self.new_k_grid[2])
            y = float(res // self.new_k_grid[2])
            z = float(res % self.new_k_grid[2])

            x /= self.new_k_grid[0]
            y /= self.new_k_grid[1]
            z /= self.new_k_grid[2]


            x = (x + self.random_offset[0]) % 1
            y = (y + self.random_offset[1]) % 1
            z = (z + self.random_offset[2]) % 1

            self.kpts[i, :] = [x, y, z]


    def get_execution_command(self, label):
        '''
        Get the command to execute the two pw.x and epsilon.x
        calculations
        '''

        # Get the MPI command replacing NPROC
        new_mpicmd  = self.mpi_cmd.replace("NPROC", str(self.n_cpu))
        
        # Replace the NPOOL variable and the PREFIX in the binary
        binary = self.binary.replace("NPOOL", str(self.n_pool)).replace("PREFIX", label)
        binary_nscf = self.binary.replace("NPOOL", str(self.n_pool)).replace("PREFIX", label + '_nscf')
        binary_eps = self.epsilon_binary.replace("NPOOL", str(self.n_pool)).replace("PREFIX", label + '_eps')

        tmt_str = ""
        if self.use_timeout:
            tmt_str = "timeout %d " % self.timeout

        submission_txt = '''
{0} {1} {2}
{0} {1} {3}
{0} {1} {4}
'''.format(tmt_str, new_mpicmd, binary, binary_nscf, binary_eps)

        # Remove the wavefunction and other heavy data
        submission_txt += '''
rm -rf {0}.wfc* {0}.save

'''.format(label)

        return submission_txt


    def read_results(self, calc, label):
        '''
        Get the results
        '''

        results = super().read_results(calc, label)

        print('THREAD ID {} START READING THE EPSILON'.format(threading.get_native_id()))

        # Add the additional information related to the  epsilon
        prefix = label
        epsilon_data = None
        eps_real = np.loadtxt(os.path.join(self.local_workdir, 'epsr_{}.dat'.format(prefix)))
        eps_imag = np.loadtxt(os.path.join(self.local_workdir, 'epsi_{}.dat'.format(prefix)))

        nw = eps_real.shape[0]
        epsilon_data = np.zeros((nw, 3), dtype = np.double)

        epsilon_data[:,0] = eps_real[:,0]
        epsilon_data[:, 1] = np.mean(eps_real[:, 1:], axis = 1)
        epsilon_data[:, 2] = np.mean(eps_imag[:, 2:], axis = 1)

        results['epsilon'] = epsilon_data
        print('THREAD ID {} READED ALSO THE EPSILON!'.format(threading.get_native_id()))

        return results



    def prepare_input_file(self, structures, calc, labels):
        '''
        Prepare the input files for the cluster
        '''    


        # Prepare the input file
        list_of_inputs = []
        list_of_outputs = []
        for i, (label, structure) in enumerate(zip(labels, structures)):
            # Avoid thread conflict
            self.lock.acquire()
            
            try:
                calc.set_directory(self.local_workdir)
                PREFIX = label
                old_input = copy.deepcopy(calc.input_data)
                old_kpts = copy.deepcopy(calc.kpts)


                calc.input_data['control'].update({'prefix' : PREFIX})
                calc.set_label(label)
                calc.write_input(structure)

                print("[THREAD {}] LBL: {} | PREFIX: {}".format(threading.get_native_id(), label, calc.input_data["control"]["prefix"]))


                # prepare the first scf calculation
                input_file = '{}.pwi'.format(label)
                output_file = '{}.pwo'.format(label)

                list_of_inputs.append(input_file)     
                list_of_outputs.append(output_file)

                # prepare the nscf calculation
                calc.input_data['control'].update({'calculation' : 'nscf'})
                calc.input_data['system'].update({'nosym' : True})
                calc.kpts = self.kpts.copy()

                # Generate the input file
                new_label = '{}_nscf'.format(label)
                calc.set_label(new_label, override_prefix = False)
                calc.write_input(structure)
                input_file = '{}.pwi'.format(new_label)
                output_file = '{}.pwo'.format(new_label)
                list_of_inputs.append(input_file)     
                list_of_outputs.append(output_file)


                # Prepare the epsilon calculation
                eps_namelist = copy.deepcopy(self.epsilon_data)
                eps_namelist['inputpp'].update({'prefix' : calc.input_data['control']['prefix']})
                eps_lines = CC.Methods.write_namelist(eps_namelist)


                new_label = '{}_eps'.format(label)
                input_file = '{}.pwi'.format(new_label)
                output_file = '{}.pwo'.format(new_label)

                # Write the epsilon input file
                eps_in_filename = os.path.join(self.local_workdir, input_file)
                with open(eps_in_filename, 'w') as fp:
                    fp.writelines(eps_lines)
                    print('fname: {} prepared'.format(eps_in_filename))


                
                list_of_inputs.append(input_file)     
                list_of_outputs.append(output_file)

                # Append also the imaginary and real part of epsilon and sigma
                outputs = ['epsi_{}.dat', 'epsr_{}.dat', 'sigmai_{}.dat', 'sigmar_{}.dat']
                list_of_outputs += [x.format(PREFIX) for x in outputs]


                # Reset the calculator with the old data
                calc.input_data = old_input
                calc.kpts = old_kpts



            except Exception as e:
                MSG = '''
Error while writing input file {}.
'''.format(label)
                print(MSG)
                print('ERROR MSG:')
                print(e)

            # Release the lock on the threads
            self.lock.release()            

        print('THREAD: {} inputs: {} outputs: {}'.format(threading.get_native_id(), list_of_inputs, list_of_outputs))
            
        
        return list_of_inputs, list_of_outputs



def get_optical_spectrum(ensemble, w_array = None):
    """
    COMPUTE THE OPTICAL SPECTRUM
    ============================

    By averaging the dielectric properties of phonon-displaced configurations,
    we can get the phonon-renormalized optical spectrum.


    Parameters
    ----------
        ensemble : sscha.Ensemble.Ensemble
            The ensemble. It should contain the optical data.
            To do it, use the cluster calculator OpticalQECluster of
            the AdvancedCalculation module.
        w_array : ndarray, Optional
            The frequencies at which to interpolate the results.
            If not passed, the full dataset of frequencies will be returned.

    Results
    -------
        w : ndarray, dtype = float
            The frequency (eV)
        n : ndarray, dtype = complex
            The (complex) refractive index
    """

    w_data = None
    for i in range(ensemble.N):
        if not 'epsilon' in ensemble.all_properties[i]:
            ERR = """
Error, the configuration {} has no 'epsilon' data.
""".format(i)
            raise ValueError(ERR)
        
        data = np.array(ensemble.all_properties[i]['epsilon'])

        if w_data is None:
            w_data = data[:,0]
            eps_real = np.zeros_like(w_data)
            eps_imag = np.zeros_like(w_data)


        eps_real += data[:,1]
        eps_imag += data[:,2]

    eps_real /= ensemble.N
    eps_imag /= ensemble.N

    # Interpolate the data if requested
    if w_array is not None:
        f_real = scipy.interpolate.interp1d(w_data, eps_real,
            kind = 'cubic', bounds_error = False, fill_value = 'extrapolate')
        eps_real = f_real(w_array)
        f_imag = scipy.interpolate.interp1d(w_data, eps_imag,
            kind = 'cubic', bounds_error = False, fill_value = 'extrapolate')
        eps_imag = f_imag(w_array)
        w_data = w_array
    
    # Build the complex epsilon
    epsilon = eps_real + 1j * eps_imag

    # Build the refractive index
    n = np.sqrt(epsilon)

    return w_data, n


