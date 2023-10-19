# Instructions

We provide here the script `run_aiida_sscha.py`, which performs a thermal expansion calculation using SSCHA and aiida-quantumespresso. 

It is preferable to execute the example that you already have some experience with both the SSCHA and AiiDA-QuantumEspresso codes. Nevertheless, you can try following the instructions and to run the example.

## Installation

We recommend to install all the packages via `mamba` (which is based on top of `conda`). After having installed `mamba`, you can simply run the following:

```console
> mamba create -n aiida-sscha -c conda-forge python gfortran libblas lapack openmpi julia openmpi-mpicc pip numpy scipy spglib aiida-core
> pip install ase quippy-ase cellconstructor python-sscha aiida-quantumespresso aiida-pseudo
> mamba activate aiida-sscha
```

Then, you should configure an AiiDA profile in order to use the example script (see also Prerequisites section).

## Prerequisites

You need to have installed in the same environment:
- `python-sscha`
- `cellconstructor`
- `aiida-core`
- `aiida-quantumespresso`

For the AiiDA part, it is essential the dameon is running and you have:
1. Configured a computer where to run the code (e.g. on your own laptop; see `aiida-core` docs for further details)
2. Configured a code for the `pw.x` binary, related to the computer (see `aiida-quantumespresso` docs for further details)
3. Installed a pseupopotentials library via `aiida-pseudo`. E.g. `aiida-pseudo install -x PBEsol -v 1.3`

Please refer to the aiida-core and aiida-quantumespresso for further installation instruction.

## How-to run 

Open the `run_aiida_sscha.py` and change the data according to your needs and local installation. Then simply

```console
> python run_aiida_sscha.py > run_aiida_sscha.log
```

Usually an actual production run would take a while. We suggest to use instead

```console
> nohup python run_aiida_sscha.py > run_aiida_sscha.log &
```

or a submit script at glance.