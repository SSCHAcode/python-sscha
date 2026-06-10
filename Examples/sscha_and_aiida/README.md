# Instructions

We provide here some scripts to run SSCHA using AiiDA, FLARE machine-learning potential, and a combination of the two as on-the-fly active learning.

* `run_aiida_sscha.py`, which performs a thermal expansion calculation using SSCHA and aiida-quantumespresso. 
* `run_flare_sscha.py`, which performs a thermal expansion calculation using SSCHA and FLARE machine-learning interatomic potential. 
* `run_aiida_flare_sscha.py`, which performs an on-the-fly active learning SSCHA calculation using aiida-quantumespresso for DFT and FLARE as the ML potential. 

It is preferable to execute the example that you already have some experience with both the SSCHA and AiiDA-QuantumESPRESSO codes. Nevertheless, you can try following the instructions and to run the example.

## Installation

We recommend to install all the packages via `mamba` (which is based on `conda`). After having installed `mamba`, you can simply run the following:

```console
> mamba create -n sscha-aiida -c conda-forge python gfortran "blas=*=openblas" openblas lapack julia pip numpy scipy spglib pkg-config aiida-core
> pip install ase cellconstructor python-sscha aiida-quantumespresso aiida-pseudo
> mamba activate sscha-aiida
```
Then, you should configure an AiiDA profile and the pw.x code in order to use the example script (see also Prerequisites section).

For on-the-fly active learning you also need the FLARE package:

```console
> git clone --depth 1 https://github.com/mir-group/flare.git
> cd flare
> pip install .
```

## Prerequisites

You need to have installed in the same environment:
- `python-sscha`
- `cellconstructor`
- `aiida-core`
- `aiida-quantumespresso`
- (optional, for active learning) `mir-flare` 

For the AiiDA part, it is essential the dameon is running and you have:
1. Configured a computer where to run the code (e.g. on your own laptop; see `aiida-core` docs for further details)
2. Configured a code for the `pw.x` binary, related to the computer (see `aiida-quantumespresso` docs for further details)
3. Installed a pseupopotentials library via `aiida-pseudo`. E.g. `aiida-pseudo install -x PBEsol -v 1.3`

Please refer to the aiida-core and aiida-quantumespresso for further installation instruction.

## How-to run 

### SSCHA with the aiida-quantumespresso interface 

Open the `run_aiida_sscha.py` and change the data according to your needs and local installation. Then simply

```console
> python run_aiida_sscha.py > run_aiida_sscha.log
```

Usually an actual production run would take a while. We suggest to use instead

```console
> nohup python run_aiida_sscha.py > run_aiida_sscha.log &
```

or a submit script at glance.

### On-the-fly active learning SSCHA with aiida-quantumespresso and FLARE interface 

Open the `run_aiida_flare_sscha.py` and change the data according to your needs and local installation. Then simply

```console
> python run_aiida_flare_sscha.py > run_aiida_sscha.log
```

Usually an actual production run would take a while. We suggest to use instead

```console
> nohup python run_aiida_sscha.py > run_aiida_sscha.log &
```

or a submit script at glance.