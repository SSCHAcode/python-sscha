#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate base
export OMP_NUM_THREADS=4

python run_aiida_flare_sscha.py > log2

python run_flare_sscha.py > log3
