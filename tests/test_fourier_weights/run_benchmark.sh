#!/bin/bash

NCFGS=(16 32 64 128 256 512 1024 2048)
CELLS=(1 2 3 4 5 6)

for CELL in ${CELLS[@]}; do
	cell_dir="benchmark/${CELL}x${CELL}x${CELL}"
	mkdir -p $cell_dir
	sed -i "s/supercell = .*/supercell = (${CELL}, ${CELL}, ${CELL})/g" benchmark.py
	for NCFG in ${NCFGS[@]}; do
		echo "Running benchmark with NCFG=${NCFG} and CELL=${CELL}"
		sed -i "s/n_configs = .*/n_configs = ${NCFG}/g" benchmark.py
		python benchmark.py > $cell_dir/N_${NCFG}.out 
	done
done
