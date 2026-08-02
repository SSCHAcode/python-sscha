# -*- coding: utf-8 -*-
"""Execute a pickled ASE calculator on a JSON-serialized structure.

Invoked by sscha.ClusterCalculators.ASEFileCalculator on the (local or
remote) cluster as:

    python -m sscha.ASEClusterRunner --calc calculator.pkl \
        --input PREFIX_input.json --output PREFIX.json

The output JSON contains energy [eV], forces [eV/Angstrom], the computed
Atoms and, when the calculator supports it, stress [ASE Voigt
(xx, yy, zz, yz, xz, xy), eV/Angstrom^3].
"""
import argparse
import pickle

import ase.io.jsonio


def main():
    """Run the calculation and dump the results."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--calc", required=True,
                        help="Path to the pickled ASE calculator.")
    parser.add_argument("--input", required=True,
                        help="Path to the JSON-serialized ASE Atoms input.")
    parser.add_argument("--output", required=True,
                        help="Path of the JSON file where results are written.")
    args = parser.parse_args()

    with open(args.calc, "rb") as handle:
        calc = pickle.load(handle)
    with open(args.input, "r") as handle:
        atoms = ase.io.jsonio.decode(handle.read())

    atoms.calc = calc
    payload = {
        "energy": float(atoms.get_potential_energy()),  # eV
        "forces": atoms.get_forces(),                   # eV / Angstrom
        "atoms": atoms,                                 # structure actually computed
    }
    try:
        payload["stress"] = atoms.get_stress()          # ASE Voigt, -eV/Angstrom^3
    except Exception:
        pass  # calculator without stress: call compute_ensemble with get_stress=False
    with open(args.output, "w") as handle:
        handle.write(ase.io.jsonio.encode(payload))


if __name__ == "__main__":
    main()
