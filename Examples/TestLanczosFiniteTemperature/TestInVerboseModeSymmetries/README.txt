Test for checking the correct symmetrization.
The test_apply_L_FT.py generates with and without symmetrization
in the Lanczos the output of output_sym and output_nosym, toghether with the two
d3 in the polarization base (symmetrized and not).

Then the isolate_subspace.sh script get the d3 from the output and its symmetrized counterpart,
as computed in the C function and printed in the DEB = 1 flag.

This is performed on a degenerate subspace given by the mode 17,31,41, 
that correspond to the maximum of the absolute value of the d3. 
The result is written into analisys_sym.dat, and should be compared with the d3_modes_sym.npy.

If they match, then the symmetrization is performed correctly inside the C function.
The application of the full Lanczos operator is compared at the end of output_sym and
output_nosym. If it is zero, then the C function is working properly.

This is the detailed test to benchmark the C implementation of the Lanczos

The check_symmetries.py performes a check on the symmetry matrix as they are read from
input, to assure that the symmetrization procedure is occurring correclty.
The otuptu dynamical matrices must be compared with the symmetries as parsed by the Lanczos.

To produce the cleaned output that is parsed faster by the check_symmetries.py script one can use the grep method as follows:

>>> grep 'Mode\|Deg space\|IN_VEC_OUT_DYN' output_sym  > cleaned_output_sym 

Where output_sym is the output of the test_apply_L_FT.py.
