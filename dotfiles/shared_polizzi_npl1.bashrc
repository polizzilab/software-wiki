# Polizzi lab shared .bashrc for workstation

# By default, python's numerical linear algebra libraries will use
# the number of threads on the entire machine. If multiple users
# run COMBS simultaneously this can easily overwhelm the processers.
# To prevent this, we limit the default number of threads used by
# these libraries.

export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export VECLIB_MAXIMUM_THREADS=4
export NUMEXPR_NUM_THREADS=4
export NUMBA_NUM_THREADS=4
