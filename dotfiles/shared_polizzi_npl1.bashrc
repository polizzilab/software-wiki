# Polizzi lab shared .bashrc for workstation
# Copy this file and put it under ~/.bashrc

# This is where COMBS is
export PYTHONPATH="/nfs/polizzi/shared/programs/design/Combs2/"

# Restrict multithreading to 4 threads
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export VECLIB_MAXIMUM_THREADS=4
export NUMEXPR_NUM_THREADS=4
export NUMBA_NUM_THREADS=4
