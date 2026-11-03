#!/bin/bash
#SBATCH --job-name=haa
#SBATCH --partition=qdisk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=%x-%j.out
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=1 QISKIT_NUM_THREADS=1
export PYTHONPATH=$SLURM_SUBMIT_DIR/../.pylib${PYTHONPATH:+:$PYTHONPATH}
PY=${PY:-python}
$PY -c "import numpy, scipy, pyscf, qiskit, qiskit_nature, qiskit_algorithms as a; print(numpy.__version__, scipy.__version__, pyscf.__version__, qiskit.__version__, qiskit_nature.__version__, a.__version__)" || exit 1
$PY "$@"
