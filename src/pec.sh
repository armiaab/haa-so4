#!/bin/bash
#SBATCH --job-name=haa-pec
#SBATCH --partition=qdisk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --array=1-1000
#SBATCH --output=%x-%A_%a.out
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=1 QISKIT_NUM_THREADS=1
export PYTHONPATH=$SLURM_SUBMIT_DIR/../.pylib${PYTHONPATH:+:$PYTHONPATH}
PY=${PY:-python}
$PY -c "import numpy, scipy, pyscf, qiskit, qiskit_nature, qiskit_algorithms as a; print(numpy.__version__, scipy.__version__, pyscf.__version__, qiskit.__version__, qiskit_nature.__version__, a.__version__)" || exit 1
for i in $(seq $SLURM_ARRAY_TASK_ID 1000 $(wc -l < tugas_pec.txt)); do
  $PY pec.py $(sed -n "${i}p" tugas_pec.txt)
done
