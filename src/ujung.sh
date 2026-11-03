#!/bin/bash
#SBATCH --job-name=haa-ujung
#SBATCH --partition=qdisk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --array=1-1000
#SBATCH --output=%x-%A_%a.out
# Exact rule-2 endpoints and rule-7 checks (ujung.py). Same environment as uji.sh; finished tasks are skipped,
# so the job can be resubmitted after a timeout.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=1 QISKIT_NUM_THREADS=1
export PYTHONPATH=$SLURM_SUBMIT_DIR/../.pylib${PYTHONPATH:+:$PYTHONPATH}
PY=${PY:-python}
$PY -c "import numpy, scipy, pyscf, qiskit, qiskit_nature, qiskit_algorithms as a; print(numpy.__version__, scipy.__version__, pyscf.__version__, qiskit.__version__, qiskit_nature.__version__, a.__version__)" || exit 1
for i in $(seq $SLURM_ARRAY_TASK_ID 1000 $($PY ujung.py)); do
  $PY ujung.py $i
done
