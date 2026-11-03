#!/bin/bash
#SBATCH --job-name=haa-ujung
#SBATCH --partition=qdisk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --array=1-1000
#SBATCH --output=%x-%A_%a.out
# Exact rule-2 endpoints and rule-7 checks (ujung_cluster.py, 1358 tasks).
# Put ujung_cluster.py and this file in one folder and run:  sbatch ujung_cluster.sh
# Array task j runs tasks j, j+1000, ... up to 1358. Finished tasks are skipped, so the job can be resubmitted.
# Output goes to $OUT (default: ./hasil_ujung). Packages come from $PYLIB if set (default ../.pylib, as in uji.sh).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=1 QISKIT_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
PYLIB=${PYLIB:-$SLURM_SUBMIT_DIR/../.pylib}
[ -d "$PYLIB" ] && export PYTHONPATH=$PYLIB${PYTHONPATH:+:$PYTHONPATH}
PY=${PY:-python}
$PY -c "import numpy, scipy, pyscf, qiskit, qiskit_nature, qiskit_algorithms as a; print(numpy.__version__, scipy.__version__, pyscf.__version__, qiskit.__version__, qiskit_nature.__version__, a.__version__)" || exit 1
N=$($PY ujung_cluster.py)
for i in $(seq $SLURM_ARRAY_TASK_ID 1000 $N); do
  $PY ujung_cluster.py $i
done
