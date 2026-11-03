#!/usr/bin/env bash
# Run ujung.py on a local machine (WSL/NixOS) instead of the cluster.
#
# Environment: the research flake's Python (~/research-env, nix develop) plus the paper's exact package versions
# installed once into ~/haa-pylib:
#   cd ~/research-env && nix develop -c python -m pip install --target ~/haa-pylib \
#       numpy==2.5.3 scipy==1.18.1 qiskit==2.5.2 qiskit-algorithms==0.4.0 qiskit-nature==0.8.0 pyscf==2.14.0
# This reproduced results/uji/jarak_LiH_so4_na0_R4.000_L2_r0.csv bit for bit (2026-11-02).
#
# Usage (from Windows):  wsl -- bash -lc 'cd ~/research-env && nix develop -c bash "<path to src>/ujung_wsl.sh" 8'
# Finished tasks are skipped, so the script can be restarted after an interruption.
set -euo pipefail
cd "$(dirname "$0")"
WORKERS=${1:-8}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=1 QISKIT_NUM_THREADS=1
export PYTHONPATH=$HOME/haa-pylib${PYTHONPATH:+:$PYTHONPATH}
export PYTHONWARNINGS=ignore
N=$(python ujung.py)
echo "$(date '+%F %T') ujung: $N tasks, $WORKERS workers"
seq 1 "$N" | xargs -P "$WORKERS" -I{} sh -c 'python ujung.py {} >/dev/null 2>>"$HOME/ujung_errors.log" || echo "task {} failed" >>"$HOME/ujung_errors.log"'
echo "$(date '+%F %T') ujung: done"
python ujung.py cek
python uji.py tabel >/dev/null
echo "$(date '+%F %T') summaries written (results/tabel/ujung_cek.csv, results/tabel/uji_lstar.csv)"
