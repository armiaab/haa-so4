# HAA-SO(4): code and data

Code and results for the article *Real Two-Qubit Blocks Reduce Gate Costs in the Hardware Adaptable Ansatz*
(A. Abdurohman, in preperation to the Journal of Chemical Theory and Computation).

The study replaces the complex canonical (CAN) two-qubit block of the Hardware Adaptable Ansatz (HAA) of
Zeng et al. (J. Chem. Theory Comput. 2023, 19, 8587) with a real SO(4) block and measures, for H4, LiH, HF, and
BeH2 (STO-3G, frozen core, qubit tapering), how many CNOTs each ansatz needs to reach chemical accuracy at every
bond length, against HAA-CAN and against hardware-efficient ansätze.

## Contents

| Path | Contents |
|---|---|
| `data/` | All results in six gzipped CSV files (see below) |
| `src/haa.py` | Core library: geometries, Hamiltonians with tapering and FCI, ansatz builders, compiled cost, energy and gradient model |
| `src/pec.py`, `pec.sh`, `tugas_pec.txt` | VQE runs (12 starts) per molecule, geometry, ansatz, and number of layers (`pec` = potential energy curve) |
| `src/uji.py`, `uji.sh`, `tugas_pelengkap.txt` | Runs for the pre-registered decision rules (U3 HEA, RealAmplitudes, controls, extended ranges, start tests); `python uji.py tabel` summarizes them |
| `src/ujung.py`, `ujung_cluster.py`, `*.sh` | Exact registered endpoints (LiH 3.9875 Å, BeH2 3.315 Å) and the geometry checks of rule 7 |
| `src/analisis.py` | Summary tables from the VQE runs |
| `src/beh2.py`, `h4.py`, `nol.py`, `m0.py`, `derau.py`, `rute.py`, `sistem.py`, ... | BeH2 one-factor test, H4 deep dive, gradients at zero, M0 floor, noise simulations, routing, molecular descriptors |
| `src/angka.py`, `paper_data.py` | Every number quoted in the article, from the results (pure Python) |
| `src/gambar_hasil.py`, `gambar_skema.py` | Result plots and schematics (Matplotlib) |
| `src/desain_teori.py`, `gambar_photocraft.py` | Theory figures; need the PhotoCraft and VectorCraft command-line tools |
| `src/kumpulkan_hasil.py`, `pulihkan.py` | Collect `results/` into `data/`, and rebuild `results/` from `data/` |

## Quick start

```bash
python src/pulihkan.py               # rebuilds results/ (26,401 CSV files) from data/
cd src && python angka.py angka.tex  # every number of the article, standard library only
```

`pulihkan.py` reproduces the original result files byte for byte, and `angka.py` then reproduces the article's
numbers exactly.

## Rerunning the calculations

The calculations used Python 3.13.5 with the packages pinned in `requirements.txt` (Qiskit 2.5.2,
Qiskit Algorithms 0.4.0, Qiskit Nature 0.8.0, PySCF 2.14.0, SciPy 1.18.1, NumPy 2.5.3). The `*.sh` files are
SLURM array jobs; each task index selects one line of the task list. Results are bit-for-bit reproducible only
when the thread counts are pinned to one:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=1 QISKIT_NUM_THREADS=1
```

## Data

| File | Rows | Contents |
|---|---|---|
| `vqe_runs.csv.gz` | 42,747 | Every VQE run, one row per start: molecule, bond length `R` (Å), block, ancillas `na` (0 = brick, 1 = star), layers `L`, CNOTs, iterations, energies (hartree), error |
| `other_results.csv.gz` | 78,038 | Summary tables, noise simulations, gradients, starting-point tests, geometry checks, molecular descriptors |
| `optimization_histories.csv.gz` | 1,418,926 | Error at every function evaluation of the recorded optimizations |
| `h4_traces.csv.gz` | 224,069 | H4 optimization traces (gradients, even-Y angles) |
| `h4_energy_slices.csv.gz` | 65,340 | H4 energy along single-parameter slices |
| `manifest.csv.gz` | 26,401 | Original file of every row and its header |

Every file starts with a `source` column, the path of the original result file. Column names are Indonesian:
galat = error, lapis = layer, blok = block (`zeng` = published HAA-CAN, `can` = CAN in R_Y layers, `so4`, `real`,
`hea` = R_Y R_Z HEA, `heau3` = U3 HEA, `ra` = RealAmplitudes), sudut = angle, jarak = range, setimbang =
equilibrium, regang = stretched, derau = noise, diam = idle, awal = start, acak = random, nol = zero, lantai =
floor, hemat = saving, biaya = cost, kedalaman = depth, rute = routing, uji = test, kontrol = control, terbaik =
best, sukses = success, impas = break-even, peta = qubit mapping, riwayat = history, jejak = trace, potong = slice,
kena_batas = iteration limit reached, buang = removed.
