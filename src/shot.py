import os
import time
from multiprocessing import Pool

import numpy as np
from qiskit.primitives import StatevectorEstimator
from qiskit_algorithms.optimizers import SPSA

import haa

SIGMA = 1e-3
RESTART = 12
MAXITER = 3000
TUGAS = [(b['molekul'], max(haa.JARAK[b['molekul']]), b['blok'], int(b['L']), int(b['topologi'] == 'star'))
         for b in haa.baca('tabel/tab4_pembanding.csv')]


def jalan(tugas):
    mol, R, ansatz, L, na = tugas
    m = haa.model(mol, R, ansatz, L, na)
    P = len(m['lab'])
    baris = []
    for s in range(RESTART):
        bising = StatevectorEstimator(default_precision=SIGMA, seed=np.random.default_rng(1000 + s))

        def E(x):
            return float(bising.run([(m['qc'], m['H'], x)]).result()[0].data.evs) + m['offset']

        t0 = time.perf_counter()
        hasil = SPSA(maxiter=MAXITER).minimize(E, np.random.default_rng(s).uniform(-1, 1, P))
        e = m['E'](hasil.x)
        baris.append(dict(molekul=mol, R=R, blok=ansatz, L=L, na=na, parameter=P, sigma=SIGMA, restart=s,
                          maxiter=MAXITER, nfev=hasil.nfev, detik=round(time.perf_counter() - t0, 1), E_eksak=e,
                          E_fci=m['e_fci'], galat_mHa=(e - m['e_fci']) * 1000))
    return baris


if __name__ == '__main__':
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))) as pool:
        haa.tulis('shot.csv', sum(pool.map(jalan, TUGAS), []))
    print(len(TUGAS))
