import os
import time
from multiprocessing import Pool

import numpy as np
from qiskit_algorithms.optimizers import L_BFGS_B

import haa

MOL, R = 'H4', 1.5
TUGAS = [(nama, na, L, s) for nama in ('can', 'zeng', 'real', 'so4') for na in (1, 0) for L in (2, 3)
         for s in range(4 if (nama, na) == ('can', 1) else 2)]


def jalan(tugas):
    nama, na, L, s = tugas
    m = haa.model(MOL, R, nama, L, na, konstanta=False)
    P = len(m['lab'])
    t0 = time.perf_counter()
    hasil = L_BFGS_B(maxiter=3000, maxfun=10 ** 8).minimize(m['E'], np.random.default_rng(s).uniform(-1, 1, P), jac=m['G'])
    e0 = haa.dasar(m['H_sys'])
    baris = dict(molekul=MOL, R=R, blok=nama, na=na, L=L, parameter=P, dim_bola_real=2 ** (m['ns'] + na) - 1, restart=s,
                 galat_mHa=(hasil.fun - e0) * 1000, **haa.arah(m['qc'], hasil.x), nit=hasil.nit,
                 detik=round(time.perf_counter() - t0, 1))
    return baris, [dict(molekul=MOL, R=R, blok=nama, na=na, L=L, restart=s, k=k, theta=t) for k, t in enumerate(hasil.x)]


if __name__ == '__main__':
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))) as pool:
        hasil = pool.map(jalan, TUGAS)
    haa.tulis('fasa.csv', [h[0] for h in hasil])
    haa.tulis('fasa_sudut.csv', sum((h[1] for h in hasil), []))
    print(len(TUGAS))
