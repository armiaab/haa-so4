import os
import time
from multiprocessing import Pool

import numpy as np
from qiskit_algorithms.optimizers import L_BFGS_B

import haa

RESTART = 12
TUGAS = [(b['molekul'], R, b['blok'], int(b['L']), int(b['na'])) for b in haa.baca('tabel/tab2_biaya.csv')
         if b['molekul'] in ('LiH', 'HF') for R in haa.JARAK[b['molekul']]]


def jalan(tugas):
    mol, R, ansatz, L, na = tugas
    m = haa.model(mol, R, ansatz, L, na)
    P = len(m['lab'])
    riwayat, ringkas = [], []
    for s in range(RESTART):
        catat = []

        def E(x):
            catat.append(m['E'](x))
            return catat[-1]

        t0 = time.perf_counter()
        hasil = L_BFGS_B(maxiter=3000, maxfun=10 ** 8).minimize(E, np.random.default_rng(s).uniform(-1, 1, P), jac=m['G'])
        kunci = dict(molekul=mol, R=R, blok=ansatz, na=na, L=L, restart=s)
        ringkas.append(dict(kunci, ns=m['ns'], P=P, cnot=haa.biaya(ansatz, m['ns'], na, L)[0], e_hf=m['e_hf'],
                            e_fci=m['e_fci'], galat_awal=(catat[0] - m['e_fci']) * 1000,
                            galat_akhir=(hasil.fun - m['e_fci']) * 1000, nit=hasil.nit, nfev=len(catat),
                            detik=round(time.perf_counter() - t0, 1)))
        galat = (np.minimum.accumulate(catat) - m['e_fci']) * 1000
        riwayat += [dict(kunci, evaluasi=i + 1, galat_mHa=round(g, 7)) for i, g in enumerate(galat)]
    return riwayat, ringkas


if __name__ == '__main__':
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))) as pool:
        hasil = pool.map(jalan, TUGAS)
    haa.tulis('riwayat.csv', sum((h[0] for h in hasil), []))
    haa.tulis('riwayat_ringkas.csv', sum((h[1] for h in hasil), []))
    print(len(TUGAS))
