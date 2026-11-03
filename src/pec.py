import sys
import time

import numpy as np
from qiskit_algorithms.optimizers import L_BFGS_B

import haa

MOL, R, ANSATZ, L, NA = sys.argv[1], float(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
PETA = sys.argv[6] if len(sys.argv) > 6 else 'taper'
RESTART = 12
MAXITER = 3000

m = haa.model(MOL, R, ANSATZ, L, NA, PETA)
P = len(m['lab'])
cnot, dalam = haa.biaya(ANSATZ, m['ns'], NA, L)
e_dasar = haa.dasar(m['H_sys']) + m['offset']

baris = []
for s in range(RESTART):
    t0 = time.perf_counter()
    x0 = np.random.default_rng(s).uniform(-1, 1, P)
    hasil = L_BFGS_B(maxiter=MAXITER, maxfun=10 ** 8).minimize(m['E'], x0, jac=m['G'])
    baris.append(dict(molekul=MOL, R=R, blok=ANSATZ, L=L, peta=PETA, ns=m['ns'], na=NA, parameter=P, cnot=cnot,
                      kedalaman=dalam, restart=s, nit=hasil.nit, nfev=hasil.nfev, kena_batas=int(hasil.nit >= MAXITER),
                      detik=round(time.perf_counter() - t0, 1), E_vqe=hasil.fun, E_hf=m['e_hf'], E_fci=m['e_fci'],
                      E_dasar=e_dasar, galat_mHa=(hasil.fun - m['e_fci']) * 1000))

haa.tulis(f'pec/{MOL}_{PETA}_na{NA}_R{R:.3f}_{ANSATZ}_L{L}.csv', baris)
print(min(b['galat_mHa'] for b in baris))
