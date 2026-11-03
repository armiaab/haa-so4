import numpy as np
from qiskit_algorithms.optimizers import L_BFGS_B

import haa

MOL, R, NA, L = 'H4', 1.5, 1, 2
THETA = np.linspace(-np.pi, np.pi, 121)
M = {nama: haa.model(MOL, R, nama, L, NA, konstanta=False) for nama in ('can', 'so4')}
E0 = haa.dasar(M['can']['H_sys'])
BLOK = {nama: np.flatnonzero(m['lab'] != 'ry') for nama, m in M.items()}


def galat(nama, x):
    return (M[nama]['E'](x) - E0) * 1000


def titik(benih):
    x = np.random.default_rng(benih).uniform(-1, 1, M['can']['lab'].size)
    x[BLOK['can']] = 0
    y = np.zeros(M['so4']['lab'].size)
    y[M['so4']['lab'] == 'ry'] = x[M['can']['lab'] == 'ry']
    return {'can': x, 'so4': y}


potong, gradien, jejak = [], [], []
for nama in M:
    x = titik(3)[nama]
    g = M[nama]['G'](x)
    for k in BLOK[nama]:
        y = x.copy()
        for t in THETA:
            y[k] = t
            potong.append(dict(blok=nama, k=k, pembangkit=M[nama]['lab'][k], gradien=g[k], theta=t, galat_mHa=galat(nama, y)))
for j in range(20):
    for nama, x in titik(100 + j).items():
        gradien += [dict(blok=nama, titik=j, k=k, abs_gradien=abs(g)) for k, g in zip(BLOK[nama], M[nama]['G'](x)[BLOK[nama]])]
for benih in range(12):
    for nama, x in titik(benih).items():
        catat = [galat(nama, x)]
        L_BFGS_B(maxiter=300, maxfun=10 ** 8, callback=lambda x: catat.append(galat(nama, x))).minimize(M[nama]['E'], x,
                                                                                                         jac=M[nama]['G'])
        jejak += [dict(blok=nama, benih=benih, iterasi=i, galat_mHa=e) for i, e in enumerate(catat)]
haa.tulis('intro_potong.csv', potong)
haa.tulis('intro_gradien.csv', gradien)
haa.tulis('intro_jejak.csv', jejak)
print(len(potong), len(gradien), len(jejak))
