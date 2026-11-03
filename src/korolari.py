import numpy as np
from qiskit_algorithms.optimizers import L_BFGS_B

import haa

MOL, R = 'H4', 1.5

baris = []
for nama, na, L in [('can', 1, 2), ('can', 0, 3), ('zeng', 1, 2), ('so4', 1, 2)]:
    m = haa.model(MOL, R, nama, L, na, konstanta=False)
    e0 = haa.dasar(m['H_sys'])
    genap = np.isin(m['lab'], haa.GENAP)
    x = np.random.default_rng(3).uniform(-1, 1, len(m['lab']))
    x[genap] = 0
    g = np.abs(m['G'](x))
    hasil = L_BFGS_B(maxiter=300, maxfun=10 ** 8).minimize(m['E'], x, jac=m['G'])
    baris.append(dict(blok=nama, na=na, L=L, grad_genap=g[genap].max() if genap.any() else '', grad_ganjil=g[~genap].max(),
                      sudut_genap_akhir=np.abs(hasil.x[genap]).max() if genap.any() else '',
                      galat_mHa=(hasil.fun - e0) * 1000, korelasi_mHa=(m['e_hf'] - m['e_fci']) * 1000))
haa.tulis('korolari.csv', baris)
print(baris)
