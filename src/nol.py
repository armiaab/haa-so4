import numpy as np

import haa

DESAIN = {'can-star': ('zeng', 1), 'can-brick': ('zeng', 0), 'so4-star': ('so4', 1), 'so4-brick': ('so4', 0), 'hea': ('hea', 0)}

baris = []
for mol in ('H4', 'LiH', 'HF', 'BeH2'):
    for R in haa.JARAK[mol]:
        for d, (blok, na) in DESAIN.items():
            m = haa.model(mol, R, blok, 2, na)
            g = np.abs(m['G'](np.zeros(len(m['lab']))))
            genap = np.isin(m['lab'], haa.GENAP)
            baris.append(dict(molekul=mol, R=R, desain=d, L=2, grad_maks=g.max(), grad_ganjil_maks=g[~genap].max(),
                              grad_genap_maks=g[genap].max() if genap.any() else 0.0, stasioner=int(g.max() < 1e-5)))
haa.tulis('nol.csv', baris)
print(len(baris))
