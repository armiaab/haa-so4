import glob
import os
import sys
import time
from collections import defaultdict

import numpy as np
from qiskit_algorithms.optimizers import ADAM, L_BFGS_B

import haa

MOL = 'BeH2'
RESTART = 12
MAXITER = {'lbfgsb': 3000, 'adam': 10000}
JAUH = [0.7, 2.5, 3.0, 3.3]
TUGAS = []
for R in JAUH:
    TUGAS += [('jarak', R, 'hea', L, 0, 'taper', 'lbfgsb') for L in range(1, 13)]
    TUGAS += [('jarak', R, b, L, na, 'taper', 'lbfgsb') for b, na in (('so4', 1), ('so4', 0), ('zeng', 1)) for L in range(1, 9)]
for R in haa.JARAK[MOL]:
    TUGAS += [('adam', R, 'hea', L, 0, 'taper', 'adam') for L in range(4, 11)]
    TUGAS += [('u3', R, 'heau3', L, 0, 'taper', 'lbfgsb') for L in range(1, 13)]
for R in (1.0, 1.326, 2.0):
    TUGAS += [('peta', R, 'hea', L, 0, 'jw', 'lbfgsb') for L in (2, 4, 7, 10, 15, 20, 25)]
TUGAS = [t + (k,) for t in TUGAS for k in range(RESTART)]


def jalan(i):
    faktor, R, blok, L, na, peta, opt, k = TUGAS[i - 1]
    nama = f'beh2/{faktor}_{blok}_na{na}_{peta}_R{R:.3f}_L{L}_r{k}.csv'
    if os.path.exists(os.path.join(haa.HASIL, nama)):
        return
    m = haa.model(MOL, R, blok, L, na, peta)
    P = len(m['lab'])
    cnot, dalam = haa.biaya(blok, m['ns'], na, L)
    x0 = np.random.default_rng(k).uniform(-1, 1, P)
    optimizer = L_BFGS_B(maxiter=MAXITER[opt], maxfun=10 ** 8) if opt == 'lbfgsb' else ADAM(maxiter=MAXITER[opt])
    t0 = time.perf_counter()
    hasil = optimizer.minimize(m['E'], x0, jac=m['G'])
    detik = round(time.perf_counter() - t0, 1)
    e = m['E'](hasil.x)
    nit = hasil.nit or hasil.nfev
    haa.tulis(nama, [dict(faktor=faktor, molekul=MOL, R=R, blok=blok, L=L, na=na, peta=peta, opt=opt, ns=m['ns'], parameter=P,
                          cnot=cnot, kedalaman=dalam, restart=k, nit=nit, nfev=hasil.nfev,
                          kena_batas=int(nit >= MAXITER[opt]), detik=detik, E_vqe=e, E_hf=m['e_hf'],
                          E_fci=m['e_fci'], E_dasar=haa.dasar(m['H_sys']) + m['offset'], galat_mHa=(e - m['e_fci']) * 1000)])
    print(nama, (e - m['e_fci']) * 1000)


def tabel():
    data = defaultdict(list)
    for f in glob.glob(os.path.join(haa.HASIL, 'pec', f'{MOL}_taper_*.csv')):
        for b in haa.baca(f):
            if (b['blok'], b['na']) in (('hea', '0'), ('so4', '0'), ('so4', '1'), ('zeng', '1')):
                data[('dasar', b['blok'], b['na'], 'taper', 'lbfgsb', float(b['R']), int(b['L']))].append(b)
    for f in glob.glob(os.path.join(haa.HASIL, 'beh2', '*.csv')):
        for b in haa.baca(f):
            data[(b['faktor'], b['blok'], b['na'], b['peta'], b['opt'], float(b['R']), int(b['L']))].append(b)
    titik = []
    for (faktor, blok, na, peta, opt, R, L), bs in sorted(data.items()):
        g = np.array([float(b['galat_mHa']) for b in bs])
        titik.append(dict(faktor=faktor, blok=blok, na=na, peta=peta, opt=opt, R=R, L=L, n=len(g), cnot=bs[0]['cnot'],
                          terbaik_mHa=g.min(), sukses=int((g < 1e3 * haa.CHEM).sum()),
                          kena_batas=sum(int(b['kena_batas']) for b in bs),
                          dasar_minus_fci_mHa=1e3 * (float(bs[0]['E_dasar']) - float(bs[0]['E_fci']))))
    haa.tulis('tabel/beh2.csv', titik)
    terbaik = {(t['faktor'], t['blok'], t['na'], t['peta'], t['opt'], t['R'], t['L']): t for t in titik}
    ringkas = []
    for faktor, blok, na, peta, opt in sorted({k[:5] for k in terbaik}):
        asal = [faktor, 'dasar'] if faktor == 'jarak' else [faktor]
        Rs = sorted({k[5] for k in terbaik if k[0] in asal and k[1:5] == (blok, na, peta, opt)})
        Ls = sorted({k[6] for k in terbaik if k[0] == faktor and k[1:5] == (blok, na, peta, opt)})
        status = {}
        for L in Ls:
            g = [min(terbaik[(a, blok, na, peta, opt, R, L)]['terbaik_mHa'] for a in asal if (a, blok, na, peta, opt, R, L) in terbaik)
                 for R in Rs if any((a, blok, na, peta, opt, R, L) in terbaik for a in asal)]
            status[L] = 'gagal' if any(x >= 1e3 * haa.CHEM for x in g) else ('lolos' if len(g) == len(Rs) else 'kurang')
        lolos = [L for L in Ls if status[L] == 'lolos']
        L_bintang = lolos[0] if lolos else (f'>{Ls[-1]}' if set(status.values()) == {'gagal'} else 'kurang')
        ringkas.append(dict(faktor=faktor, blok=blok, na=na, peta=peta, opt=opt, R_min=Rs[0], R_maks=Rs[-1], n_R=len(Rs),
                            L_diuji=' '.join(map(str, Ls)), L_bintang=L_bintang,
                            status=' '.join(f'{L}:{status[L][0]}' for L in Ls)))
        print(*ringkas[-1].values())
    haa.tulis('tabel/beh2_lstar.csv', ringkas)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print(len(TUGAS))
    elif sys.argv[1] == 'tabel':
        tabel()
    else:
        jalan(int(sys.argv[1]))
