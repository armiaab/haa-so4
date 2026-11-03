import functools
import glob
import os
import sys
import time
from collections import defaultdict

import numpy as np
from qiskit.primitives import StatevectorEstimator
from qiskit_algorithms.optimizers import COBYLA, L_BFGS_B, SPSA
from qiskit_algorithms.utils import algorithm_globals

import haa

MOL = ('H4', 'LiH', 'HF', 'BeH2')
RESTART = 12
TAMBAH = {'H4': [1.8, 2.0, 2.25], 'LiH': [3.6, 4.0], 'HF': [1.9, 2.2, 2.5], 'BeH2': [0.7, 2.5, 3.0, 3.3]}
# Exact rule-2 endpoints (2.5 R0) that TAMBAH only approximated; computed by ujung.py.
UJUNG = {'LiH': [3.9875], 'BeH2': [3.315]}
# Layer caps of rule 5 (aturan_keputusan.md); the R_Y R_Z HEA on BeH2 has 12.
BATAS = {'so4': 6, 'zeng': 6, 'hea': 10, 'heau3': 10, 'ra': 10}


def batas(blok, mol):
    return 12 if (blok, mol) == ('hea', 'BeH2') else BATAS.get(blok)


def daftar(mol):
    """Registered range of rule 2: original grid plus added points up to exactly 2.5 R0."""
    ujung = round(2.5 * haa.SETIMBANG[mol], 4)
    tambah = TAMBAH[mol] + UJUNG.get(mol, [])
    return sorted(set(haa.JARAK[mol]) | {R for R in tambah if min(haa.JARAK[mol]) <= R <= ujung + 1e-9})
HAA = [('so4', 1), ('so4', 0), ('zeng', 1), ('zeng', 0)]
HEA = ['hea', 'heau3', 'ra']
KONTROL = [('canu3', 1), ('canu3', 0), ('so4u3', 1), ('so4u3', 0)]
H4 = {'can-star': ('zeng', 1, 2), 'can-brick': ('zeng', 0, 3), 'so4-star': ('so4', 1, 2), 'so4-brick': ('so4', 0, 2),
      'hea': ('hea', 0, 4)}
BLOKGEN = ('XX', 'YY', 'ZZ', 'XY', 'YX', 'YZ', 'ZY')
CHEM = 1e3 * haa.CHEM

TUGAS = []
for baris in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tugas_pelengkap.txt')).read().split('\n'):
    if baris:
        mol, R, b, na, L = baris.split()
        TUGAS.append(('pelengkap', mol, float(R), b, int(na), int(L)))
for mol in ('H4', 'LiH', 'HF'):
    TUGAS += [('heau3', mol, R, 'heau3', 0, L) for R in haa.JARAK[mol] for L in range(1, 11)]
for mol in MOL:
    TUGAS += [('kontrol', mol, R, b, na, L) for R in haa.JARAK[mol] for b, na in KONTROL for L in range(1, 5)]
for mol in MOL:
    TUGAS += [('jarak', mol, R, b, na, L) for R in TAMBAH[mol] for b, na in ([('zeng', 0)] if mol == 'BeH2' else HAA)
              for L in range(1, 7)]
    TUGAS += [('jarak', mol, R, b, 0, L) for R in TAMBAH[mol] for b in (HEA[1:] if mol == 'BeH2' else HEA)
              for L in range(1, 11)]
for mol in MOL:
    TUGAS += [('ra', mol, R, 'ra', 0, L) for R in haa.JARAK[mol] for L in range(1, 11)]
TUGAS = [t + (k,) for t in TUGAS for k in range(RESTART)]
for d in H4:
    for R in (0.9, 1.5):
        TUGAS += [('awal', d, R, 'kecil', e, 'lbfgsb', k) for e in (0.01, 0.1) for k in range(RESTART)]
        TUGAS += [('awal', d, R, 'nol', 0, o, 0) for o in ('lbfgsb', 'cobyla')]
        TUGAS += [('awal', d, R, 'm0', 0, o, k) for o in ('cobyla', 'spsa', 'spsa-shot') for k in range(RESTART)]
        TUGAS += [('awal', d, R, 'nol', 0, 'spsa', k) for k in range(RESTART)]


def pec(bagian, mol, R, blok, na, L, k):
    nama = f'uji/{bagian}_{mol}_{blok}_na{na}_R{R:.3f}_L{L}_r{k}.csv'
    if os.path.exists(os.path.join(haa.HASIL, nama)):
        return
    m = haa.model(mol, R, blok, L, na)
    P = len(m['lab'])
    cnot, dalam = haa.biaya(blok, m['ns'], na, L)
    x0 = np.random.default_rng(k).uniform(-1, 1, P)
    t0 = time.perf_counter()
    hasil = L_BFGS_B(maxiter=3000, maxfun=10 ** 8).minimize(m['E'], x0, jac=m['G'])
    e = m['E'](hasil.x)
    haa.tulis(nama, [dict(bagian=bagian, molekul=mol, R=R, blok=blok, L=L, na=na, ns=m['ns'], parameter=P, cnot=cnot,
                          kedalaman=dalam, restart=k, nit=hasil.nit, nfev=hasil.nfev, kena_batas=int(hasil.nit >= 3000),
                          detik=round(time.perf_counter() - t0, 1), E_vqe=e, E_hf=m['e_hf'], E_fci=m['e_fci'],
                          E_dasar=haa.dasar(m['H_sys']) + m['offset'], galat_mHa=(e - m['e_fci']) * 1000)])
    print(nama, (e - m['e_fci']) * 1000)


def awal(d, R, jenis, eps, opt, k):
    nama = f'uji/awal_{d}_R{R:.3f}_{jenis}{eps:g}_{opt}_r{k}.csv'
    if os.path.exists(os.path.join(haa.HASIL, nama)):
        return
    blok, na, L = H4[d]
    m = haa.model('H4', R, blok, L, na)
    P = len(m['lab'])
    if jenis == 'nol':
        x0 = np.zeros(P)
    elif jenis == 'm0':
        x0 = np.random.default_rng(k).uniform(-1, 1, P)
        x0[np.isin(m['lab'], haa.GENAP + BLOKGEN)] = 0
    else:
        x0 = np.random.default_rng(k).uniform(-eps, eps, P)
    E = m['E']
    if opt == 'spsa-shot':
        bising = StatevectorEstimator(default_precision=1e-3, seed=np.random.default_rng(1000 + k))
        E = lambda x: float(bising.run([(m['qc'], m['H'], x)]).result()[0].data.evs) + m['offset']
    algorithm_globals.random_seed = k
    t0 = time.perf_counter()
    if opt == 'lbfgsb':
        hasil = L_BFGS_B(maxiter=3000, maxfun=10 ** 8).minimize(E, x0, jac=m['G'])
    elif opt == 'cobyla':
        hasil = COBYLA(maxiter=10000).minimize(E, x0)
    else:
        hasil = SPSA(maxiter=5000).minimize(E, x0)
    genap = np.isin(m['lab'], haa.GENAP)
    e = m['E'](hasil.x)
    haa.tulis(nama, [dict(desain=d, R=R, awal=jenis, eps=eps, opt=opt, restart=k, galat_awal_mHa=(m['E'](x0) - m['e_fci']) * 1000,
                          galat_mHa=(e - m['e_fci']) * 1000, nit=hasil.nit or '', nfev=hasil.nfev,
                          sudut_genap=float(np.abs(hasil.x[genap]).max()) if genap.any() else 0.0,
                          detik=round(time.perf_counter() - t0, 1))])
    print(nama, (e - m['e_fci']) * 1000)


def jalan(i):
    t = TUGAS[i - 1]
    awal(*t[1:]) if t[0] == 'awal' else pec(*t)


@functools.lru_cache(None)
def biaya(blok, ns, na, L):
    return haa.biaya(blok, ns, na, L)


def tabel():
    data = defaultdict(list)
    sumber = glob.glob(os.path.join(haa.HASIL, 'pec', '*_taper_*.csv')) + glob.glob(os.path.join(haa.HASIL, 'beh2', '*.csv'))
    sumber += [f for f in glob.glob(os.path.join(haa.HASIL, 'uji', '*.csv')) if not os.path.basename(f).startswith(('awal_', 'cek_'))]
    ns = {}
    for f in sumber:
        for b in haa.baca(f):
            if b.get('peta', 'taper') != 'taper' or b.get('opt', 'lbfgsb') != 'lbfgsb' or b['molekul'] not in MOL:
                continue
            k = (b['molekul'], b['blok'], int(b['na']))
            ns[k] = int(b['ns'])
            data[k + (float(b['R']), int(b['L']))].append(float(b['galat_mHa']))
    # Rule 7 (aturan_keputusan.md): a geometry whose RHF did not converge, or whose tapered ground state lies more
    # than 1e-8 Ha from E_FCI, is removed for all designs. The checks come from ujung.py / ujung_cluster.py.
    buang = {(b['molekul'], float(b['R'])) for f in glob.glob(os.path.join(haa.HASIL, 'uji', 'cek_*.csv'))
             for b in haa.baca(f) if b['buang'] == '1'}
    for mol, R in sorted(buang):
        print('rule 7: removed for all designs:', mol, R)
    ringkas = []
    for k in sorted({k[:3] for k in data}):
        mol = k[0]
        luas = sorted(set(haa.JARAK[mol]) | {R for R in TAMBAH[mol] if R >= min(haa.JARAK[mol])})
        cap = batas(k[1], mol)
        for rentang, Rs in (('asli', haa.JARAK[mol]), ('luas', luas), ('daftar', daftar(mol))):
            Rs = [R for R in Rs if (mol, R) not in buang]
            titik = [x for x in data if x[:3] == k and x[3] in Rs]
            tambahan = set(Rs) - set(haa.JARAK[mol])
            if rentang != 'asli' and not any(x[3] in tambahan for x in titik):
                continue
            Ls = sorted({x[4] for x in titik if cap is None or x[4] <= cap})
            baris = dict(molekul=mol, blok=k[1], na=k[2], rentang=rentang, R_maks=max(Rs), n_R=len(Rs), L_diuji=' '.join(map(str, Ls)),
                         n_min=min(len(data[x]) for x in titik))
            for nama, syarat in (('terbaik', 1), ('tiga', 3)):
                status = {}
                for L in Ls:
                    ada = [R for R in Rs if k + (R, L) in data]
                    gagal = any(sum(g < CHEM for g in data[k + (R, L)]) < syarat for R in ada)
                    status[L] = 'gagal' if gagal else ('lolos' if len(ada) == len(Rs) else 'kurang')
                lolos = [L for L in Ls if status[L] == 'lolos']
                baris[f'L_{nama}'] = lolos[0] if lolos else (f'>{Ls[-1]}' if set(status.values()) == {'gagal'} else 'kurang')
                baris[f'cnot_{nama}'], baris[f'kedalaman_{nama}'] = biaya(k[1], ns[k], k[2], lolos[0]) if lolos else ('', '')
                # Worst-geometry error (mHa) at the last tested layer count for a design that fails: the best start
                # for 'terbaik', the third-best start for 'tiga'.
                Lc = Ls[-1] if Ls else None
                lengkap = Lc is not None and all(k + (R, Lc) in data and len(data[k + (R, Lc)]) >= syarat for R in Rs)
                baris[f'galat_batas_{nama}'] = (max(sorted(data[k + (R, Lc)])[syarat - 1] for R in Rs)
                                                if lengkap and not lolos else '')
            ringkas.append(baris)
            print(*baris.values())
    haa.tulis('tabel/uji_lstar.csv', ringkas)
    awal_ = defaultdict(list)
    for f in glob.glob(os.path.join(haa.HASIL, 'uji', 'awal_*.csv')):
        for b in haa.baca(f):
            awal_[(b['desain'], float(b['R']), b['awal'], float(b['eps']), b['opt'])].append(float(b['galat_mHa']))
    for b in haa.baca('h4_optim.csv'):
        if b['awal'] in ('m0', 'acak') and float(b['R']) in (0.9, 1.5) and int(b['L']) == H4[b['desain']][2]:
            awal_[(b['desain'], float(b['R']), b['awal'], 0.0, 'lbfgsb')].append(1e3 * float(b['galat']))
    hasil = [dict(desain=d, R=R, awal=a, eps=e, opt=o, n=len(g), sukses=sum(x < CHEM for x in g), terbaik_mHa=min(g),
                  median_mHa=float(np.median(g))) for (d, R, a, e, o), g in sorted(awal_.items())]
    haa.tulis('tabel/uji_awal.csv', hasil)
    for h in hasil:
        print(*h.values())


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print(len(TUGAS))
    elif sys.argv[1] == 'tabel':
        tabel()
    else:
        jalan(int(sys.argv[1]))
