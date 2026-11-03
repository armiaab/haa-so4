import glob
import os

import numpy as np

import haa

MOL = ['H4', 'LiH', 'HF', 'BeH2']
BLOK = ['can', 'zeng', 'real', 'so4']
TOPO = {1: 'star', 0: 'brick'}
CHEM = haa.CHEM * 1000
AMBANG = [CHEM, 1.0, 0.1]
KS = [1, 2, 3, 4, 6, 8, 10, 12]

data = {}
for p in glob.glob(f'{haa.HASIL}/pec/*_taper_*.csv'):
    b = haa.baca(os.path.relpath(p, haa.HASIL))
    data[(b[0]['molekul'], int(b[0]['na']), b[0]['blok'], int(b[0]['L']), round(float(b[0]['R']), 3))] = b


def seri(mol, na, blok):
    t = {(k[3], k[4]): np.array([float(x['galat_mHa']) for x in v]) for k, v in data.items() if k[:3] == (mol, na, blok)}
    return t, sorted({L for L, _ in t}), sorted({R for _, R in t})


def lengkap(t, Rs, L):
    return all((L, R) in t for R in Rs)


def L_bintang(t, Ls, Rs, ambang=CHEM):
    for L in Ls:
        if lengkap(t, Rs, L) and max(t[(L, R)].min() for R in Rs) < ambang:
            return L, str(L) if L == 1 or L - 1 in Ls else f'≤{L}'
    return None, f'>{max(Ls)}'


biaya = []
for mol in MOL:
    for na in (1, 0):
        for blok in BLOK + ['hea']:
            t, Ls, Rs = seri(mol, na, blok)
            if not t:
                continue
            L, teks = L_bintang(t, Ls, Rs)
            b0 = data[(mol, na, blok, L, Rs[0])][0]
            N2q, D2q = haa.biaya(blok, int(b0['ns']), na, L)
            batas = [np.mean([int(x.get('kena_batas') or 0) for x in data[(mol, na, blok, L, R)]]) for R in Rs]
            biaya.append(dict(molekul=mol, topologi=TOPO[na], na=na, blok=blok, L=L, L_teks=teks, P=int(b0['parameter']),
                              N2q=N2q, D2q=D2q, qubit=int(b0['ns']) + na, terburuk_mHa=max(t[(L, R)].min() for R in Rs),
                              sukses=np.mean(np.concatenate([t[(L, R)] for R in Rs]) < CHEM), capped=np.mean(batas)))
haa.tulis('tabel/tab2_biaya.csv', biaya)
d = {(b['molekul'], b['na'], b['blok']): b for b in biaya}

impas = []
for mol in MOL:
    for blok in BLOK:
        b0, b1 = d[(mol, 0, blok)], d[(mol, 1, blok)]
        ns = b1['qubit'] - 1
        rasio = b0['L'] / b1['L']
        impas.append(dict(molekul=mol, ns=ns, blok=blok, L0=b0['L'], L1=b1['L'], rasio=rasio, impas_N2q=ns / (ns - 1),
                          impas_D2q=ns / 2, N2q_0=b0['N2q'], N2q_1=b1['N2q'], D2q_0=b0['D2q'], D2q_1=b1['D2q'],
                          ancilla_menang_N2q=b1['N2q'] < b0['N2q'], ancilla_menang_D2q=b1['D2q'] < b0['D2q'],
                          ramalan_N2q=rasio > ns / (ns - 1), ramalan_D2q=rasio > ns / 2))
haa.tulis('tabel/tab3_impas.csv', impas)

pembanding = []
for mol in MOL:
    can = min((d[(mol, na, b)] for na in (1, 0) for b in ('can', 'zeng')), key=lambda b: (b['N2q'], b['D2q'], -b['sukses']))
    for peran, b in (('hea', d[(mol, 0, 'hea')]), ('can_terbaik', can), ('so4_star', d[(mol, 1, 'so4')]),
                     ('so4_brick', d[(mol, 0, 'so4')])):
        pembanding.append(dict(molekul=mol, peran=peran, blok=b['blok'], topologi=b['topologi'], L=b['L'], P=b['P'],
                               N2q=b['N2q'], D2q=b['D2q'], qubit=b['qubit'], terburuk_mHa=b['terburuk_mHa'], sukses=b['sukses']))
haa.tulis('tabel/tab4_pembanding.csv', pembanding)

ketahanan, pergeometri, start, hampir = [], [], [], []
for mol in MOL:
    for na in (1, 0):
        for blok in BLOK:
            t, Ls, Rs = seri(mol, na, blok)
            baris = dict(molekul=mol, topologi=TOPO[na], blok=blok)
            for a in AMBANG:
                L, teks = L_bintang(t, Ls, Rs, a)
                baris[f'L_{a:g}'] = teks
                baris[f'N2q_{a:g}'] = haa.biaya(blok, d[(mol, na, blok)]['qubit'] - na, na, L)[0] if L else ''
            ketahanan.append(baris)
rng = np.random.default_rng(0)
for b in biaya:
    mol, na, blok, L = b['molekul'], b['na'], b['blok'], b['L']
    t, Ls, Rs = seri(mol, na, blok)
    for R in Rs:
        f = data[(mol, na, blok, L, R)]
        pergeometri.append(dict(molekul=mol, topologi=TOPO[na], blok=blok, L=L, R=R, terbaik_mHa=t[(L, R)].min(),
                                median_mHa=np.median(t[(L, R)]), sukses=int((t[(L, R)] < CHEM).sum()), n=len(f),
                                capped=sum(int(x.get('kena_batas') or 0) for x in f),
                                detik_median=np.median([float(x['detik']) for x in f])))
    Lp = [x for x in Ls if lengkap(t, Rs, x)]
    n = min(t[(x, R)].size for x in Lp for R in Rs)
    baris = dict(molekul=mol, topologi=TOPO[na], blok=blok, L12=L)
    for k in KS:
        Lk = np.full(4000, np.inf)
        for x in sorted(Lp, reverse=True):
            e = np.array([t[(x, R)][:n] for R in Rs])
            pilih = rng.random((4000, e.shape[0], n)).argsort(axis=2)[:, :, :k]
            terbaik = np.take_along_axis(np.broadcast_to(e, (4000,) + e.shape), pilih, 2).min(2)
            Lk[(terbaik < CHEM).all(axis=1)] = x
        baris.update({f'sama_k{k}': np.mean(Lk == L), f'lebih_k{k}': np.mean(Lk > L), f'kurang_k{k}': np.mean(Lk < L)})
    start.append(baris)
    if L - 1 in Ls and lengkap(t, Rs, L - 1):
        gagal = [t[(L - 1, R)].min() for R in Rs if t[(L - 1, R)].min() >= CHEM]
        hampir.append(dict(molekul=mol, topologi=TOPO[na], blok=blok, L=L, n_gagal=len(gagal), n_R=len(Rs), terburuk=max(gagal),
                           terdekat=min(gagal)))
haa.tulis('tabel/tab5_ketahanan.csv', ketahanan)
haa.tulis('tabel/tabS1_pergeometri.csv', pergeometri)
haa.tulis('tabel/start.csv', start)
haa.tulis('tabel/hampir.csv', hampir)

konfigurasi = {(k['molekul'], round(float(k['R']), 3)): k for k in haa.baca('sistem_konfigurasi.csv') if k['urutan'] == '1'}
eksitasi = []
for g in haa.baca('sistem.csv'):
    mol, R = g['molekul'], round(float(g['R']), 3)
    if mol not in MOL:
        continue
    baris = dict(molekul=mol, R=R, korelasi_mHa=g['korelasi_mHa'], w_hf=g['w_hf'], n99=g['n99'], k_ca=g['k_ca'],
                 dominan=konfigurasi[(mol, R)]['balik'], c_dominan=konfigurasi[(mol, R)]['c'])
    for na in (1, 0):
        for blok in BLOK + ['hea']:
            t, Ls, Rs = seri(mol, na, blok)
            if t:
                baris[f'L_{blok}_{TOPO[na]}'] = next((L for L in Ls if (L, R) in t and t[(L, R)].min() < CHEM), '')
    eksitasi.append(baris)
haa.tulis('tabel/tabS3_eksitasi.csv', eksitasi)
print(len(biaya))
