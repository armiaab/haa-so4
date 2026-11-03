import numpy as np
from qiskit.quantum_info import Statevector
from scipy.special import betainc

import haa

SAMPEL = 3000
DESAIN = [('real', 0), ('real', 1), ('so4', 0), ('so4', 1), ('ra', 0), ('can', 0), ('can', 1)]


def kelas(nama, ns, na, L):
    lab = haa.label(nama, ns, na, L)
    if nama == 'zeng':
        return ['u3-awal'] * (3 * ns) + lab[3 * ns:-3 * ns] + ['u3-akhir'] * (3 * ns)
    ry = ['ry-sis'] * ns + ['ry-anc'] * na
    if nama in haa.SATU:
        return ry * (L + 1)
    return (ry + [g for _ in haa.pasangan(ns, na) for g in haa.BLOK[nama]]) * L


def gradien(desain):
    rng = np.random.default_rng(7)
    baris = []
    for mol in ('H4', 'BeH2'):
        R = max(haa.JARAK[mol])
        for nama, na in desain:
            for L in (2, 4):
                m = haa.model(mol, R, nama, L, na)
                lab = np.array(kelas(nama, m['ns'], na, L))
                P = lab.size
                var = np.array([m['G'](rng.uniform(-1, 1, P)) for _ in range(200)]).var(axis=0, ddof=1)
                g0 = np.abs(m['G'](np.zeros(P)))
                for k in ['semua'] + list(dict.fromkeys(lab)):
                    s = np.ones(P, bool) if k == 'semua' else lab == k
                    baris.append(dict(molekul=mol, R=R, blok=nama, na=na, L=L, parameter=P, kelas=k, jumlah=int(s.sum()),
                                      var_rata=var[s].mean(), var_median=np.median(var[s]), g0_maks=g0[s].max(),
                                      nol_di_hf=np.mean(g0[s] < 1e-10), sampel=200, lebar=1.0))
    return baris


rng = np.random.default_rng(7)
lantai = np.mean([haa.dkl(rng.random(SAMPEL)) for _ in range(200)])
ekspresi = []
for mol in ('H4', 'BeH2'):
    H_sys, offset, hf, e_hf, e_fci = haa.hamiltonian(mol, haa.SETIMBANG[mol])
    ns = H_sys.num_qubits
    for nama, na in DESAIN:
        for L in (1, 2, 3, 4):
            qc = haa.ansatz(nama, ns, na, L, hf)
            P, N = qc.num_parameters, 2 ** (ns + na)
            s1, s2 = [np.array([Statevector(qc.assign_parameters(x)).data for x in rng.uniform(-np.pi, np.pi, (SAMPEL, P))])
                      for _ in range(2)]
            F = np.clip(np.abs(np.sum(s1.conj() * s2, axis=1)) ** 2, 0, 1)
            real = nama != 'can'
            r = haa.selang(betainc(0.5, (N - 1) / 2, F), rng) if real else (np.nan,) * 3
            c = haa.selang(1 - (1 - F) ** (N - 1), rng)
            ekspresi.append(dict(molekul=mol, blok=nama, na=na, L=L, ns=ns, qubit=ns + na, N=N, parameter=P, sampel=SAMPEL,
                                 bin=25, real=int(real), maks_imag=max(np.abs(s1.imag).max(), np.abs(s2.imag).max()),
                                 dkl_real=r[0], dkl_real_bawah=r[1], dkl_real_atas=r[2], dkl_kompleks=c[0],
                                 dkl_kompleks_bawah=c[1], dkl_kompleks_atas=c[2], lantai=lantai))
haa.tulis('ekspresi.csv', ekspresi)
haa.tulis('gradien.csv', gradien(DESAIN) + gradien([('zeng', 1), ('zeng', 0)]))
print(len(ekspresi))
