import functools
import os
import time
from multiprocessing import Pool

import numpy as np
import scipy.sparse.linalg as sla
from qiskit.quantum_info import Statevector
from qiskit_algorithms.optimizers import L_BFGS_B
from scipy.special import betainc

import haa

MOL = 'H4'
DESAIN = {'hea': ('hea', 0, 4), 'can-star': ('zeng', 1, 2), 'can-brick': ('zeng', 0, 3), 'so4-star': ('so4', 1, 2),
          'so4-brick': ('so4', 0, 2)}
LAPIS = {'hea': range(1, 7), 'can-star': range(1, 7), 'can-brick': range(1, 7), 'so4-star': range(1, 5), 'so4-brick': range(1, 5)}
R_JEJAK = (0.9, 1.5)
R_RAPAT = sorted({round(r, 2) for r in np.arange(0.5, 2.501, 0.05)} | set(haa.JARAK[MOL]))
EPS = (1e-8, 1e-4, 1e-2)
BLOKGEN = ('XX', 'YY', 'ZZ', 'XY', 'YX', 'YZ', 'ZY')
THETA = np.linspace(-np.pi, np.pi, 121)
RESTART = 12
MAXITER = 3000


@functools.lru_cache(None)
def sistem(R):
    H_sys, offset, hf, e_hf, e_fci = haa.hamiltonian(MOL, R)
    w, v = sla.eigsh(H_sys.to_matrix(sparse=True), k=4, which='SA')
    w, v = w[np.argsort(w)], v[:, np.argsort(w)]
    return dict(H_sys=H_sys, e_hf=e_hf, e_fci=e_fci, e0=w[0] + offset, celah=w[1] - w[0], dasar=v[:, w < w[0] + 1e-8],
                i_hf=int(np.argmax(np.abs(Statevector(hf).data))))


@functools.lru_cache(None)
def model(d, R, L):
    blok, na, _ = DESAIN[d]
    m = haa.model(MOL, R, blok, L, na)
    m['genap'] = np.isin(m['lab'], haa.GENAP)
    m['blokgen'] = np.isin(m['lab'], BLOKGEN)
    return m


def awal(m, jenis, k):
    x = np.zeros(len(m['lab'])) if jenis == 'nol' else np.random.default_rng(k).uniform(-1, 1, len(m['lab']))
    if jenis == 'm0':
        x[m['genap'] | m['blokgen']] = 0
    return x


def maks(v, s):
    return float(np.abs(v[s]).max()) if s.any() else 0.0


def keadaan(m, s, x):
    psi = Statevector(m['qc'].assign_parameters(x)).data
    M = psi.reshape(2 ** m['na'], 2 ** m['ns'])
    rho = M.T @ M.conj()
    fasa = psi * np.exp(-1j * np.angle(psi[np.argmax(np.abs(psi))]))
    return dict(fidelitas=np.linalg.norm(M @ s['dasar'].conj()) ** 2, kemurnian=np.trace(rho @ rho).real,
                bobot_hf=rho[s['i_hf'], s['i_hf']].real, bobot_imajiner=np.sum(fasa.imag ** 2))


def optimasi(m, x0, jejak):
    hit, memo = [0], {}

    def f(x):
        hit[0] += 1
        memo.setdefault(x.tobytes(), {})['E'] = m['E'](x)
        return memo[x.tobytes()]['E']

    def g(x):
        memo.setdefault(x.tobytes(), {})['G'] = m['G'](x)
        return memo[x.tobytes()]['G']

    def catat(x):
        d = memo.get(np.asarray(x).tobytes(), {})
        memo.clear()
        e = d['E'] if 'E' in d else m['E'](x)
        gr = d['G'] if 'G' in d else m['G'](x)
        jejak.append(dict(iterasi=len(jejak), galat=e - m['e_fci'], norma_gradien=np.linalg.norm(gr),
                          grad_genap=maks(gr, m['genap']), grad_ganjil=maks(gr, ~m['genap']),
                          sudut_genap=maks(np.asarray(x), m['genap']), nfev=hit[0]))

    catat(x0)
    t0 = time.perf_counter()
    hasil = L_BFGS_B(maxiter=MAXITER, maxfun=10 ** 8, callback=catat).minimize(f, x0, jac=g)
    gr = m['G'](hasil.x)
    return hasil.x, dict(galat=hasil.fun - m['e_fci'], galat_awal=m['E'](x0) - m['e_fci'], nit=hasil.nit, nfev=hasil.nfev,
                         njac=hasil.njev, kena_batas=int(hasil.nit >= MAXITER), detik=round(time.perf_counter() - t0, 2),
                         grad_genap=maks(gr, m['genap']), grad_ganjil=maks(gr, ~m['genap']),
                         sudut_genap=maks(hasil.x, m['genap']), sudut_blok=maks(hasil.x, m['blokgen']))


def kurva(R):
    s = sistem(R)
    baris = dict(R=R, E_hf=s['e_hf'], E_fci=s['e_fci'], E0=s['e0'], celah=s['celah'],
                 bobot_hf=np.sum(np.abs(s['dasar'][s['i_hf']]) ** 2), degenerasi=s['dasar'].shape[1])
    if R in haa.JARAK[MOL]:
        lab = s['H_sys'].paulis.to_labels()
        nY = np.array([q.count('Y') for q in lab])
        bobot = np.array([len(q) - q.count('I') for q in lab])
        baris.update(suku=len(lab), suku_Y_ganjil=int((nY % 2).sum()), koef_imag_maks=np.abs(s['H_sys'].coeffs.imag).max(),
                     suku_ber_Y=int((nY > 0).sum()), **{f'suku_bobot{b}': int((bobot == b).sum()) for b in np.unique(bobot)})
    return {'h4_kurva': [baris]}


def optim(d, R, L, jenis, k, eps):
    m = model(d, R, L)
    simpan = R in R_JEJAK and L == DESAIN[d][2]
    x0 = awal(m, jenis, k)
    if jenis == 'lolos':
        x0 = optimasi(m, awal(m, 'm0', k), [])[0] + eps * np.random.default_rng(1000 + k).uniform(-1, 1, len(x0))
    jejak = []
    x, baris = optimasi(m, x0, jejak)
    kunci = dict(desain=d, R=R, L=L, awal=jenis, benih=k, eps=eps)
    cnot, dalam = haa.biaya(DESAIN[d][0], m['ns'], m['na'], L)
    baris = dict(kunci, **baris, P=len(x), cnot=cnot, dalam=dalam, **keadaan(m, sistem(R), x))
    if not simpan:
        return {'h4_optim': [baris]}
    a = haa.arah(m['qc'], x)
    baris.update(bobot_imajiner_fasa=a.pop('bobot_imajiner'), **a)
    H = np.array([(m['G'](x + 1e-4 * e) - m['G'](x - 1e-4 * e)) / 2e-4 for e in np.eye(len(x))])
    H = (H + H.T) / 2
    g = m['genap']
    eig = {'eig': np.linalg.eigvalsh(H)}
    if g.any() and (~g).any():
        eig.update(eig_genap=np.linalg.eigvalsh(H[np.ix_(g, g)]), eig_ganjil=np.linalg.eigvalsh(H[np.ix_(~g, ~g)]))
        baris['silang_maks'] = np.abs(H[np.ix_(g, ~g)]).max()
    return {'h4_optim': [baris],
            'h4_hessian': [dict(kunci, bagian=b, indeks=i, eig=e) for b, v in eig.items() for i, e in enumerate(v)],
            'h4_sudut': [dict(kunci, k=i, kelas=m['lab'][i], sudut=t) for i, t in enumerate(x)],
            'h4_jejak': [dict(kunci, **j) for j in jejak]}


def potong(d, titik):
    m = model(d, 1.5, DESAIN[d][2])
    x = awal(m, titik, 3)
    gr = m['G'](x)
    e = m['E'](x) - m['e_fci']
    baris = []
    for k in range(len(x)):
        y = x.copy()
        for t in THETA:
            y[k] = t
            baris.append(dict(desain=d, titik=titik, galat_titik=e, k=k, kelas=m['lab'][k], sudut=x[k], gradien=gr[k], theta=t,
                              galat=m['E'](y) - m['e_fci']))
    return {'h4_potong': baris}


def gradien(d, titik):
    m = model(d, 1.5, DESAIN[d][2])
    return {'h4_gradien': [dict(desain=d, titik=titik, sampel=j, k=k, kelas=m['lab'][k], abs_gradien=abs(g))
                           for j in range(20) for k, g in enumerate(m['G'](awal(m, titik, 100 + j)))]}


def varians(d, L):
    m = model(d, 1.5, L)
    rng = np.random.default_rng(7)
    P = len(m['lab'])
    var = np.array([m['G'](rng.uniform(-1, 1, P)) for _ in range(200)]).var(axis=0, ddof=1)
    g0 = np.abs(m['G'](np.zeros(P))).max()
    baris = []
    for k in ['semua'] + list(dict.fromkeys(m['lab'])):
        s = np.ones(P, bool) if k == 'semua' else m['lab'] == k
        baris.append(dict(desain=d, L=L, P=P, g0_maks=g0, kelas=k, jumlah=int(s.sum()), var_rata=var[s].mean(),
                          var_median=np.median(var[s])))
    return {'h4_varians': baris}


def ekspresi(d, L):
    m = model(d, 0.9, L)
    rng = np.random.default_rng(11)
    P = len(m['lab'])
    s1, s2 = [np.array([Statevector(m['qc'].assign_parameters(x)).data for x in rng.uniform(-np.pi, np.pi, (3000, P))])
              for _ in range(2)]
    N = s1.shape[1]
    F = np.clip(np.abs(np.sum(s1.conj() * s2, axis=1)) ** 2, 0, 1)
    im = max(np.abs(s1.imag).max(), np.abs(s2.imag).max())
    lantai = np.mean([haa.dkl(rng.random(3000)) for _ in range(200)])
    c = haa.selang(1 - (1 - F) ** (N - 1), rng)
    r = haa.selang(betainc(0.5, (N - 1) / 2, F), rng) if im < 1e-9 else ('', '', '')
    return {'h4_ekspresi': [dict(desain=d, L=L, P=P, N=N, sampel=3000, bin=25, maks_imag=im, real=im < 1e-9, lantai=lantai,
                                 dkl_kompleks=c[0], dkl_kompleks_bawah=c[1], dkl_kompleks_atas=c[2], dkl_real=r[0],
                                 dkl_real_bawah=r[1], dkl_real_atas=r[2])]}


TUGAS = [(kurva, R) for R in R_RAPAT]
for R in haa.JARAK[MOL]:
    for d, (_, _, L) in DESAIN.items():
        TUGAS += [(optim, d, R, L, 'nol', 0, 0.0)] + [(optim, d, R, L, 'm0', k, 0.0) for k in range(RESTART)]
        if R in R_JEJAK:
            TUGAS += [(optim, d, R, L, 'acak', k, 0.0) for k in range(RESTART)]
            TUGAS += [(optim, d, R, L, 'lolos', k, e) for e in EPS for k in range(4) if d.startswith('can')]
for d, Ls in LAPIS.items():
    TUGAS += [(optim, d, 1.5, L, 'm0', k, 0.0) for L in Ls if L != DESAIN[d][2] for k in range(RESTART)]
TUGAS += [(f, d, t) for f in (potong, gradien) for d in DESAIN for t in ('m0', 'acak')]
TUGAS += [(varians, d, L) for d in DESAIN for L in (1, 2, 3, 4, 6)]
TUGAS += [(ekspresi, d, L) for d in DESAIN for L in (1, 2, 3, 4)]


def jalan(t):
    return t[0](*t[1:])


if __name__ == '__main__':
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))) as pool:
        hasil = pool.map(jalan, TUGAS, chunksize=1)
    tabel = {}
    for h in hasil:
        for nama, baris in h.items():
            tabel.setdefault(nama, []).extend(baris)
    for nama, baris in tabel.items():
        haa.tulis(f'{nama}.csv', baris)
    print(len(TUGAS))
