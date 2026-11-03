import functools
import itertools
import os
import time
from multiprocessing import Pool

import numpy as np
from qiskit import transpile
from qiskit.converters import circuit_to_dag
from qiskit.quantum_info import DensityMatrix, Kraus, Operator, Pauli, Statevector, SuperOp
from qiskit_algorithms.optimizers import L_BFGS_B

import haa

P_GRID = [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
CHEM = haa.CHEM * 1000
KONFIGURASI = [(b['molekul'], R, geometri, b['peran'], b['blok'], int(b['L']), int(b['topologi'] == 'star'))
               for b in haa.baca('tabel/tab4_pembanding.csv')
               for geometri, R in (('setimbang', haa.SETIMBANG[b['molekul']]), ('regang', max(haa.JARAK[b['molekul']])))]


@functools.lru_cache(None)
def depol(p, k):
    pauli = [Pauli(''.join(s)).to_matrix() for s in itertools.product('IXYZ', repeat=k)]
    return SuperOp(Kraus([np.sqrt(1 - p + p / 4 ** k) * pauli[0]] + [np.sqrt(p / 4 ** k) * a for a in pauli[1:]]))


def jadwal(qc):
    lapis = []
    for lap in circuit_to_dag(qc).layers():
        gerbang = [(SuperOp(Operator(node.op)), [qc.find_bit(q).index for q in node.qargs]) for node in lap['graph'].op_nodes()]
        sibuk = {q for _, qs in gerbang for q in qs}
        dua = any(len(qs) == 2 for _, qs in gerbang)
        lapis.append((gerbang, [q for q in range(qc.num_qubits) if q not in sibuk] if dua else []))
    return lapis


def energi(lapis, H, p, diam):
    rho = DensityMatrix.from_int(0, 2 ** H.num_qubits)
    for gerbang, bebas in lapis:
        for op, qs in gerbang:
            rho = rho.evolve(op.compose(depol(p if len(qs) == 2 else p / 10, len(qs))), qs)
        for q in bebas if diam else []:
            rho = rho.evolve(depol(p / 10, 1), [q])
    return float(np.real(rho.expectation_value(H)))


def p_bintang(f):
    a, b = np.log(1e-8), np.log(1e-1)
    for _ in range(30):
        c = (a + b) / 2
        a, b = (c, b) if f(np.exp(c)) < CHEM else (a, c)
    return float(np.exp((a + b) / 2))


def jalan(k):
    mol, R, geometri, peran, ansatz, L, na = k
    t0 = time.perf_counter()
    terbaik = min(haa.baca(f'pec/{mol}_taper_na{na}_R{R:.3f}_{ansatz}_L{L}.csv'), key=lambda b: float(b['galat_mHa']))
    m = haa.model(mol, R, ansatz, L, na)
    P = len(m['lab'])
    x0 = np.random.default_rng(int(terbaik['restart'])).uniform(-1, 1, P)
    hasil = L_BFGS_B(maxiter=3000, maxfun=10 ** 8).minimize(m['E'], x0, jac=m['G'])
    qc = transpile(haa.ansatz(ansatz, m['ns'], na, L, m['hf'], hasil.x, kompilasi=True), basis_gates=['cx', 'u'],
                   optimization_level=1, seed_transpiler=7)
    lapis = jadwal(qc)

    def galat(p, diam=False):
        return (energi(lapis, m['H'], p, diam) + m['offset'] - m['e_fci']) * 1000

    baris = dict(molekul=mol, R=R, geometri=geometri, peran=peran, blok=ansatz, topologi='star' if na else 'brick', L=L,
                 parameter=P, cnot=qc.count_ops().get('cx', 0), cnot_rumus=haa.biaya(ansatz, m['ns'], na, L)[0],
                 satu_qubit=qc.count_ops().get('u', 0), kedalaman2q=qc.depth(lambda i: i.operation.num_qubits == 2),
                 galat_produksi=float(terbaik['galat_mHa']), galat_optimasi=(hasil.fun - m['e_fci']) * 1000,
                 galat_p0=galat(0.0),
                 galat_sv=(float(np.real(Statevector(qc).expectation_value(m['H']))) + m['offset'] - m['e_fci']) * 1000)
    for p in P_GRID:
        baris[f'galat_p{p:g}'] = galat(p)
        baris[f'galat_diam_p{p:g}'] = galat(p, True)
    for model in ('', '_diam'):
        baris[f'kepekaan{model}'] = (baris[f'galat{model}_p1e-05'] - baris['galat_p0']) / 1e-5
        baris[f'p_derau{model}'] = CHEM / baris[f'kepekaan{model}']
    baris['p_bintang'] = p_bintang(galat)
    baris['p_bintang_diam'] = p_bintang(lambda p: galat(p, True))
    baris['detik'] = round(time.perf_counter() - t0, 1)
    sudut = [dict(molekul=mol, R=R, geometri=geometri, peran=peran, blok=ansatz, na=na, L=L, restart=terbaik['restart'],
                  k=i, theta=t) for i, t in enumerate(hasil.x)]
    return baris, sudut


if __name__ == '__main__':
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))) as pool:
        hasil = pool.map(jalan, KONFIGURASI)
    haa.tulis('derau.csv', [h[0] for h in hasil])
    haa.tulis('derau_sudut.csv', sum((h[1] for h in hasil), []))
    print(len(KONFIGURASI))
