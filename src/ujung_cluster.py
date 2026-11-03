#!/usr/bin/env python
"""HAA-SO(4): exact rule-2 endpoints and rule-7 checks. Self-contained (no project imports) for a SLURM array job.

The extended grids end at rounded endpoints for two molecules (LiH 4.0 A instead of 2.5 R0 = 3.9875 A, BeH2 3.3 A
instead of 3.315 A). This script computes the exact endpoints for every design of tab:lstar, with the layer caps of
rule 5 and 12 starts, and checks every geometry of the registered range for RHF convergence and for the distance
between the lowest eigenvalue of the tapered Hamiltonian and E_FCI (rule 7).

Usage
    python ujung_cluster.py          number of tasks (1358)
    python ujung_cluster.py list     all tasks, one per line
    python ujung_cluster.py <i>      run task i (1-based); does nothing if its output file already exists

Output: $OUT/uji/ujung_*.csv (one per start) and $OUT/uji/cek_*.csv (one per geometry), with the same columns as
results/uji/ from uji.py. $OUT defaults to ./hasil_ujung next to this file. Copy both kinds of file into
results/uji/ of the project afterwards.

The Hamiltonian, ansatz, cost and optimizer code is copied verbatim from src/haa.py and src/uji.py, so the results
match the earlier runs bit for bit with the same package versions: numpy 2.5.3, scipy 1.18.1, qiskit 2.5.2,
qiskit-algorithms 0.4.0, qiskit-nature 0.8.0, pyscf 2.14.0.
"""
import csv
import os
import sys
import time

import numpy as np
import scipy.sparse.linalg as sla
from pyscf import gto, mcscf, scf
from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector
from qiskit.circuit.library import CXGate
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import Operator, SparsePauliOp
from qiskit.synthesis import TwoQubitBasisDecomposer
from qiskit_algorithms.gradients import ReverseEstimatorGradient
from qiskit_algorithms.optimizers import L_BFGS_B
from qiskit_nature.second_q.circuit.library import HartreeFock
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.transformers import FreezeCoreTransformer

OUT = os.environ.get('OUT', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hasil_ujung'))

# ---- from src/haa.py ------------------------------------------------------------------------------------------
GEOMETRI = {
    'H4': lambda R: '; '.join(f'H 0 0 {i * R:.6f}' for i in range(4)),
    'LiH': lambda R: f'Li 0 0 0; H 0 0 {R}',
    'HF': lambda R: f'H 0 0 0; F 0 0 {R}',
    'BeH2': lambda R: f'Be 0 0 0; H 0 0 {R}; H 0 0 {-R}',
}
JARAK = {
    'H4': [0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5],
    'LiH': [1.0, 1.3, 1.595, 2.0, 2.6, 3.2],
    'HF': [0.7, 0.8, 0.9, 1.0, 1.1, 1.3, 1.6],
    'BeH2': [1.0, 1.2, 1.326, 1.5, 1.7, 2.0],
}
SETIMBANG = {'H4': 0.9, 'LiH': 1.595, 'HF': 1.0, 'BeH2': 1.326}
BLOK = {'zeng': ('ZZ', 'XX', 'YY'), 'so4': ('XY', 'YX', 'YZ', 'ZY')}
SATU = {'hea': ('ry', 'rz'), 'ra': ('ry',), 'heau3': ('rz', 'ry', 'rz')}
decom = TwoQubitBasisDecomposer(CXGate())


def tulis(nama, baris):
    p = os.path.join(OUT, nama)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for b in baris for k in b)), restval='')
        w.writeheader()
        w.writerows(baris)


def hamiltonian(mol, R, peta='taper'):
    atom = GEOMETRI[mol](R)
    problem = FreezeCoreTransformer().transform(PySCFDriver(atom=atom, basis='sto3g').run())
    mapper = problem.get_tapered_mapper(JordanWignerMapper()) if peta == 'taper' else JordanWignerMapper()
    H_sys = mapper.map(problem.hamiltonian.second_q_op())
    offset = sum(problem.hamiltonian.constants.values())
    hf = HartreeFock(problem.num_spatial_orbitals, problem.num_particles, mapper)
    mf = scf.RHF(gto.M(atom=atom, basis='sto3g', verbose=0)).run()
    cas = mcscf.CASCI(mf, problem.num_spatial_orbitals, problem.num_particles)
    cas.verbose = 0
    return H_sys, offset, hf, float(mf.e_tot), float(cas.kernel()[0])


def dasar(H_sys):
    m = H_sys.to_matrix(sparse=True)
    if m.shape[0] <= 4:
        return float(np.linalg.eigvalsh(m.toarray())[0])
    return float(sla.eigsh(m, k=1, which='SA', maxiter=100000)[0][0])


def pasangan(ns, na):
    if na:
        return [(i, ns + i % na) for i in range(ns)]
    return [(i, i + 1) for i in range(0, ns - 1, 2)] + [(i, i + 1) for i in range(1, ns - 1, 2)]


def rot(qc, basis, sudut, a, b):
    for q, s in zip((a, b), basis):
        if s == 'X':
            qc.h(q)
        if s == 'Y':
            qc.rx(np.pi / 2, q)
    qc.cx(a, b)
    qc.rz(sudut, b)
    qc.cx(a, b)
    for q, s in zip((a, b), basis):
        if s == 'X':
            qc.h(q)
        if s == 'Y':
            qc.rx(-np.pi / 2, q)


def label(nama, ns, na, L):
    n = ns + na
    if nama in SATU:
        return [g for _ in range(L + 1) for g in SATU[nama] for _ in range(n)]
    blok = [g for _ in pasangan(ns, na) for g in BLOK[nama]]
    if nama == 'zeng':
        return ['rz', 'ry', 'rz'] * ns + blok * L + ['rz', 'ry', 'rz'] * ns
    if nama.endswith('u3'):
        return (['rz', 'ry', 'rz'] * n + blok) * L
    return (['ry'] * n + blok) * L


def u3(qc, t, k, ns):
    for q in range(ns):
        qc.rz(t[k], q)
        qc.ry(t[k + 1], q)
        qc.rz(t[k + 2], q)
        k += 3
    return k


def ansatz(nama, ns, na, L, hf=None, sudut=None, kompilasi=False):
    t = ParameterVector('t', len(label(nama, ns, na, L))) if sudut is None else sudut
    n = ns + na
    qc = QuantumCircuit(n)
    if hf is not None:
        qc.compose(hf, range(ns), inplace=True)
    k = 0
    if nama in SATU:
        for lapis in range(L + 1):
            for g in SATU[nama]:
                for q in range(n):
                    getattr(qc, g)(t[k], q)
                    k += 1
            if lapis < L:
                for q in range(n - 1):
                    qc.cx(q, q + 1)
        return qc
    if nama == 'zeng':
        k = u3(qc, t, k, ns)
    for _ in range(L):
        if nama.endswith('u3'):
            k = u3(qc, t, k, n)
        elif nama != 'zeng':
            for q in range(n):
                qc.ry(t[k], q)
                k += 1
        for a, b in pasangan(ns, na):
            sub = QuantumCircuit(2)
            for basis in BLOK[nama]:
                rot(sub, basis, t[k], 0, 1)
                k += 1
            qc.compose(decom(Operator(sub)) if kompilasi else sub, [a, b], inplace=True)
    if nama == 'zeng':
        u3(qc, t, k, ns)
    return qc


def biaya(nama, ns, na, L):
    qc = ansatz(nama, ns, na, L, sudut=np.full(len(label(nama, ns, na, L)), 0.7), kompilasi=True)
    return qc.count_ops().get('cx', 0), qc.depth(lambda i: i.operation.num_qubits == 2)


def model(mol, R, nama, L, na, peta='taper', konstanta=True):
    H_sys, offset, hf, e_hf, e_fci = hamiltonian(mol, R, peta)
    ns = H_sys.num_qubits
    qc = ansatz(nama, ns, na, L, hf)
    H = SparsePauliOp('I' * na).tensor(H_sys)
    est, rev = StatevectorEstimator(), ReverseEstimatorGradient()

    def E(x):
        return float(est.run([(qc, H, x)]).result()[0].data.evs) + (offset if konstanta else 0)

    def G(x):
        return np.asarray(rev.run([qc], [H], [x]).result().gradients[0], dtype=float)

    return dict(H_sys=H_sys, H=H, offset=offset, hf=hf, e_hf=e_hf, e_fci=e_fci, ns=ns, na=na, qc=qc, E=E, G=G,
                lab=np.array(label(nama, ns, na, L)))


# ---- from src/uji.py ------------------------------------------------------------------------------------------
MOL = ('H4', 'LiH', 'HF', 'BeH2')
RESTART = 12
TAMBAH = {'H4': [1.8, 2.0, 2.25], 'LiH': [3.6, 4.0], 'HF': [1.9, 2.2, 2.5], 'BeH2': [0.7, 2.5, 3.0, 3.3]}
UJUNG = {'LiH': [3.9875], 'BeH2': [3.315]}
BATAS = {'so4': 6, 'zeng': 6, 'hea': 10, 'heau3': 10, 'ra': 10}
HAA = [('so4', 1), ('so4', 0), ('zeng', 1), ('zeng', 0)]
HEA = ['hea', 'heau3', 'ra']


def batas(blok, mol):
    return 12 if (blok, mol) == ('hea', 'BeH2') else BATAS.get(blok)


def daftar(mol):
    """Registered range of rule 2: original grid plus added points up to exactly 2.5 R0."""
    ujung = round(2.5 * SETIMBANG[mol], 4)
    tambah = TAMBAH[mol] + UJUNG.get(mol, [])
    return sorted(set(JARAK[mol]) | {R for R in tambah if min(JARAK[mol]) <= R <= ujung + 1e-9})


def pec(bagian, mol, R, blok, na, L, k):
    nama = f'uji/{bagian}_{mol}_{blok}_na{na}_R{R:.3f}_L{L}_r{k}.csv'
    if os.path.exists(os.path.join(OUT, nama)):
        return
    m = model(mol, R, blok, L, na)
    P = len(m['lab'])
    cnot, dalam = biaya(blok, m['ns'], na, L)
    x0 = np.random.default_rng(k).uniform(-1, 1, P)
    t0 = time.perf_counter()
    hasil = L_BFGS_B(maxiter=3000, maxfun=10 ** 8).minimize(m['E'], x0, jac=m['G'])
    e = m['E'](hasil.x)
    tulis(nama, [dict(bagian=bagian, molekul=mol, R=R, blok=blok, L=L, na=na, ns=m['ns'], parameter=P, cnot=cnot,
                      kedalaman=dalam, restart=k, nit=hasil.nit, nfev=hasil.nfev, kena_batas=int(hasil.nit >= 3000),
                      detik=round(time.perf_counter() - t0, 1), E_vqe=e, E_hf=m['e_hf'], E_fci=m['e_fci'],
                      E_dasar=dasar(m['H_sys']) + m['offset'], galat_mHa=(e - m['e_fci']) * 1000)])
    print(nama, (e - m['e_fci']) * 1000)


# ---- rule-7 check -----------------------------------------------------------------------------------------------
def cek(mol, R):
    nama = f'uji/cek_{mol}_R{R:.4f}.csv'
    if os.path.exists(os.path.join(OUT, nama)):
        return
    atom = GEOMETRI[mol](R)
    mf = scf.RHF(gto.M(atom=atom, basis='sto3g', verbose=0))
    mf.kernel()
    driver = PySCFDriver(atom=atom, basis='sto3g')
    driver.run()
    calc = getattr(driver, '_calc', None)
    H_sys, offset, _, e_hf, e_fci = hamiltonian(mol, R)
    e0 = dasar(H_sys) + offset
    tulis(nama, [dict(molekul=mol, R=R, qubit=H_sys.num_qubits, rhf_konvergen=int(bool(mf.converged)),
                      driver_konvergen=int(bool(calc.converged)) if calc is not None else '',
                      E_hf=e_hf, E_fci=e_fci, E_dasar=e0, selisih=abs(e0 - e_fci),
                      buang=int(abs(e0 - e_fci) > 1e-8 or not mf.converged
                                or (calc is not None and not calc.converged)))])
    print(nama, abs(e0 - e_fci), mf.converged)


# ---- task list (same order as src/ujung.py) ---------------------------------------------------------------------
TUGAS = []
for mol, Rs in UJUNG.items():
    for R in Rs:
        TUGAS += [('ujung', mol, R, b, na, L) for b, na in HAA for L in range(1, batas(b, mol) + 1)]
        TUGAS += [('ujung', mol, R, b, 0, L) for b in HEA for L in range(1, batas(b, mol) + 1)]
TUGAS = [t + (k,) for t in TUGAS for k in range(RESTART)]
TUGAS += [('cek', mol, R) for mol in MOL for R in daftar(mol)]

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print(len(TUGAS))
    elif sys.argv[1] == 'list':
        for i, t in enumerate(TUGAS, 1):
            print(i, *t)
    else:
        t = TUGAS[int(sys.argv[1]) - 1]
        cek(*t[1:]) if t[0] == 'cek' else pec(*t)
