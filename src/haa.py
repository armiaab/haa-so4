import csv
import os

import numpy as np
import scipy.sparse.linalg as sla
from pyscf import gto, mcscf, scf
from qiskit import QuantumCircuit
from qiskit.circuit import ParameterVector
from qiskit.circuit.library import CXGate
from qiskit.primitives import StatevectorEstimator
from qiskit.quantum_info import Operator, SparsePauliOp, Statevector
from qiskit.synthesis import TwoQubitBasisDecomposer
from qiskit_algorithms.gradients import ReverseEstimatorGradient
from qiskit_nature.second_q.circuit.library import HartreeFock
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.transformers import FreezeCoreTransformer

GEOMETRI = {
    'H2': lambda R: f'H 0 0 0; H 0 0 {R}',
    'H4': lambda R: '; '.join(f'H 0 0 {i * R:.6f}' for i in range(4)),
    'LiH': lambda R: f'Li 0 0 0; H 0 0 {R}',
    'HF': lambda R: f'H 0 0 0; F 0 0 {R}',
    'BeH2': lambda R: f'Be 0 0 0; H 0 0 {R}; H 0 0 {-R}',
}
JARAK = {
    'H2': [0.5, 0.735, 1.0, 1.5, 2.0, 2.5],
    'H4': [0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5],
    'LiH': [1.0, 1.3, 1.595, 2.0, 2.6, 3.2],
    'HF': [0.7, 0.8, 0.9, 1.0, 1.1, 1.3, 1.6],
    'BeH2': [1.0, 1.2, 1.326, 1.5, 1.7, 2.0],
}
SETIMBANG = {'H2': 0.735, 'H4': 0.9, 'LiH': 1.595, 'HF': 1.0, 'BeH2': 1.326}
BLOK = {'can': ('ZZ', 'XX', 'YY'), 'zeng': ('ZZ', 'XX', 'YY'), 'real': ('XY', 'YX'), 'so4': ('XY', 'YX', 'YZ', 'ZY'),
        'canu3': ('ZZ', 'XX', 'YY'), 'so4u3': ('XY', 'YX', 'YZ', 'ZY')}
SATU = {'hea': ('ry', 'rz'), 'ra': ('ry',), 'heau3': ('rz', 'ry', 'rz')}
GENAP = ('rz', 'ZZ', 'XX', 'YY')
CHEM = 1.5936e-3
HASIL = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results'))
decom = TwoQubitBasisDecomposer(CXGate())


def baca(nama):
    with open(os.path.join(HASIL, nama)) as f:
        return list(csv.DictReader(f))


def tulis(nama, baris):
    p = os.path.join(HASIL, nama)
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


def dkl(u, bin=25):
    q = np.histogram(u, bins=bin, range=(0, 1))[0] / len(u)
    q = q[q > 0]
    return float(np.sum(q * np.log(q * bin)))


def selang(u, rng):
    d = dkl(u)
    bawah, atas = np.percentile([dkl(rng.choice(u, len(u))) for _ in range(200)], [2.5, 97.5])
    return d, max(0.0, 2 * d - atas), 2 * d - bawah


def peringkat(M, skala):
    return int((np.linalg.svd(M, compute_uv=False) > 1e-7 * skala).sum())


def arah(qc, x):
    psi = Statevector(qc.assign_parameters(x)).data
    J = np.array([0.5 * Statevector(qc.assign_parameters(x + np.pi * e)).data for e in np.eye(len(x))])
    a, b = psi.real, psi.imag
    w, v = np.linalg.eigh([[a @ a, a @ b], [a @ b, b @ b]])
    fase = np.exp(1j * np.arctan2(-v[1, 1], v[0, 1]))
    u, Ju = psi * fase, J * fase
    r = u.real / np.linalg.norm(u.real)
    Re = Ju.real - np.outer(Ju.real @ r, r)
    Im = Ju.imag - np.outer(Ju.imag @ r, r)
    semua = np.hstack([Re, Im])
    skala = np.linalg.svd(semua, compute_uv=False)[0]
    nR, nI = np.linalg.norm(Re, axis=1), np.linalg.norm(Im, axis=1)
    return dict(bobot_imajiner=w[0] / (w[0] + w[1]), D_R=peringkat(Re, skala), D_I=peringkat(Im, skala),
                D_C=peringkat(semua, skala), hidup_R=int((nR > 1e-6).sum()), hidup_I_saja=int(((nR <= 1e-6) & (nI > 1e-6)).sum()))
