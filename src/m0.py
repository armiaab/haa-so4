from functools import reduce

import numpy as np
from qiskit.quantum_info import Statevector
from qiskit_nature.second_q.circuit.library import HartreeFock
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.transformers import FreezeCoreTransformer
from scipy.optimize import minimize

import haa

MOL = ('H4', 'LiH', 'HF', 'BeH2')

lantai = []
for mol in MOL:
    for R in haa.JARAK[mol]:
        H_sys, offset, hf, e_hf, e_fci = haa.hamiltonian(mol, R)
        ns = H_sys.num_qubits
        H = H_sys.to_matrix(sparse=True)
        b = int(np.argmax(np.abs(Statevector(hf).data)))

        def psi(t):
            v = [np.array([-np.sin(t[q] / 2), np.cos(t[q] / 2)]) if (b >> q) & 1 else np.array([np.cos(t[q] / 2), np.sin(t[q] / 2)])
                 for q in range(ns)]
            return reduce(np.kron, v[::-1])

        E = lambda t: float(np.real(psi(t) @ (H @ psi(t)))) + offset
        x0s = [np.zeros(ns)] + [np.random.default_rng(k).uniform(-np.pi, np.pi, ns) for k in range(40)]
        e = min(minimize(E, x0, method='L-BFGS-B').fun for x0 in x0s)
        lantai.append(dict(molekul=mol, R=R, galat_hf_mHa=1e3 * (e_hf - e_fci), galat_lantai_mHa=1e3 * (e - e_fci),
                           sisa_korelasi=(e - e_fci) / (e_hf - e_fci)))
haa.tulis('m0_lantai.csv', lantai)

flip = []
for mol in MOL:
    R = haa.SETIMBANG[mol]
    problem = FreezeCoreTransformer().transform(PySCFDriver(atom=haa.GEOMETRI[mol](R), basis='sto3g').run())
    jw = JordanWignerMapper()
    tm = problem.get_tapered_mapper(jw)
    op = problem.hamiltonian.second_q_op()
    Hf, Ht, z2 = jw.map(op).to_matrix(sparse=True), tm.map(op).to_matrix(sparse=True), tm.z2symmetries
    n, sq, val = Hf.shape[0].bit_length() - 1, list(z2.sq_list), list(z2.tapering_values)
    sisa = [q for q in range(n) if q not in sq]
    C = [c.to_matrix(sparse=True) for c in z2.cliffords]
    hf_f = int(np.argmax(np.abs(Statevector(HartreeFock(problem.num_spatial_orbitals, problem.num_particles, jw)).data)))
    hf_t = int(np.argmax(np.abs(Statevector(HartreeFock(problem.num_spatial_orbitals, problem.num_particles, tm)).data)))
    N = bin(hf_f).count('1')

    def penuh(b):
        v = np.zeros(2 ** n, complex)
        for pola in range(2 ** len(sq)):
            x, a = 0, 1.0
            for j, q in enumerate(sisa):
                x |= ((b >> j) & 1) << q
            for i, q in enumerate(sq):
                x |= ((pola >> i) & 1) << q
                a *= (val[i] if (pola >> i) & 1 else 1) / np.sqrt(2)
            v[x] += a
        for c in C:
            v = c @ v
        return v

    for j in range(len(sisa)):
        phi = penuh(hf_t ^ (1 << j))
        bobot = {}
        for x in np.nonzero(np.abs(phi) > 1e-12)[0]:
            ne = bin(int(x)).count('1')
            k = f'eksitasi{bin(int(x) ^ hf_f).count("1") // 2}' if ne == N else 'elektron_lain'
            bobot[k] = bobot.get(k, 0) + abs(phi[x]) ** 2
        flip.append(dict(molekul=mol, R=R, qubit=j, elemen=float(np.real(Ht[hf_t ^ (1 << j), hf_t])), **bobot))
haa.tulis('flip_taper.csv', flip)
print(len(lantai), len(flip))
