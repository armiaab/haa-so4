import numpy as np
import scipy.sparse.linalg as sla
from qiskit.quantum_info import Statevector
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.transformers import FreezeCoreTransformer

import haa


def ganjil_y(H):
    return sum(1 for lab, c in zip(H.paulis.to_labels(), H.coeffs) if abs(c) > 1e-12 and lab.count('Y') % 2)


def k_energi(m, psi, e0, ambang):
    urut = np.argsort(-np.abs(psi))
    for k in range(1, len(psi) + 1):
        v = np.zeros_like(psi)
        v[urut[:k]] = psi[urut[:k]]
        v /= np.linalg.norm(v)
        if np.real(np.vdot(v, m @ v)) - e0 < ambang:
            return k


sistem, konfigurasi, simetri = [], [], []
for mol in ('H2', 'H4', 'LiH', 'HF', 'BeH2'):
    penuh = PySCFDriver(atom=haa.GEOMETRI[mol](haa.SETIMBANG[mol]), basis='sto3g').run()
    aktif = FreezeCoreTransformer().transform(penuh)
    z = aktif.get_tapered_mapper(JordanWignerMapper()).z2symmetries
    simetri.append(dict(molekul=mol, R=haa.SETIMBANG[mol], simetri=' '.join(s.to_label() for s in z.symmetries),
                        qubit=' '.join(map(str, z.sq_list)), sektor=' '.join(map(str, z.tapering_values))))
    for R in haa.JARAK[mol]:
        H, offset, hf, e_hf, e_fci = haa.hamiltonian(mol, R)
        Hjw = haa.hamiltonian(mol, R, 'jw')[0]
        m = H.to_matrix(sparse=True)
        mjw = Hjw.to_matrix(sparse=True)
        w, v = np.linalg.eigh(m.toarray()) if m.shape[0] <= 4 else sla.eigsh(m, k=1, which='SA', maxiter=100000)
        e0, psi = w[0], v[:, 0]
        hf_sv = np.abs(Statevector(hf).data)
        i_hf = int(np.argmax(hf_sv))
        psi = np.real(psi * np.conj(psi[i_hf]) / abs(psi[i_hf]))
        psi /= np.linalg.norm(psi)
        p = psi ** 2
        ns = H.num_qubits
        sistem.append(dict(molekul=mol, orbital=penuh.num_spatial_orbitals, aktif=aktif.num_spatial_orbitals,
                           beku=penuh.num_spatial_orbitals - aktif.num_spatial_orbitals, elektron_alfa=aktif.num_particles[0],
                           elektron_beta=aktif.num_particles[1], qubit_jw=2 * aktif.num_spatial_orbitals, R=R, ns=ns,
                           pauli=len(H), maks_imag=np.abs(H.coeffs.imag).max(), E_hf=e_hf, E_fci=e_fci, E0=e0 + offset,
                           korelasi_mHa=(e_hf - e_fci) * 1000, w_hf=p[i_hf], PR=1 / np.sum(p ** 2),
                           n99=int(np.searchsorted(np.cumsum(np.sort(p)[::-1]), 0.99) + 1),
                           k_ca=k_energi(m, psi, e0, haa.CHEM), k_01=k_energi(m, psi, e0, 1e-4), hf_bit=format(i_hf, f'0{ns}b'),
                           n_jw=Hjw.num_qubits, n_taper=ns, suku_jw=len(Hjw), suku_taper=len(H), ganjilY_jw=ganjil_y(Hjw),
                           ganjilY_taper=ganjil_y(H), imag_jw=np.abs(mjw.data.imag).max(), imag_taper=np.abs(m.data.imag).max(),
                           hf_basis=int(np.isclose(hf_sv.max(), 1))))
        for urutan, i in enumerate(np.argsort(-p)[:8]):
            konfigurasi.append(dict(molekul=mol, R=R, urutan=urutan, bit=format(int(i), f'0{ns}b'), c=psi[i],
                                    balik=' '.join(str(q) for q in range(ns) if ((int(i) ^ i_hf) >> q) & 1)))
haa.tulis('sistem.csv', sistem)
haa.tulis('sistem_konfigurasi.csv', konfigurasi)
haa.tulis('simetri.csv', simetri)
print(len(sistem))
