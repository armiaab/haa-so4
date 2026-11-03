import os
import time
from multiprocessing import Pool

import numpy as np
from qiskit import transpile
from qiskit.transpiler import CouplingMap

import haa

TUGAS = [(b['molekul'], b['topologi'], b['blok'], int(b['L']), int(b['na']), int(b['qubit']) - int(b['na']))
         for b in haa.baca('tabel/tab2_biaya.csv')]


def jalan(tugas):
    mol, topologi, ansatz, L, na, ns = tugas
    t0 = time.perf_counter()
    sudut = np.random.default_rng(1).uniform(-np.pi, np.pi, len(haa.label(ansatz, ns, na, L)))
    qc = haa.ansatz(ansatz, ns, na, L, sudut=sudut, kompilasi=True)
    n = qc.num_qubits
    baris = dict(molekul=mol, topologi=topologi, blok=ansatz, L=L, qubit=n, N2q_rumus=haa.biaya(ansatz, ns, na, L)[0])
    for nama, peta in (('penuh', None), ('garis', CouplingMap.from_line(n)), ('heavyhex', CouplingMap.from_heavy_hex(3))):
        hasil = [transpile(qc, basis_gates=['cx', 'u'], coupling_map=peta, optimization_level=3, seed_transpiler=s)
                 for s in range(8)]
        baris[f'N2q_{nama}'], baris[f'D2q_{nama}'] = min((t.count_ops().get('cx', 0),
                                                          t.depth(lambda i: i.operation.num_qubits == 2)) for t in hasil)
    baris['detik'] = round(time.perf_counter() - t0, 1)
    return baris


if __name__ == '__main__':
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', os.cpu_count()))) as pool:
        haa.tulis('tabel/rute.csv', pool.map(jalan, TUGAS))
    print(len(TUGAS))
