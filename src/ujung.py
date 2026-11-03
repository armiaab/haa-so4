"""Exact registered endpoints (rule 2: up to 2.5 R0) and the rule-7 geometry checks.

The extended grids of uji.py end at rounded endpoints for two molecules (LiH 4.0 A instead of 3.9875 A, BeH2
3.3 A instead of 3.315 A). This job adds the exact endpoints for every design that enters tab:lstar, with the
registered caps (rule 5) and 12 starts, and records for every geometry of the registered range whether RHF
converged and how far the lowest eigenvalue of the tapered Hamiltonian lies from E_FCI (rule 7).

Cluster:   cd src && sbatch ujung.sh          (array job; one task = one start of one configuration, or one check)
           python ujung.py                    (prints the number of tasks)
           python ujung.py cek                (summarizes the checks into results/tabel/ujung_cek.csv)
           python uji.py tabel                (adds the 'daftar' range to results/tabel/uji_lstar.csv)
Local:     copy results/uji/ujung_*.csv, results/uji/cek_*.csv and results/tabel/*.csv back, then
           cd src && python angka.py ../angka.tex
"""
import glob
import os
import sys

from pyscf import gto, scf
from qiskit_nature.second_q.drivers import PySCFDriver

import haa
import uji

TUGAS = []
for mol, Rs in uji.UJUNG.items():
    for R in Rs:
        TUGAS += [('ujung', mol, R, b, na, L) for b, na in uji.HAA for L in range(1, uji.batas(b, mol) + 1)]
        TUGAS += [('ujung', mol, R, b, 0, L) for b in uji.HEA for L in range(1, uji.batas(b, mol) + 1)]
TUGAS = [t + (k,) for t in TUGAS for k in range(uji.RESTART)]
TUGAS += [('cek', mol, R) for mol in uji.MOL for R in uji.daftar(mol)]


def cek(mol, R):
    nama = f'uji/cek_{mol}_R{R:.4f}.csv'
    if os.path.exists(os.path.join(haa.HASIL, nama)):
        return
    atom = haa.GEOMETRI[mol](R)
    mf = scf.RHF(gto.M(atom=atom, basis='sto3g', verbose=0))
    mf.kernel()
    driver = PySCFDriver(atom=atom, basis='sto3g')
    driver.run()
    calc = getattr(driver, '_calc', None)
    H_sys, offset, _, e_hf, e_fci = haa.hamiltonian(mol, R)
    e0 = haa.dasar(H_sys) + offset
    haa.tulis(nama, [dict(molekul=mol, R=R, qubit=H_sys.num_qubits, rhf_konvergen=int(bool(mf.converged)),
                          driver_konvergen=int(bool(calc.converged)) if calc is not None else '',
                          E_hf=e_hf, E_fci=e_fci, E_dasar=e0, selisih=abs(e0 - e_fci),
                          buang=int(abs(e0 - e_fci) > 1e-8 or not mf.converged
                                    or (calc is not None and not calc.converged)))])
    print(nama, abs(e0 - e_fci), mf.converged)


def ringkas():
    baris = [b for f in sorted(glob.glob(os.path.join(haa.HASIL, 'uji', 'cek_*.csv'))) for b in haa.baca(f)]
    haa.tulis('tabel/ujung_cek.csv', baris)
    for b in baris:
        print(b['molekul'], b['R'], 'RHF', b['rhf_konvergen'], 'driver', b['driver_konvergen'], 'selisih', b['selisih'],
              'buang', b['buang'])
    print(len(baris), 'geometries;', sum(int(b['buang']) for b in baris), 'to remove under rule 7')


def jalan(i):
    t = TUGAS[i - 1]
    cek(*t[1:]) if t[0] == 'cek' else uji.pec(*t)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print(len(TUGAS))
    elif sys.argv[1] == 'cek':
        ringkas()
    else:
        jalan(int(sys.argv[1]))
