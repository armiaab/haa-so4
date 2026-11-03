"""Collect the ~28,000 result CSVs in results/ into a few files for the data deposit.

  vqe_runs.csv                 every VQE run (results/pec, results/uji, results/beh2), one row per start
  other_results.csv            every other table and analysis, one row per original row
  optimization_histories.csv   results/riwayat.csv (energy at every function evaluation)
  h4_traces.csv                results/h4_jejak.csv (H4 optimization traces)
  h4_energy_slices.csv         results/h4_potong.csv (H4 energy slices)
  manifest.csv                 original file -> collected file and original header (used by pulihkan.py)
  HAA-SO4_results.xlsx         the same tables for reading, one sheet per table, with a README sheet

Each collected CSV starts with a `source` column, the path of the original file relative to results/.
results/hasil_ujung (a copy of the cluster output already in results/uji) and results/shot.csv (seed bug, not
used in the paper) are left out. The .xlsx needs openpyxl; the CSVs need only the standard library.

usage: python src/kumpulkan_hasil.py <output dir>
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

HASIL = Path(__file__).resolve().parents[1] / 'results'
BESAR = {'riwayat.csv': 'optimization_histories.csv', 'h4_jejak.csv': 'h4_traces.csv',
         'h4_potong.csv': 'h4_energy_slices.csv'}
LEWATI = ('hasil_ujung/', 'shot.csv')
KETERANGAN = {
    'vqe_runs': 'Every VQE run: molecule, bond length R, design (blok, na = ancillas), layers L, CNOTs, starts '
                '(restart = seed), optimizer counts, energies, error (galat) in hartree.',
    'uji_awal': 'H4 starting-point tests (start type, kick size eps, optimizer).',
    'uji_cek': 'Rule-7 checks for every registered geometry (RHF convergence, tapered ground state vs FCI).',
    'derau': 'Density-matrix noise simulations (kappa, noise allowance).',
    'derau_sudut': 'Optimized angles used in the noise simulations.',
    'ekspresi': 'Expressibility (KL divergence) of each design.',
    'fasa': 'Ranks of the real and imaginary tangent directions.',
    'fasa_sudut': 'Angles of the direction-rank runs.',
    'flip_taper': 'What each tapered-qubit flip encodes (excitations).',
    'gradien': 'Gradient variance by parameter class.',
    'h4_ekspresi': 'H4: expressibility.', 'h4_gradien': 'H4: gradient samples.',
    'h4_hessian': 'H4: Hessian eigenvalues at the end of trapped runs.', 'h4_kurva': 'H4: potential curve and Hamiltonian data.',
    'h4_optim': 'H4: optimizations by start type.', 'h4_sudut': 'H4: final angles.', 'h4_varians': 'H4: gradient variance.',
    'intro_gradien': 'Data for the introduction figure (gradients).', 'intro_jejak': 'Data for the introduction figure (traces).',
    'intro_potong': 'Data for the introduction figure (energy slices).', 'korolari': 'Numerical check of the corollary.',
    'm0_lantai': 'Lowest energy over real product states (M0 floor).', 'nol': 'Gradients at theta = 0 for every design.',
    'riwayat_ringkas': 'Summary of the optimization histories.', 'simetri': 'Z2 symmetries used for tapering.',
    'sistem': 'Molecular descriptors (orbitals, qubits, energies).', 'sistem_konfigurasi': 'Leading configurations.',
    'versi': 'Software versions.',
}
GLOSARIUM = ('galat = error; lapis = layer; blok = block (zeng = published HAA-CAN, can = CAN in R_Y layers, so4, '
             'real, hea, heau3, ra = RealAmplitudes); na = number of ancillas (0 = brick, 1 = star); sudut = angle; '
             'jarak = distance/range; setimbang = equilibrium; regang = stretched; derau = noise; diam = idle; '
             'awal = start; acak = random; nol = zero; lantai = floor; hemat = saving; biaya = cost; kedalaman = depth; '
             'rute = routing; uji = test; kontrol = control; terbaik = best; sukses = success; impas = break-even; '
             'peta = qubit mapping; eksitasi = excitation; riwayat = history; jejak = trace; potong = slice; '
             'ekspresi = expressibility; fasa = phase; pelengkap = fill-in; kena_batas = hit the iteration limit; '
             'rentang = range (asli = original, luas = extended, daftar = registered); buang = removed.')


def kelompok(rel):
    if rel.startswith(LEWATI) or rel in LEWATI:
        return None
    if rel in BESAR:
        return BESAR[rel]
    bagian = rel.split('/')
    if bagian[0] in ('pec', 'beh2') or (bagian[0] == 'uji' and not bagian[1].startswith(('awal_', 'cek_'))):
        return 'vqe_runs.csv'
    return 'other_results.csv'


def lembar(rel):
    """Sheet of the .xlsx for an original file."""
    bagian = rel[:-4].split('/')
    if bagian[0] in ('pec', 'beh2') or (bagian[0] == 'uji' and not bagian[1].startswith(('awal_', 'cek_'))):
        return 'vqe_runs'
    if bagian[0] == 'uji':
        return 'uji_' + bagian[1].split('_')[0]
    return '_'.join(bagian)[:31]


def angka(v):
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return v


def main():
    keluar = Path(sys.argv[1])
    keluar.mkdir(parents=True, exist_ok=True)
    data, kolom, manifest = defaultdict(list), defaultdict(dict), []
    for p in sorted(HASIL.rglob('*.csv')):
        rel = p.relative_to(HASIL).as_posix()
        k = kelompok(rel)
        if k is None:
            continue
        with open(p, newline='', encoding='utf-8') as f:
            r = csv.reader(f)
            kepala = next(r)
            baris = list(r)
        manifest.append((rel, k, ','.join(kepala)))
        kolom[k].update(dict.fromkeys(kepala))
        data[k] += [(rel, dict(zip(kepala, b))) for b in baris]
    for k, isi in data.items():
        with open(keluar / k, 'w', newline='', encoding='utf-8') as f:
            w = csv.writer(f)
            w.writerow(['source'] + list(kolom[k]))
            w.writerows([rel] + [b.get(c, '') for c in kolom[k]] for rel, b in isi)
        print(f'{k}: {len(isi)} rows')
    with open(keluar / 'manifest.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['source', 'collected_file', 'header'])
        w.writerows(manifest)
    print(f'manifest.csv: {len(manifest)} files')
    if '--tanpa-xlsx' not in sys.argv:
        tulis_xlsx(keluar / 'HAA-SO4_results.xlsx', data)


def tulis_xlsx(nama, data):
    from openpyxl import Workbook
    from openpyxl.styles import Font
    lembaran = defaultdict(list)
    for k in ('vqe_runs.csv', 'other_results.csv'):
        for rel, b in data[k]:
            lembaran[lembar(rel)].append((rel, b))
    wb = Workbook()
    ws = wb.active
    ws.title = 'README'
    ws.append(['HAA-SO(4) study: results in one workbook. Large per-iteration data are in the CSV files '
               '(optimization_histories.csv, h4_traces.csv, h4_energy_slices.csv). Energies in hartree unless a '
               'column name ends in mHa; R in angstrom.'])
    ws.append(['Column names are Indonesian: ' + GLOSARIUM])
    ws.append([])
    ws.append(['sheet', 'rows', 'source files', 'description'])
    for c in ws[4]:
        c.font = Font(bold=True)
    for nama_l, isi in lembaran.items():
        sumber = sorted({rel for rel, _ in isi})
        dasar = nama_l.removeprefix('tabel_')
        ket = KETERANGAN.get(nama_l) or KETERANGAN.get(dasar) or ('Summary table of the paper.' if nama_l.startswith('tabel_') else '')
        ws.append([nama_l, len(isi), sumber[0] if len(sumber) == 1 else f'{len(sumber)} files ({sumber[0]}, ...)', ket])
        kol = list(dict.fromkeys(c for _, b in isi for c in b))
        lb = wb.create_sheet(nama_l)
        lb.append(['source'] + kol if len(sumber) > 1 else kol)
        for c in lb[1]:
            c.font = Font(bold=True)
        for rel, b in isi:
            nilai = [angka(b.get(c, '')) if b.get(c, '') != '' else None for c in kol]
            lb.append([rel] + nilai if len(sumber) > 1 else nilai)
        lb.freeze_panes = 'A2'
    wb.save(nama)
    print(f'{nama.name}: {len(lembaran)} sheets')


if __name__ == '__main__':
    main()
