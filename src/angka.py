import glob
import math
import os
import statistics
import sys

import paper_data as haa

MOL = ('H4', 'LiH', 'HF', 'BeH2')
keluar = sys.argv[1] if len(sys.argv) > 1 else os.path.join(haa.HASIL, 'angka.tex')
nilai = {}


def ilmiah(x, d=1):
    mantisa, eksponen = f'{x:.{d}e}'.split('e')
    return rf'{mantisa}\times10^{{{int(eksponen)}}}'


def dua(x):
    return f'{x:.{max(0, 1 - math.floor(math.log10(abs(x))))}f}'


def hitung(baris):
    return f"{sum(float(b['galat']) < haa.CHEM for b in baris)}/{len(baris)}"


t = {(b['molekul'], b['topologi'], b['blok']): b for b in haa.baca('tabel/tab2_biaya.csv')}
for m in MOL:
    for top in ('star', 'brick'):
        for blok, nama in (('zeng', 'can'), ('so4', 'so4')):
            nilai[f'N-{m}-{nama}-{top}'] = t[(m, top, blok)]['N2q']
            nilai[f'L-{m}-{nama}-{top}'] = t[(m, top, blok)]['L']
    nilai[f'N-{m}-hea'] = t[(m, 'brick', 'hea')]['N2q']
    nilai[f'L-{m}-hea'] = t[(m, 'brick', 'hea')]['L']
    can = min(int(t[(m, top, 'zeng')]['N2q']) for top in ('star', 'brick'))
    nilai[f'N-{m}-can-terbaik'] = can
    nilai[f'hemat-{m}'] = round(100 * (1 - int(t[(m, 'brick', 'so4')]['N2q']) / int(t[(m, 'brick', 'zeng')]['N2q'])))
hemat = [nilai[f'hemat-{m}'] for m in MOL]
nilai['hemat-min'], nilai['hemat-maks'] = min(hemat), max(hemat)
nilai['n-can-lebih-hea'] = sum(nilai[f'N-{m}-can-terbaik'] > int(nilai[f'N-{m}-hea']) for m in MOL)
nilai['n-can-lapis-kurang'] = sum(min(int(t[(m, top, 'zeng')]['L']) for top in ('star', 'brick')) < int(nilai[f'L-{m}-hea']) for m in MOL)
nilai['n-so4-brick-le-hea'] = sum(int(nilai[f'N-{m}-so4-brick']) <= int(nilai[f'N-{m}-hea']) for m in MOL)
banding = [int(t[(m, top, 'so4')]['L']) - int(t[(m, top, b)]['L']) for m in MOL for top in ('star', 'brick') for b in ('can', 'zeng')]
nilai['n-banding-lapis'] = len(banding)
nilai['n-so4-lapis-kurang'] = sum(x < 0 for x in banding)
nilai['n-so4-lapis-sama'] = sum(x == 0 for x in banding)
nilai['n-so4-lapis-lebih'] = sum(x > 0 for x in banding)
nol = [b for b in haa.baca('nol.csv') if b['molekul'] == 'H4' and b['desain'] in ('can-star', 'can-brick', 'so4-star')]
nilai['nol-H4-maks'] = ilmiah(max(float(b['grad_maks']) for b in nol))
sektor = max(abs(float(b['E_dasar']) - float(b['E_fci'])) for f in glob.glob(os.path.join(haa.HASIL, 'pec', '*_taper_*.csv'))
             for b in haa.baca(f))
nilai['sektor-maks'] = ilmiah(sektor)
nilai['n-pec-taper'] = len(glob.glob(os.path.join(haa.HASIL, 'pec', '*_taper_*.csv')))

# Results chapter. Block comparisons (La-, Na-) stay on the original ranges; uji.py adds no geometries for them.
NAMA = {'can': 'kontrol', 'zeng': 'can', 'real': 'real', 'so4': 'so4', 'hea': 'hea'}
sel = [(m, top) for m in MOL for top in ('star', 'brick')]
for (m, top, blok), b in t.items():
    nama = NAMA[blok] if blok == 'hea' else f'{NAMA[blok]}-{top}'
    nilai[f'La-{m}-{nama}'], nilai[f'Na-{m}-{nama}'], nilai[f'D-{m}-{nama}'] = b['L'], b['N2q'], b['D2q']
    nilai[f'sukses-{m}-{nama}'] = round(100 * float(b['sukses']))


def lapis(m, top, blok):
    return int(t[(m, top, blok)]['L'])


def cnot(m, top, blok):
    return int(t[(m, top, blok)]['N2q'])


nilai['n-sel'] = len(sel)
nilai['n-kontrol-lapis-lebih'] = sum(lapis(m, top, 'can') > lapis(m, top, 'zeng') for m, top in sel)
for blok in ('can', 'zeng'):
    nilai[f'n-so4-{NAMA[blok]}-kurang'] = sum(lapis(m, top, 'so4') < lapis(m, top, blok) for m, top in sel)
    nilai[f'n-so4-{NAMA[blok]}-sama'] = sum(lapis(m, top, 'so4') == lapis(m, top, blok) for m, top in sel)
    nilai[f'n-so4-{NAMA[blok]}-murah'] = sum(cnot(m, top, 'so4') < cnot(m, top, blok) for m, top in sel)
nilai['n-real-lapis-lebih'] = sum(lapis(m, top, 'real') > lapis(m, top, 'so4') for m, top in sel)
nilai['n-real-lapis-sama'] = sum(lapis(m, top, 'real') == lapis(m, top, 'so4') for m, top in sel)
nilai['n-real-kontrol-murah'] = sum(cnot(m, top, 'real') < cnot(m, top, 'can') for m, top in sel)
nilai['n-hemat-sepertiga'] = sum(3 * cnot(m, top, 'so4') == 2 * cnot(m, top, 'zeng') for m, top in sel)
hemat = [round(100 * (1 - cnot(m, 'star', 'so4') / cnot(m, 'star', 'zeng'))) for m in MOL]
nilai['hemat-star-min'], nilai['hemat-star-maks'] = min(hemat), max(hemat)
nilai['n-so4-brick-dangkal'] = sum(int(t[(m, 'brick', 'so4')]['D2q']) < int(t[(m, 'brick', 'hea')]['D2q']) for m in MOL)
sukses = [nilai[f'sukses-{m}-so4-{top}'] for m in MOL for top in ('star', 'brick')]
nilai['sukses-so4-min'] = min(sukses)
ketat = {(b['molekul'], b['topologi'], b['blok']): b for b in haa.baca('tabel/tab5_ketahanan.csv')}


def cnot_ketat(m, top, blok, ambang):
    x = ketat[(m, top, blok)][f'N2q_{ambang}']
    return int(x) if x else math.inf


for ambang, nama in (('1', '1'), ('0.1', '01')):
    nilai[f'n-ketat-{nama}'] = sum(cnot_ketat(m, top, 'so4', ambang) < cnot_ketat(m, top, 'zeng', ambang) for m, top in sel)
nilai['n-can-gagal-01'] = sum(not ketat[(m, top, 'zeng')]['N2q_0.1'] for m, top in sel)

impas = haa.baca('tabel/tab3_impas.csv')
nilai['n-impas'] = len(impas)
nilai['n-ancilla-cnot'] = sum(b['ancilla_menang_N2q'] == 'True' for b in impas)
nilai['n-ancilla-dalam'] = sum(b['ancilla_menang_D2q'] == 'True' for b in impas)
nilai['rasio-impas-maks'] = f"{max(float(b['rasio']) for b in impas):g}"
nilai['impas-dalam-min'] = f"{min(float(b['impas_D2q']) for b in impas):g}"
rute = haa.baca('tabel/rute.csv')
for b in rute:
    if b['topologi'] == 'star' and b['blok'] in ('zeng', 'so4'):
        nilai[f'Nhh-{b["molekul"]}-{NAMA[b["blok"]]}-star'] = b['N2q_heavyhex']
lurus = [b for b in rute if b['topologi'] == 'brick']
nilai['n-rute-brick'] = len(lurus)
nilai['n-rute-brick-sama'] = sum(b['N2q_heavyhex'] == b['N2q_garis'] == b['N2q_rumus'] for b in lurus)

derau = haa.baca('derau.csv')
for b in derau:
    kunci = f'{b["molekul"]}-{b["geometri"]}-{b["peran"].replace("_", "-")}'
    nilai[f'pn-{kunci}'] = dua(float(b['p_derau']) / 1e-5)
    nilai[f'ps-{kunci}'] = dua(float(b['p_bintang']) / 1e-5)
    nilai[f'opt-{kunci}'] = f'{float(b["galat_optimasi"]):.2f}'
kepekaan = [math.log(float(b['kepekaan'])) for b in derau]
log_cnot = [math.log(int(b['cnot'])) for b in derau]
log_dalam = [math.log(int(b['kedalaman2q'])) for b in derau]
nilai['n-derau'] = len(derau)
nilai['korelasi-N'] = f'{statistics.correlation(kepekaan, log_cnot):.2f}'
nilai['korelasi-D'] = f'{statistics.correlation(kepekaan, log_dalam):.2f}'
nilai['korelasi-ND'] = f'{statistics.correlation(log_cnot, log_dalam):.2f}'
turun = {i: 100 * (1 - float(b['p_derau_diam']) / float(b['p_derau'])) for i, b in enumerate(derau)}
nilai['diam-median'], nilai['diam-maks'] = round(statistics.median(turun.values())), round(max(turun.values()))
for top in ('star', 'brick'):
    x = [turun[i] for i, b in enumerate(derau) if b['topologi'] == top]
    nilai[f'diam-{top}-min'], nilai[f'diam-{top}-maks'] = round(min(x)), round(max(x))
PERAN = ('hea', 'can_terbaik', 'so4_star', 'so4_brick')
GEO = ('setimbang', 'regang')


def toleransi(model):
    p = {(b['molekul'], b['geometri'], b['peran']): float(b[f'p_derau{model}']) for b in derau}
    return p, {(m, g): sorted(PERAN, key=lambda r: p[(m, g, r)]) for m in MOL for g in GEO}


p_derau, urutan = toleransi('')
urutan_diam = toleransi('_diam')[1]
nilai['n-titik-derau'] = len(urutan)
nilai['n-urutan-diam-sama'] = sum(urutan[x] == urutan_diam[x] for x in urutan)
for model in ('', '_diam'):
    p = toleransi(model)[0]
    so4 = {(m, g): max(p[(m, g, 'so4_star')], p[(m, g, 'so4_brick')]) for m in MOL for g in GEO}
    rasio = [so4[x] / p[x + ('can_terbaik',)] for x in so4]
    kunci = 'rasio-pn-can' + model.replace('_', '-')
    nilai[f'{kunci}-min'], nilai[f'{kunci}-maks'] = f'{min(rasio):.1f}', f'{max(rasio):.1f}'
so4 = {(m, g): max(p_derau[(m, g, 'so4_star')], p_derau[(m, g, 'so4_brick')]) for m in MOL for g in GEO}
for m in MOL:
    rasio = [so4[(m, g)] / p_derau[(m, g, 'hea')] for g in GEO]
    nilai[f'rasio-pn-hea-{m}-min'], nilai[f'rasio-pn-hea-{m}-maks'] = f'{min(rasio):.1f}', f'{max(rasio):.1f}'


def gamma(p, d):
    return (1 + p * (d * d - 2) / d ** 2) / (1 - p)


for p, nama in ((3e-3, '3e-3'), (1e-2, '1e-2')):
    C = {(b['molekul'], b['geometri'], b['peran']): gamma(p, 4) ** (2 * int(b['cnot'])) *
         gamma(p / 10, 2) ** (2 * int(b['satu_qubit'])) for b in derau}
    nilai[f'pec-min-{nama}'], nilai[f'pec-maks-{nama}'] = f'{min(C.values()):.3g}', f'{max(C.values()):.3g}'
    for peran in PERAN:
        nilai[f'pec-BeH2-{peran.replace("_", "-")}-{nama}'] = f"{C[('BeH2', 'setimbang', peran)]:.3g}"

nol = haa.baca('nol.csv')
for m in MOL:
    for desain in ('can-star', 'can-brick', 'so4-star', 'so4-brick', 'hea'):
        baris = [b for b in nol if b['molekul'] == m and b['desain'] == desain]
        nilai[f'g0-{m}-{desain}'] = ilmiah(max(float(b['grad_maks']) for b in baris))
        nilai[f'stas-{m}-{desain}'] = f"{sum(int(b['stasioner']) for b in baris)}/{len(baris)}"
nilai['n-stas-campur'] = sum(nilai[k].split('/')[0] not in ('0', nilai[k].split('/')[1]) for k in nilai if k.startswith('stas-'))
nilai['nol-genap-maks'] = ilmiah(max(float(b['grad_genap_maks']) for b in nol), 0)
for b in haa.baca('flip_taper.csv'):
    if b['eksitasi2']:
        nilai[f'flip-{b["molekul"]}'] = ilmiah(abs(float(b['elemen'])), 2)
        nilai[f'flip-qubit-{b["molekul"]}'] = ' and '.join(filter(None, (nilai.get(f'flip-qubit-{b["molekul"]}'), b['qubit'])))
lantai = haa.baca('m0_lantai.csv')
for m in MOL:
    pulih = {float(b['R']): 100 * (1 - float(b['sisa_korelasi'])) for b in lantai if b['molekul'] == m}
    nilai[f'lantai-{m}-min'], nilai[f'lantai-{m}-maks'] = round(max(min(pulih.values()), 0)), round(max(pulih.values()))
    nilai[f'lantai-{m}-n'] = sum(x > 0.5 for x in pulih.values())
    nilai[f'lantai-{m}-R'] = f'{max(pulih, key=pulih.get):g}'
nilai['lantai-galat-min'] = f"{min(float(b['galat_lantai_mHa']) for b in lantai):.1f}"

optim = haa.baca('h4_optim.csv')
LSTAR = {'hea': 4, 'can-star': 2, 'can-brick': 3, 'so4-star': 2, 'so4-brick': 2}
for desain, L in LSTAR.items():
    for jenis in ('m0', 'acak', 'nol'):
        nilai[f'awal-{desain}-{jenis}'] = hitung([b for b in optim if b['desain'] == desain and int(b['L']) == L
                                                 and b['awal'] == jenis])
    for eps, nama in (('1e-08', 'e8'), ('0.0001', 'e4'), ('0.01', 'e2')):
        baris = [b for b in optim if b['desain'] == desain and b['awal'] == 'lolos' and b['eps'] == eps]
        if baris:
            nilai[f'tendang-{desain}-{nama}'] = hitung(baris)
            nilai[f'lepas-{desain}-{nama}'] = f"{sum(float(b['sudut_genap']) > 1e-3 for b in baris)}/{len(baris)}"
for eps, nama in (('1e-08', 'e8'), ('0.0001', 'e4'), ('0.01', 'e2')):
    baris = [b for b in optim if b['desain'].startswith('can') and b['awal'] == 'lolos' and b['eps'] == eps]
    nilai[f'tendang-can-{nama}'] = hitung(baris)
    nilai[f'lepas-can-{nama}'] = f"{sum(float(b['sudut_genap']) > 1e-3 for b in baris)}/{len(baris)}"
m0 = [b for b in optim if b['desain'].startswith('can') and b['awal'] == 'm0' and int(b['L']) == LSTAR[b['desain']]]
nilai['m0-sudut-genap-maks'] = ilmiah(max(float(b['sudut_genap']) for b in m0), 0)
for desain in ('can-star', 'can-brick'):
    for L in range(1, 7):
        nilai[f'bocor-{desain}-L{L}'] = hitung([b for b in optim if b['desain'] == desain and b['awal'] == 'm0'
                                               and int(b['L']) == L and b['R'] == '1.5'])
bocor = [b for b in optim if b['desain'].startswith('can') and b['awal'] == 'm0' and float(b['sudut_genap']) > 1e-3]
nilai['bocor-nit-min'], nilai['bocor-nit-maks'] = min(int(b['nit']) for b in bocor), max(int(b['nit']) for b in bocor)
nilai['m0-nit-maks'] = max(int(b['nit']) for b in m0)
eig = {}
for b in haa.baca('h4_hessian.csv'):
    if b['bagian'] == 'eig_genap' and b['awal'] == 'm0' and b['desain'].startswith('can'):
        k = (b['desain'], b['R'], b['benih'])
        eig[k] = min(eig.get(k, math.inf), float(b['eig']))
sadel = [x for x in eig.values() if x < -1e-2]
nilai['n-sadel'], nilai['n-hessian'] = len(sadel), len(eig)
nilai['lambda-min'], nilai['lambda-maks'] = f'{min(sadel):.2f}', f'{max(sadel):.2f}'
nilai['lebar-min'], nilai['lebar-maks'] = ilmiah(1e-5 / abs(min(sadel)), 0), ilmiah(1e-5 / abs(max(sadel)), 0)

for b in haa.baca('tabel/beh2_lstar.csv'):
    if b['L_bintang'] == 'kurang':
        continue
    k = f'beh2-{b["faktor"]}-{b["blok"]}-na{b["na"]}'
    nilai[k] = b['L_bintang'].replace('>', r'\ensuremath{>}')
    nilai[f'{k}-nR'] = b['n_R']
    if b['L_bintang'].isdigit():
        L, ns, na = int(b['L_bintang']), 7, int(b['na'])
        c = {'so4': 2, 'zeng': 3}.get(b['blok'], 1)
        nilai[f'{k}-N'] = c * (ns if na else ns - 1) * L
        nilai[f'{k}-D'] = c * (ns if na else 2) * L if b['blok'] in ('so4', 'zeng') else ns + 2 * L - 3
# Decision-rule runs (uji.py tabel). tab:lstar and every HAA-CAN or HEA comparison use the range of rule 2 ('luas');
# the block controls stay on the original ranges, for which uji.py adds no geometries.
uji = {(b['molekul'], b['blok'], int(b['na']), b['rentang']): b for b in haa.baca('tabel/uji_lstar.csv')}
ns = {m: int(t[(m, 'brick', 'so4')]['qubit']) for m in MOL}
BATAS = {'zeng': 6, 'so4': 6, 'hea': 10, 'heau3': 10, 'ra': 10}
DESAIN = {'can-star': ('zeng', 1), 'can-brick': ('zeng', 0), 'so4-star': ('so4', 1), 'so4-brick': ('so4', 0),
          'hea': ('hea', 0), 'heau3': ('heau3', 0), 'ra': ('ra', 0)}
# Main range: the registered rule-2 range ('daftar', exact 2.5 R0 endpoints, rule-7 geometries removed) once every
# design has complete data there; until then the sampled extended grids ('luas').
UTAMA = 'daftar' if all(uji.get((m, b, na, 'daftar')) and 'kurang' not in (uji[(m, b, na, 'daftar')]['L_terbaik'],
                                                                       uji[(m, b, na, 'daftar')]['L_tiga'])
                        for m in MOL for b, na in DESAIN.values()) else 'luas'
print('main range for tab:lstar and the HEA/HAA-CAN comparisons:', UTAMA)
nilai['rentang-daftar'] = int(UTAMA == 'daftar')


def bintang(m, blok, na, krit='terbaik', rentang=UTAMA):
    b = uji.get((m, blok, na, rentang)) or uji[(m, blok, na, 'asli')]
    L = b[f'L_{krit}']
    if L.isdigit():
        return int(L), int(b[f'cnot_{krit}']), int(b[f'kedalaman_{krit}'])
    if L == 'kurang':
        L = f'>{BATAS[blok]}'
    return L, None, None


def batas_bawah(m, blok, na, krit='terbaik'):
    L, N, _ = bintang(m, blok, na, krit)
    if N is not None:
        return N, True
    cap = 12 if (blok, m) == ('hea', 'BeH2') else BATAS[blok]
    return haa.biaya(blok, ns[m], na, cap + 1)[0], False


for m in MOL:
    for nama, (blok, na) in DESAIN.items():
        L, N, D = bintang(m, blok, na)
        nilai[f'L-{m}-{nama}'] = L if isinstance(L, int) else L.replace('>', r'\ensuremath{>}')
        nilai[f'N-{m}-{nama}'] = N if N is not None else '--'
        nilai[f'D-{m}-{nama}'] = D if D is not None else '--'
    nilai[f'hemat-{m}'] = round(100 * (1 - nilai[f'N-{m}-so4-brick'] / nilai[f'N-{m}-can-brick']))
hemat = [nilai[f'hemat-{m}'] for m in MOL]
nilai['hemat-min'], nilai['hemat-maks'] = min(hemat), max(hemat)
hemat = [round(100 * (1 - nilai[f'N-{m}-so4-star'] / nilai[f'N-{m}-can-star'])) for m in MOL]
nilai['hemat-star-min'], nilai['hemat-star-maks'] = min(hemat), max(hemat)
luas = [(m, top) for m in MOL for top in ('star', 'brick')]
nilai['n-so4-can-kurang'] = sum(nilai[f'L-{m}-so4-{top}'] < nilai[f'L-{m}-can-{top}'] for m, top in luas)
nilai['n-so4-can-sama'] = sum(nilai[f'L-{m}-so4-{top}'] == nilai[f'L-{m}-can-{top}'] for m, top in luas)
nilai['n-so4-can-murah'] = sum(nilai[f'N-{m}-so4-{top}'] < nilai[f'N-{m}-can-{top}'] for m, top in luas)
nilai['n-hemat-sepertiga'] = sum(3 * nilai[f'N-{m}-so4-{top}'] == 2 * nilai[f'N-{m}-can-{top}'] for m, top in luas)
banding = [lapis(m, top, 'so4') - lapis(m, top, 'can') for m, top in sel]
banding += [nilai[f'L-{m}-so4-{top}'] - nilai[f'L-{m}-can-{top}'] for m, top in luas]
nilai['n-banding-lapis'] = len(banding)
nilai['n-so4-lapis-kurang'] = sum(x < 0 for x in banding)
nilai['n-so4-lapis-sama'] = sum(x == 0 for x in banding)
nilai['n-so4-lapis-lebih'] = sum(x > 0 for x in banding)
for hea in ('heau3', 'hea'):
    akhir = '' if hea == 'hea' else '-u3'
    nilai[f'n-can-lapis-kurang{akhir}'] = sum(min(nilai[f'L-{m}-can-star'], nilai[f'L-{m}-can-brick']) < nilai[f'L-{m}-{hea}']
                                             for m in MOL)
    nilai[f'n-can-lebih{akhir}'] = sum(min(nilai[f'N-{m}-can-star'], nilai[f'N-{m}-can-brick']) > nilai[f'N-{m}-{hea}']
                                      for m in MOL)
    nilai[f'n-so4-brick-dangkal{akhir}'] = sum(nilai[f'D-{m}-so4-brick'] < nilai[f'D-{m}-{hea}'] for m in MOL)
nilai['n-can-lebih-hea'] = nilai['n-can-lebih']
nilai['n-so4-brick-le-hea'] = sum(nilai[f'N-{m}-so4-brick'] <= nilai[f'N-{m}-hea'] for m in MOL)


def menang(m, krit):
    so4, ok = batas_bawah(m, 'so4', 0, krit)
    hea, ok_hea = batas_bawah(m, 'heau3', 0, krit)
    if ok and ok_hea:
        return so4 <= hea
    if ok and not ok_hea:
        return so4 <= hea
    return False


W = {krit: [m for m in MOL if menang(m, krit)] for krit in ('terbaik', 'tiga')}
nilai['W'] = min(len(W['terbaik']), len(W['tiga']))
nilai['W-kalah'] = ', '.join(m for m in MOL if m not in W['terbaik'] or m not in W['tiga'])
for top in ('canu3', 'so4u3'):
    for m in MOL:
        for na, nama in ((1, 'star'), (0, 'brick')):
            L, N, _ = bintang(m, top, na, rentang='asli')
            nilai[f'La-{m}-{top}-{nama}'], nilai[f'Na-{m}-{top}-{nama}'] = L, N
kontrol = [(m, top) for m in MOL for top in ('star', 'brick')]
nilai['n-u3-can-murah'] = sum(nilai[f'Na-{m}-canu3-{top}'] < nilai[f'Na-{m}-so4u3-{top}'] for m, top in kontrol)
nilai['n-u3-so4-lapis-kurang'] = sum(nilai[f'La-{m}-so4u3-{top}'] < nilai[f'La-{m}-canu3-{top}'] for m, top in kontrol)
nilai['n-u3-lapis-sama'] = sum(nilai[f'La-{m}-so4u3-{top}'] == nilai[f'La-{m}-canu3-{top}'] for m, top in kontrol)
sel_lstar = [(m, d) for m in MOL for d in DESAIN]
nilai['n-tiga-sel'] = len(sel_lstar)
nilai['n-tiga-berubah'] = sum(bintang(m, *DESAIN[d])[0] != bintang(m, *DESAIN[d], krit='tiga')[0] for m, d in sel_lstar)
for m, d in sel_lstar:
    L3, N3, D3 = bintang(m, *DESAIN[d], krit='tiga')
    nilai[f'D3-{m}-{d}'] = D3 if D3 is not None else '--'
nilai['n-tiga-banding-berubah'] = 0 if W['terbaik'] == W['tiga'] else 1
eig = {}
# Sampled grids only; the exact endpoints and rule-7 checks of ujung.py (ujung_*, cek_*) are counted separately.
for f in glob.glob(os.path.join(haa.HASIL, 'pec', '*_taper_*.csv')) + \
        [g for g in glob.glob(os.path.join(haa.HASIL, 'uji', '[!a]*.csv'))
         if not os.path.basename(g).startswith(('ujung_', 'cek_'))]:
    for b in haa.baca(f):
        k = (b['molekul'], round(float(b['R']), 3))
        eig[k] = max(eig.get(k, 0.0), abs(float(b['E_dasar']) - float(b['E_fci'])))
nilai['n-geometri'] = len(eig)
nilai['n-geometri-buang'] = sum(x > 1e-8 for x in eig.values())
nilai['geometri-selisih-maks'] = ilmiah(max(eig.values()))
for (m, blok, na, rentang), b in uji.items():
    if m == 'BeH2' and rentang == UTAMA:
        nilai[f'tiga-{m}-{blok}-na{na}'] = (b['L_tiga'] if b['L_tiga'] != 'kurang' else f'>{BATAS[blok]}').replace('>', r'\ensuremath{>}')
jauh = {}
for f in glob.glob(os.path.join(haa.HASIL, 'pec', '*BeH2*_taper_*.csv')) + glob.glob(os.path.join(haa.HASIL, 'beh2', '*.csv')) + \
        [g for g in glob.glob(os.path.join(haa.HASIL, 'uji', '*_BeH2_*.csv')) if not os.path.basename(g).startswith('cek_')]:
    for b in haa.baca(f):
        if b.get('peta', 'taper') == 'taper' and b.get('opt', 'lbfgsb') == 'lbfgsb' and b['molekul'] == 'BeH2':
            jauh.setdefault((b['blok'], int(b['na']), round(float(b['R']), 3), int(b['L'])), []).append(float(b['galat_mHa']))
RJ = sorted(set(haa.JARAK['BeH2']) | {0.7, 2.5, 3.0, 3.3})
for blok, na in (('hea', 0), ('heau3', 0)):
    for syarat, akhir in ((1, ''), (3, '-tiga')):
        Ls = sorted({k[3] for k in jauh if k[:2] == (blok, na)})
        lolos = [L for L in Ls if all((blok, na, R, L) in jauh and sum(g < 1e3 * haa.CHEM for g in jauh[(blok, na, R, L)]) >= syarat
                                      for R in RJ)]
        nilai[f'beh2-jarak-{blok}-na{na}{akhir}'] = lolos[0] if lolos else rf'\ensuremath{{>}}{Ls[-1]}'
        if syarat == 1 and lolos:
            nilai[f'beh2-jarak-{blok}-na{na}-N'] = haa.biaya(blok, 7, na, lolos[0])[0]
            nilai[f'beh2-jarak-{blok}-na{na}-nR'] = len(RJ)
            gagal = [str(L) for L in Ls if L > lolos[0] and L not in lolos]
            nilai[f'beh2-jarak-{blok}-na{na}-gagal'] = ', '.join(gagal)
            # Prose lists: "10, 11, or 12" and "10, 11, and 12".
            for kata in ('atau', 'dan'):
                hubung = {'atau': 'or', 'dan': 'and'}[kata]
                nilai[f'beh2-jarak-{blok}-na{na}-gagal-{kata}'] = (
                    gagal[0] if len(gagal) == 1 else f'{gagal[0]} {hubung} {gagal[1]}' if len(gagal) == 2
                    else ', '.join(gagal[:-1]) + f', {hubung} {gagal[-1]}') if gagal else ''
awal_ = haa.baca('tabel/uji_awal.csv')
for grup, desain in (('can', ('can-star', 'can-brick')), ('so4', ('so4-star', 'so4-brick')), ('hea', ('hea',))):
    for jenis in ('m0', 'nol', 'kecil'):
        for opt in ('lbfgsb', 'cobyla', 'spsa', 'spsa-shot'):
            bs = [b for b in awal_ if b['desain'] in desain and b['awal'] == jenis and b['opt'] == opt]
            if bs:
                nilai[f'uji-awal-{grup}-{jenis}-{opt}'] = f"{sum(int(b['sukses']) for b in bs)}/{sum(int(b['n']) for b in bs)}"
                g = [float(b['terbaik_mHa']) for b in bs]
                nilai[f'uji-awal-{grup}-{jenis}-{opt}-min'] = f'{min(g):.1f}'
                nilai[f'uji-awal-{grup}-{jenis}-{opt}-maks'] = f'{max(g):.1f}'

KATA_ = ('zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten')

# Numbers that were typed into the text (open record m9 of the referee report).
h4 = [b for b in haa.baca('sistem.csv') if b['molekul'] == 'H4']
assert len({b['suku_jw'] for b in h4}) == 1 and len({b['suku_taper'] for b in h4}) == 1
nilai['h4-suku-jw'], nilai['h4-suku-taper'] = int(h4[0]['suku_jw']), int(h4[0]['suku_taper'])
nilai['h4-sektor-maks'] = ilmiah(max(abs(float(b['E0']) - float(b['E_fci'])) for b in h4))
p_, d_ = 1e-3, 4
nilai['pec-gamma-cnot'] = f'{(1 + p_ * (d_ * d_ - 2) / (d_ * d_)) / (1 - p_):.5f}'
peta = [b for b in haa.baca('tabel/beh2.csv') if b['faktor'] == 'peta']
nilai['beh2-peta-galat-min'] = f"{min(max(float(b['terbaik_mHa']) for b in peta if b['L'] == L) for L in {b['L'] for b in peta}):.0f}"


# SI table tab:dua: both success criteria, with the rule-4 lower bound (CNOTs one layer above the last tested
# count) and the worst-geometry error at that count for designs that fail. 't2' = sampled extended grids ('luas');
# 't2d' = registered range of rule 2 ('daftar', exact 2.5 R0 endpoints from ujung.py), written once complete.
def dua_kriteria(rentang, awalan):
    lengkap = True
    for m in MOL:
        for d, (blok, na) in DESAIN.items():
            b = uji.get((m, blok, na, rentang))
            for krit in ('terbaik', 'tiga'):
                k = f'{awalan}-{m}-{d}-{krit}'
                L = b[f'L_{krit}'] if b else 'kurang'
                if L.isdigit():
                    nilai[k + '-L'], nilai[k + '-N'] = int(L), int(b[f'cnot_{krit}'])
                    nilai[k + '-D'], nilai[k + '-galat'] = int(b[f'kedalaman_{krit}']), '--'
                elif L.startswith('>'):
                    nilai[k + '-L'] = rf'\ensuremath{{>}}{L[1:]}'
                    nilai[k + '-N'] = rf'\ensuremath{{\geq}}{haa.biaya(blok, ns[m], na, int(L[1:]) + 1)[0]}'
                    g = b.get(f'galat_batas_{krit}', '')
                    nilai[k + '-D'], nilai[k + '-galat'] = '--', f'{float(g):.1f}' if g else '--'
                else:
                    lengkap = False
                    for s in ('-L', '-N', '-D', '-galat'):
                        nilai[k + s] = 'n/a'
    return lengkap


dua_kriteria(UTAMA, 't2')

# Rule-7 checks of ujung.py / ujung_cluster.py on every geometry of the registered range.
cek = [b for f in glob.glob(os.path.join(haa.HASIL, 'uji', 'cek_*.csv')) for b in haa.baca(f)]
if cek:
    buang = sorted((b['molekul'], float(b['R'])) for b in cek if b['buang'] == '1')
    nilai['n-geometri-daftar'] = len(cek)
    nilai['n-buang-daftar'] = len(buang)
    nilai['n-buang-daftar-kata'] = KATA_[len(buang)] if len(buang) <= 10 else len(buang)
    nilai['buang-geometri'] = ', '.join(rf"{m.replace('H2', 'H$_2$').replace('H4', 'H$_4$')} at {R:g}~\AA"
                                        for m, R in buang) or 'none'
    nilai['n-rhf-gagal'] = sum(b['rhf_konvergen'] != '1' or b['driver_konvergen'] == '0' for b in cek)
    nilai['cek-selisih-maks'] = ilmiah(max(float(b['selisih']) for b in cek))
    nilai['buang-selisih'] = ilmiah(max(float(b['selisih']) for b in cek if b['buang'] == '1')) if buang else '--'
# Do the sampled grids (rounded endpoints) and the registered range give the same resources in every cell?
kol_ = ('L_terbaik', 'cnot_terbaik', 'kedalaman_terbaik', 'L_tiga', 'cnot_tiga', 'kedalaman_tiga')
if UTAMA == 'daftar':
    nilai['luas-sama-daftar'] = int(all(uji[(m, b, na, 'luas')][c] == uji[(m, b, na, 'daftar')][c]
                                        for m in MOL for b, na in DESAIN.values() for c in kol_))
if any(r == 'daftar' for (_, _, _, r) in uji) and dua_kriteria('daftar', 't2d'):
    # Verdicts on the registered range (rules 3 and 4), for the text once the exact endpoints exist.
    def N_atau_batas(m, blok, na, krit):
        b = uji[(m, blok, na, 'daftar')]
        L = b[f'L_{krit}']
        if L.isdigit():
            return int(b[f'cnot_{krit}']), True
        return haa.biaya(blok, ns[m], na, int(L[1:]) + 1)[0], False

    def menang_daftar(m, krit):
        (so4, ok), (hea, _) = N_atau_batas(m, 'so4', 0, krit), N_atau_batas(m, 'heau3', 0, krit)
        return ok and so4 <= hea

    Wd = {krit: [m for m in MOL if menang_daftar(m, krit)] for krit in ('terbaik', 'tiga')}
    nilai['W-daftar'] = min(len(Wd['terbaik']), len(Wd['tiga']))
    nilai['W-daftar-kata'] = KATA_[nilai['W-daftar']]
    nilai['W-daftar-kalah'] = ', '.join(m for m in MOL if m not in Wd['terbaik'] or m not in Wd['tiga'])
    nilai['daftar-sama'] = int(nilai['W-daftar'] == nilai['W'])
    for top, akhir in (('brick', ''), ('star', '-star')):
        pasang = [(nilai[f't2d-{m}-so4-{top}-terbaik-N'], nilai[f't2d-{m}-can-{top}-terbaik-N']) for m in MOL]
        if all(isinstance(s, int) and isinstance(c, int) for s, c in pasang):
            h = [round(100 * (1 - s / c)) for s, c in pasang]
            nilai[f'hemat{akhir}-daftar-min'], nilai[f'hemat{akhir}-daftar-maks'] = min(h), max(h)
    print('registered range (daftar): W =', nilai['W-daftar'], '(lost:', nilai['W-daftar-kalah'] or 'none', ')',
          '| savings', {k: v for k, v in nilai.items() if k.startswith('hemat') and 'daftar' in k})
else:
    print('registered range (daftar): incomplete; run ujung.py and uji.py tabel first')

# Spelled-out forms of small counts for prose ("three of the four molecules"), key + '-kata'.
KATA = ('zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten')
for k, v in list(nilai.items()):
    if (k == 'W' or k.startswith('n-')) and str(v).isdigit() and int(v) <= 10:
        nilai[f'{k}-kata'] = KATA[int(v)]

with open(keluar, 'w') as f:
    for k, v in nilai.items():
        f.write(f'\\expandafter\\def\\csname val@{k}\\endcsname{{{v}}}\n')
print(len(nilai), keluar)
