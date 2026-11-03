import glob
import os
import sys
from collections import defaultdict

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import paper_data as haa

KELUAR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(haa.HASIL, '..', 'figures', 'jctc')
os.makedirs(KELUAR, exist_ok=True)
MOL = ('H4', 'LiH', 'HF', 'BeH2')
NAMA_MOL = {'H4': r'H$_4$', 'LiH': 'LiH', 'HF': 'HF', 'BeH2': r'BeH$_2$'}
TAMBAH = haa.TAMBAH
CHEM = haa.CHEM
MERAH, ABU, BIRU, BIRU_MUDA = '#B23A28', '#606060', '#005A8D', '#007766'
GAYA = {('zeng', 1): ('HAA-CAN, star', ABU, '-', 'o'), ('zeng', 0): ('HAA-CAN, brick', ABU, '--', 's'),
        ('so4', 1): ('HAA-SO(4), star', MERAH, '-', 'o'), ('so4', 0): ('HAA-SO(4), brick', MERAH, '--', 's'),
        ('heau3', 0): (r'HEA, $U_3$', BIRU, '-', '^'), ('hea', 0): (r'HEA, $R_YR_Z$', BIRU_MUDA, ':', 'v')}
plt.rcParams.update({'font.family': 'Arial', 'mathtext.fontset': 'stixsans', 'font.size': 9, 'axes.labelsize': 9,
                     'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'axes.spines.top': True,
                     'axes.spines.right': True, 'savefig.dpi': 600, 'pdf.fonttype': 42,
                     'svg.fonttype': 'none', 'axes.linewidth': 0.7})


def simpan(fig, nama):
    for ext in ('pdf', 'svg'):
        fig.savefig(os.path.join(KELUAR, f'{nama}.{ext}'))
    fig.savefig(os.path.join(KELUAR, f'{nama}.png'), dpi=180)
    plt.close(fig)


LANTAI = 1e-6


def sumbu_galat(ax, label=True):
    ax.set_yscale('log')
    ax.set_ylim(LANTAI * 0.7, 2)
    ax.axhline(CHEM, color='0.35', lw=0.8, ls=(0, (4, 3)), zorder=0)
    if label:
        ax.set_ylabel('energy error (hartree)')


def kumpul(cocok=lambda b: True):
    data = defaultdict(list)
    sumber = glob.glob(os.path.join(haa.HASIL, 'pec', '*_taper_*.csv')) + glob.glob(os.path.join(haa.HASIL, 'beh2', '*.csv'))
    sumber += [f for f in glob.glob(os.path.join(haa.HASIL, 'uji', '*.csv')) if not os.path.basename(f).startswith(('awal_', 'cek_'))]
    for f in sumber:
        for b in haa.baca(f):
            if b.get('opt', 'lbfgsb') != 'lbfgsb' or not cocok(b):
                continue
            data[(b['molekul'], b['blok'], int(b['na']), b.get('peta', 'taper'), round(float(b['R']), 3), int(b['L']))].append(
                float(b['galat_mHa']) / 1e3)
    return data


def kurva(data, mol, blok, na, Rs, peta='taper'):
    Ls = sorted({k[5] for k in data if k[:4] == (mol, blok, na, peta)})
    titik = [(L, max(min(data[(mol, blok, na, peta, R, L)]) for R in Rs)) for L in Ls
             if all((mol, blok, na, peta, R, L) in data for R in Rs)]
    return [t[0] for t in titik], [max(t[1], LANTAI) for t in titik]


def biaya():
    data = kumpul()
    ns = {m: int(b['qubit']) for b in haa.baca('tabel/tab2_biaya.csv') for m in [b['molekul']] if b['topologi'] == 'brick'}
    fig, axs = plt.subplots(2, 2, figsize=(6.6, 5.3), sharey=True)
    dibuang = haa.buang()
    for ax, mol in zip(axs.flat, MOL):
        # Registered range of rule 2 (exact 2.5 R0 endpoints, rule-7 geometries removed), as in tab:lstar.
        Rs = [round(R, 3) for R in haa.daftar(mol, dibuang)]
        for (blok, na), (nama, warna, garis, tanda) in GAYA.items():
            Ls, g = kurva(data, mol, blok, na, Rs)
            N = [haa.biaya(blok, ns[mol], na, L)[0] for L in Ls]
            ax.plot(N, g, color=warna, ls=garis, marker=tanda, ms=3.2, lw=1.1, label=nama)
        sumbu_galat(ax, label=ax in axs[:, 0])
        ax.set_title(f'({"abcd"[MOL.index(mol)]}) {NAMA_MOL[mol]}', fontsize=9, loc='left')
        ax.set_xlim(0, None)
        if ax in axs[1]:
            ax.set_xlabel(r'CNOTs $N_\mathrm{2q}$')
    fig.legend(*axs[0, 0].get_legend_handles_labels(), frameon=False, loc='upper center', ncol=3)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    simpan(fig, 'hasil_biaya')


def derau():
    baris = haa.baca('derau.csv')
    PERAN = {'hea': (r'HEA, $R_YR_Z$', BIRU_MUDA, 'v'), 'can_terbaik': ('HAA-CAN', ABU, 'o'),
             'so4_star': ('HAA-SO(4), star', MERAH, 'o'), 'so4_brick': ('HAA-SO(4), brick', MERAH, 's')}
    titik = [(m, g) for m in MOL for g in ('setimbang', 'regang')]
    nilai = {(b['molekul'], b['geometri'], b['peran']): b for b in baris}
    fig, ax = plt.subplots(figsize=(6.6, 2.7))
    for j, (peran, (nama, warna, tanda)) in enumerate(PERAN.items()):
        isi = 'none' if peran == 'so4_brick' else warna
        x = [i + (j - 1.5) * 0.12 for i in range(len(titik))]
        y = [float(nilai[t + (peran,)]['p_derau']) for t in titik]
        ax.scatter(x, y, s=26, marker=tanda, facecolors=isi, edgecolors=warna, lw=1.1, label=nama, zorder=3)
    for k in (1, 3):
        ax.axvspan(2 * k - 0.5, 2 * k + 1.5, color='0.95', zorder=0, lw=0)
    ax.set_yscale('log')
    ax.set_xticks(range(len(titik)))
    def jarak(R):
        s = f'{float(R):g}'
        return s if '.' in s else s + '.0'
    ax.set_xticklabels([f"{NAMA_MOL[m]}\n{jarak(nilai[(m, g, 'hea')]['R'])} Å" for m, g in titik])
    ax.set_xlim(-0.6, len(titik) - 0.4)
    ax.set_ylabel(r'noise allowance $p_\mathrm{noise}$')
    ax.legend(frameon=False, loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.16))
    fig.tight_layout()
    simpan(fig, 'hasil_derau')


def awal():
    LSTAR = {'hea': 4, 'can-star': 2, 'can-brick': 3, 'so4-star': 2, 'so4-brick': 2}
    run = defaultdict(list)
    for f in glob.glob(os.path.join(haa.HASIL, 'uji', 'awal_*.csv')):
        for b in haa.baca(f):
            jenis = {'m0': f"m0-{b['opt']}", 'kecil': 'kecil'}.get(b['awal'])
            if jenis:
                run[(b['desain'], jenis)].append(float(b['galat_mHa']) / 1e3)
    for b in haa.baca('h4_optim.csv'):
        if b['awal'] in ('m0', 'acak') and float(b['R']) in (0.9, 1.5) and int(b['L']) == LSTAR[b['desain']]:
            run[(b['desain'], 'm0-lbfgsb' if b['awal'] == 'm0' else 'acak')].append(float(b['galat']))
    JENIS = [('acak', 'random\nL-BFGS-B'), ('kecil', r'near zero' '\nL-BFGS-B'), ('m0-lbfgsb', r'in $\mathcal{M}_0$' '\nL-BFGS-B'),
             ('m0-cobyla', r'in $\mathcal{M}_0$' '\nCOBYLA'), ('m0-spsa', r'in $\mathcal{M}_0$' '\nSPSA')]
    DESAIN = [('can-star', 'HAA-CAN, star', ABU, 'o'), ('can-brick', 'HAA-CAN, brick', ABU, 's'),
              ('so4-star', 'HAA-SO(4), star', MERAH, 'o'), ('so4-brick', 'HAA-SO(4), brick', MERAH, 's'),
              ('hea', r'HEA, $R_YR_Z$', BIRU_MUDA, 'v')]
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    lebar = 0.15
    rng = np.random.default_rng(0)
    for j, (d, nama, warna, tanda) in enumerate(DESAIN):
        for i, (jenis, _) in enumerate(JENIS):
            g = np.maximum(np.array(run[(d, jenis)]), LANTAI)
            x = i + (j - 2) * lebar + rng.uniform(-0.04, 0.04, len(g))
            isi = 'none' if d.endswith('brick') else warna
            ax.scatter(x, g, s=9, marker=tanda, facecolors=isi, edgecolors=warna, lw=0.7, label=nama if i == 0 else None,
                       zorder=3)
    sumbu_galat(ax)
    ax.set_xticks(range(len(JENIS)))
    ax.set_xticklabels([t for _, t in JENIS])
    fig.legend(*ax.get_legend_handles_labels(), frameon=False, loc='upper center', ncol=3, columnspacing=1.0,
               handletextpad=0.3)
    fig.tight_layout(rect=(0, 0, 1, 0.84))
    simpan(fig, 'hasil_awal')


def beh2():
    data = kumpul(lambda b: b['molekul'] == 'BeH2')
    pendek = haa.JARAK['BeH2']
    lebar = sorted(set(pendek) | {0.7, 2.5, 3.0, 3.3})
    SERI = [('hea', 'taper', pendek, r'$R_YR_Z$, 1.0–2.0 Å', BIRU_MUDA, '-', 'o'),
            ('hea', 'taper', lebar, r'$R_YR_Z$, 0.7–3.3 Å', BIRU_MUDA, '--', 's'),
            ('heau3', 'taper', pendek, r'$U_3$, 1.0–2.0 Å', BIRU, '-', 'o'),
            ('heau3', 'taper', lebar, r'$U_3$, 0.7–3.3 Å', BIRU, '--', 's'),
            ('hea', 'jw', [1.0, 1.326, 2.0], r'$R_YR_Z$, no tapering', '0.3', ':', 'D')]
    fig, ax = plt.subplots(figsize=(3.33, 3.5))
    for blok, peta, Rs, nama, warna, garis, tanda in SERI:
        Ls, g = kurva(data, 'BeH2', blok, 0, Rs, peta)
        ax.plot(Ls, g, color=warna, ls=garis, marker=tanda, ms=3.2, lw=1.1, label=nama)
    ax.axvline(25, color=ABU, lw=0.8)
    ax.text(24.4, 0.5, 'Zeng et al.', rotation=90, ha='right', va='top', fontsize=7, color=ABU)
    sumbu_galat(ax)
    ax.set_ylim(5e-5, 5)
    ax.set_xlabel('layers $L$')
    ax.legend(frameon=False, loc='upper center', fontsize=7, ncol=2, bbox_to_anchor=(0.5, -0.2), columnspacing=0.8)
    fig.tight_layout()
    simpan(fig, 'hasil_beh2')


if __name__ == '__main__':
    for f in (biaya, derau, awal, beh2):
        f()
        print(f.__name__, 'ok')
