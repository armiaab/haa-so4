"""Lightweight access to manuscript data without importing simulation packages.

Geometry constants are read from the simulation source using literal AST values.
The cost formulas describe generic compiled blocks, not hardware routing.
"""
import ast
import csv
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HASIL = str(ROOT / 'results')


def constants(filename, names):
    found = {}
    for node in ast.parse((ROOT / 'src' / filename).read_text(encoding='utf-8')).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in names:
                    found[target.id] = ast.literal_eval(node.value)
    if set(found) != set(names):
        raise ValueError(f'Missing constants in {filename}: {set(names) - set(found)}')
    return found


globals().update(constants('haa.py', {'JARAK', 'SETIMBANG', 'CHEM'}))
TAMBAH = constants('uji.py', {'TAMBAH'})['TAMBAH']
UJUNG = constants('uji.py', {'UJUNG'})['UJUNG']


def daftar(mol, buang=()):
    """Registered range of rule 2 (as uji.daftar), without the geometries removed under rule 7."""
    ujung = round(2.5 * SETIMBANG[mol], 4)
    tambah = TAMBAH[mol] + UJUNG.get(mol, [])
    Rs = set(JARAK[mol]) | {R for R in tambah if min(JARAK[mol]) <= R <= ujung + 1e-9}
    return sorted(R for R in Rs if (mol, R) not in set(buang))


def buang():
    """(molecule, R) removed under rule 7, from the checks of ujung.py / ujung_cluster.py."""
    import glob
    return sorted((b['molekul'], float(b['R'])) for f in glob.glob(os.path.join(HASIL, 'uji', 'cek_*.csv'))
                  for b in baca(f) if b['buang'] == '1')


def baca(nama):
    with open(os.path.join(HASIL, nama), newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def biaya(nama, ns, na, L):
    if L < 1 or ns < 3 or na not in (0, 1):
        raise ValueError('Manuscript cost formulas require L >= 1, ns >= 3, na in {0,1}')
    if nama in ('hea', 'heau3', 'ra'):
        if na:
            raise ValueError('HEA baselines do not use ancillas')
        return (ns - 1) * L, ns + 2 * L - 3
    c = {'can': 3, 'zeng': 3, 'canu3': 3, 'real': 2, 'so4': 2, 'so4u3': 2}[nama]
    return c * (ns if na else ns - 1) * L, c * (ns if na else 2) * L
