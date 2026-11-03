"""Rebuild the results/ folder that the analysis scripts read from the collected files of kumpulkan_hasil.py.
The collected CSVs may also be gzipped (name.csv.gz).

usage: python src/pulihkan.py [data dir, default data/] [results dir, default results/]
"""
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'data'
HASIL = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / 'results'


def buka(nama):
    if (DATA / nama).exists():
        return open(DATA / nama, newline='', encoding='utf-8')
    return gzip.open(DATA / (nama + '.gz'), 'rt', newline='', encoding='utf-8')


with buka('manifest.csv') as f:
    manifest = list(csv.DictReader(f))
baris = defaultdict(list)
for k in sorted({m['collected_file'] for m in manifest}):
    with buka(k) as f:
        for b in csv.DictReader(f):
            baris[b['source']].append(b)
for m in manifest:
    kepala = m['header'].split(',')
    p = HASIL / m['source']
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(kepala)
        w.writerows([b[c] for c in kepala] for b in baris[m['source']])
print(f'{len(manifest)} files written to {HASIL}')
