"""Copy a vector figure with every stroke thinner than 0.5 pt raised to 0.5 pt (JCTC: "Lines should be no thinner
than 0.5 point"). Only the line-width operators of the page content change; the original file is not touched.
Valid for figures drawn at final size without scaling (CTM scale 1), as the PhotoCraft/VectorCraft exports are.

usage: python src/garis_minimum.py <in.pdf> <out.pdf> [minimum width in pt, default 0.5]
"""
import re
import sys

import pymupdf

masuk, keluar = sys.argv[1], sys.argv[2]
minimum = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
doc = pymupdf.open(masuk)
halaman = doc[0]
for cm in re.findall(rb'((?:-?\d*\.?\d+ ){6})cm', b''.join(doc.xref_stream(x) for x in halaman.get_contents())):
    a, b, c, d = (float(v) for v in cm.split()[:4])
    if abs(abs(a * d - b * c) - 1) > 1e-6:
        sys.exit(f'scaled coordinate system ({cm.decode()}); widths are not in pt, nothing written')
n = 0
for xref in halaman.get_contents():
    teks = doc.xref_stream(xref)

    def naikkan(m):
        global n
        if float(m.group(1)) < minimum:
            n += 1
            return f'{minimum:g} w'.encode()
        return m.group(0)
    doc.update_stream(xref, re.sub(rb'(?<![\d.])(\d*\.?\d+) w\b', naikkan, teks))
doc.save(keluar, garbage=3, deflate=True)
print(f'{n} line widths raised to {minimum:g} pt: {keluar}')
