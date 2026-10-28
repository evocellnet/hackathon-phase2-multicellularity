#!/usr/bin/env python3
"""Convert one sheet of an .ods spreadsheet to TSV using only the stdlib (no pandas/odfpy).
Usage: ods_to_tsv.py file.ods sheet_index out.tsv"""
import sys, zipfile, xml.etree.ElementTree as ET
NS = {'table': 'urn:oasis:names:tc:opendocument:xmlns:table:1.0',
      'text': 'urn:oasis:names:tc:opendocument:xmlns:text:1.0',
      'office': 'urn:oasis:names:tc:opendocument:xmlns:office:1.0'}
T = '{%s}' % NS['table']
def cell_text(c):
    ps = c.findall('.//text:p', NS)
    return ' '.join(''.join(p.itertext()) for p in ps).replace('\t', ' ').replace('\n', ' ')
def main(path, idx, out):
    root = ET.fromstring(zipfile.ZipFile(path).read('content.xml'))
    sheets = root.findall('.//table:table', NS)
    print('sheets:', [s.get(T + 'name') for s in sheets], file=sys.stderr)
    sh = sheets[int(idx)]
    rows = []
    for r in sh.findall('.//table:table-row', NS):
        rrep = int(r.get(T + 'number-rows-repeated', 1))
        vals = []
        for c in r:
            if c.tag not in (T + 'table-cell', T + 'covered-table-cell'):
                continue
            rep = int(c.get(T + 'number-columns-repeated', 1))
            v = cell_text(c) or c.get('{%s}value' % NS['office'], '')
            vals.extend([v] * min(rep, 500))
        while vals and vals[-1] == '':
            vals.pop()
        if not vals:
            continue
        rows.extend([vals] * min(rrep, 1))
    with open(out, 'w') as fh:
        for v in rows:
            fh.write('\t'.join(v) + '\n')
    print(len(rows), 'rows', file=sys.stderr)
if __name__ == '__main__':
    main(*sys.argv[1:4])
