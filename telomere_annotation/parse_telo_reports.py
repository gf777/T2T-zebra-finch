#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_telo_reports.py

Iterate a root folder and parse every "*.telo.report" file into a single CSV table.

Expected output columns (one row per accession):
accession,total_paths,total_gaps,total_telomeres,mean_length,median_length,min_length,max_length,two_telomeres,one_telomere,zero_telomeres,t2t,gapped_t2t,missassembled,gapped_missassembled,incomplete,gapped_incomplete,no_telomeres,gapped_no_telomeres,discordant,gapped_discordant

Notes
-----
• The "+++ Path Summary Report +++" section is variable-length across assemblies.
• If the "Chromosome Telomere/Gap Completeness" section is missing, category counts
  are computed from the per-path 'type' column of the table.
• Length stats (mean/median/min/max) are parsed if present; otherwise left blank.
• Accessions are derived from the first directory component that looks like GCA_*/GCF_*,
  or else from the filename prefix up to the first underscore.

"""
import argparse
import sys
from pathlib import Path
import re
import pandas as pd
from typing import Dict, Any, List

# --- Palette, updated names (for downstream plotting; not written to CSV) ---
PALETTE = {
    't2t'                : '#1A9641',
    'gapped_t2t'         : '#9CCF60',
    'incomplete'         : '#FFC754',
    'gapped_incomplete'  : '#FFE885',
    'missassembled'      : '#D6594C',
    'gapped_missassembled':'#F58B6D',
    'discordant'         : '#8278F4',
    'gapped_discordant'  : '#B395EB',
    'no_telomeres'       : '#C8C8C8',
    'gapped_no_telomeres': '#F0F0F0',
}

# Category keys used in output
CATEGORY_KEYS = [
    't2t', 'gapped_t2t', 'missassembled', 'gapped_missassembled',
    'incomplete', 'gapped_incomplete', 'no_telomeres', 'gapped_no_telomeres',
    'discordant', 'gapped_discordant'
]

RE_INT = re.compile(r'(\d+)')
RE_FLOAT = re.compile(r'([0-9]+(?:\.[0-9]+)?)')

def _to_key(s: str) -> str:
    """Convert category label to canonical key (e.g., 'Gapped T2T' -> 'gapped_t2t')."""
    return s.strip().lower().replace(' ', '_')

def parse_report_text(text: str) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        'total_paths': None, 'total_gaps': None, 'total_telomeres': None,
        'mean_length': None, 'median_length': None, 'min_length': None, 'max_length': None,
        'two_telomeres': None, 'one_telomere': None, 'zero_telomeres': None,
        't2t': 0, 'gapped_t2t': 0, 'missassembled': 0, 'gapped_missassembled': 0,
        'incomplete': 0, 'gapped_incomplete': 0, 'no_telomeres': 0, 'gapped_no_telomeres': 0,
        'discordant': 0, 'gapped_discordant': 0,
    }

    lines = [ln.rstrip('\n') for ln in text.splitlines()]

    # Simple "key: value" lines anywhere in file
    for ln in lines:
        ll = ln.lower()
        if ll.startswith('total paths:'):
            m = RE_INT.search(ln);          out['total_paths'] = int(m.group(1)) if m else out['total_paths']
        elif ll.startswith('total gaps:'):
            m = RE_INT.search(ln);          out['total_gaps'] = int(m.group(1)) if m else out['total_gaps']
        elif ll.startswith('total telomeres:'):
            m = RE_INT.search(ln);          out['total_telomeres'] = int(m.group(1)) if m else out['total_telomeres']
        elif ll.startswith('mean length:'):
            m = RE_FLOAT.search(ln);        out['mean_length'] = float(m.group(1)) if m else out['mean_length']
        elif ll.startswith('median length:'):
            m = RE_FLOAT.search(ln);        out['median_length'] = float(m.group(1)) if m else out['median_length']
        elif ll.startswith('min length:'):
            m = RE_FLOAT.search(ln);        out['min_length'] = float(m.group(1)) if m else out['min_length']
        elif ll.startswith('max length:'):
            m = RE_FLOAT.search(ln);        out['max_length'] = float(m.group(1)) if m else out['max_length']
        elif ll.startswith('two telomeres:'):
            m = RE_INT.search(ln);          out['two_telomeres'] = int(m.group(1)) if m else out['two_telomeres']
        elif ll.startswith('one telomere:'):
            m = RE_INT.search(ln);          out['one_telomere'] = int(m.group(1)) if m else out['one_telomere']
        elif ll.startswith('zero telomeres:'):
            m = RE_INT.search(ln);          out['zero_telomeres'] = int(m.group(1)) if m else out['zero_telomeres']
        else:
            # Category lines like "Gapped T2T: 17"
            if ':' in ln:
                lhs, rhs = ln.split(':', 1)
                cat = _to_key(lhs)
                if cat in CATEGORY_KEYS:
                    m = RE_INT.search(rhs)
                    if m: out[cat] = int(m.group(1))

    # Fallback via the Path Summary table if needed
    need_categories = any(out[k] is None or out[k] == 0 for k in [
        't2t','gapped_t2t','missassembled','gapped_missassembled','incomplete',
        'gapped_incomplete','no_telomeres','gapped_no_telomeres','discordant','gapped_discordant'
    ])
    if out['total_paths'] is None or out['total_gaps'] is None or out['total_telomeres'] is None or need_categories:
        in_table = False
        header_seen = False
        idx_type = idx_telomeres = idx_gaps = None
        total_paths_calc = total_gaps_calc = total_telos_calc = 0
        cats_calc = {k:0 for k in ['t2t','gapped_t2t','missassembled','gapped_missassembled','incomplete',
                                   'gapped_incomplete','no_telomeres','gapped_no_telomeres','discordant','gapped_discordant']}

        for ln in lines:
            if ln.strip() == '+++ Path Summary Report +++':
                in_table, header_seen = True, False
                continue
            if in_table and not header_seen:
                if ln.lower().startswith('pos\theader\t'):
                    cols = ln.split('\t')
                    idx_telomeres = cols.index('telomeres') if 'telomeres' in cols else None
                    idx_gaps = cols.index('gaps') if 'gaps' in cols else None
                    idx_type = cols.index('type') if 'type' in cols else None
                    header_seen = True
                continue
            if in_table and header_seen:
                if not ln.strip() or ln.startswith('+++') or ln.lower().startswith('two telomeres'):
                    in_table = False
                    continue
                if ln.strip() == '...':  # elided rows
                    continue
                parts = ln.split('\t')
                if len(parts) < 6:
                    continue
                total_paths_calc += 1
                if idx_telomeres is not None and idx_telomeres < len(parts):
                    try: total_telos_calc += int(parts[idx_telomeres])
                    except: pass
                if idx_gaps is not None and idx_gaps < len(parts):
                    try: total_gaps_calc += int(parts[idx_gaps])
                    except: pass
                if idx_type is not None and idx_type < len(parts):
                    key = _to_key(parts[idx_type])
                    if key in cats_calc: cats_calc[key] += 1

        if out['total_paths'] is None: out['total_paths'] = total_paths_calc or out['total_paths']
        if out['total_gaps'] is None: out['total_gaps'] = total_gaps_calc or out['total_gaps']
        if out['total_telomeres'] is None: out['total_telomeres'] = total_telos_calc or out['total_telomeres']
        for k, v in cats_calc.items():
            if (out[k] is None) or (out[k] == 0):
                out[k] = v

    return out

def infer_accession(path: Path) -> str:
    # Prefer parent dir like GCA_*/GCF_*
    for ancestor in [path.parent, *path.parents]:
        name = ancestor.name
        if name.startswith('GCA_') or name.startswith('GCF_'):
            return name
    # Else derive from filename prefix
    stem = path.name
    if '_' in stem:
        return stem.split('_', 1)[0]
    return stem.rsplit('.', 1)[0]

def find_report_files(root: Path) -> List[Path]:
    return sorted(root.rglob('*.telo.report'))

def main():
    ap = argparse.ArgumentParser(description="Parse Teloscope .telo.report files into a CSV.")
    ap.add_argument('--root', default='.', help='Root folder to search (default: current directory)')
    ap.add_argument('--out', default='teloscope_reports_compiled.csv', help='Output CSV path')
    args = ap.parse_args()

    root = Path(args.root).resolve()
    files = find_report_files(root)
    if not files:
        print(f'No .telo.report files found under: {root}', file=sys.stderr)
        sys.exit(2)

    rows = []
    for fp in files:
        try:
            text = fp.read_text(encoding='utf-8', errors='ignore')
            metrics = parse_report_text(text)
            metrics['accession'] = infer_accession(fp)
            rows.append(metrics)
        except Exception as e:
            print(f'[WARN] Failed to parse {fp}: {e}', file=sys.stderr)

    cols = [
        'accession','total_paths','total_gaps','total_telomeres',
        'mean_length','median_length','min_length','max_length',
        'two_telomeres','one_telomere','zero_telomeres',
        't2t','gapped_t2t','missassembled','gapped_missassembled',
        'incomplete','gapped_incomplete','no_telomeres','gapped_no_telomeres',
        'discordant','gapped_discordant'
    ]
    df = pd.DataFrame(rows)
    for c in cols:
        if c not in df.columns:
            df[c] = None
    df = df[cols].sort_values('accession').reset_index(drop=True)
    df.to_csv(args.out, index=False)
    print(f'Wrote: {args.out}   (n={len(df)})')

if __name__ == '__main__':
    main()
