"""Reconstruct the initial LP matrices, not unrecorded solver iterations.

Uses integer arithmetic and the Python standard library. The row order,
deduplication, and representative vectors match the original builders in
check_cyclic.py and check_code.py. Later cutting planes were not logged;
the retained exclusion certificates contain the cuts needed for their proofs.
"""

from validation import require
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import csv, gzip, itertools, json, math, time
from export_all_matrices import sources, parameter_layout, decode_certificate, encode
from noncyclic import code_generated

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'matrices'


def reduced_vector(v):
    v = tuple(F(x) for x in v)
    den = math.lcm(*(x.denominator for x in v))
    return tuple(int(den*x) for x in v), den


def build_initial(item):
    ident, name, key, n, phase, mod, gens, extra, result = item
    layout = parameter_layout(n, mod, gens, extra)
    m = len(layout)
    pairpos = [(i-1, j-1, k) for k, p in enumerate(layout)
               for i, j in p['coordinate_pairs']]
    required_extra = None if extra is None else reduced_vector(extra)
    extra_negative = (None if extra is None else
                      (tuple(-x for x in required_extra[0]), required_extra[1]))
    rows, rhs, numerators, denominators, seen = [], [], [], [], set()

    def features(b):
        fs = [0]*m
        for i, j, k in pairpos:
            fs[k] += 2*b[i]*b[j]
        return fs

    def add_reduced(b, den):
        if not any(b):
            return
        norm = sum(x*x for x in b)
        if norm == den*den and sum(bool(x) for x in b) == 1:
            return
        if (b, den) in (required_extra, extra_negative):
            return
        fs = features(b)
        row = tuple([-x for x in fs] + [den*den])
        val = norm-den*den
        if (row, val) in seen:
            return
        seen.add((row, val))
        rows.append(row)
        rhs.append(val)
        numerators.append(b)
        denominators.append(den)

    for i in range(n):
        for j in range(i+1, n):
            for sign in (-1, 1):
                v = [0]*n
                v[i], v[j] = 1, sign
                add_reduced(tuple(v), 1)
    if phase == 'cyclic':
        cosets = [tuple(k*x % mod for x in gens[0]) for k in range(1, mod)]
    else:
        cosets = [c for c in code_generated(gens, mod) if any(c)]
    for c in cosets:
        divisor = math.gcd(mod, *c)
        den = mod//divisor
        choices = [(r//divisor, (r-mod)//divisor) if r else (0,) for r in c]
        for b in itertools.product(*choices):
            add_reduced(b, den)
    # The original builders append this row even if a duplicate is present.
    rows.append(tuple([0]*m + [1]))
    rhs.append(1)
    numerators.append(None)
    denominators.append(None)
    eqs, eqrhs = [], []
    if required_extra is not None:
        b, den = required_extra
        eqs.append(tuple(features(b) + [0]))
        eqrhs.append(den*den-sum(x*x for x in b))
    rec = {
        'id': ident, 'source_file': name, 'source_entry': key,
        'origin': 'reconstructed initial problem; later solver states were not retained',
        'rank': n, 'phase': phase, 'lattice_modulus': mod,
        'lattice_generator_rows': gens, 'extra_vector': extra,
        'gram_parameters': layout,
        'column_order': [p['name'] for p in layout]+['delta'],
        'objective': 'maximize delta; all variables are unrestricted real numbers',
        'A_ineq': rows, 'b_ineq': rhs,
        'E_eq': eqs, 'h_eq': eqrhs,
        'inequality_vector_numerators': numerators,
        'inequality_vector_denominators': denominators,
    }
    if result['status'] == 'excluded':
        cert = decode_certificate(result, phase)
        present = [(tuple(row), b) in seen or
                   (tuple(row), b) == (tuple([0]*m+[1]), 1)
                   for row, b in zip(cert['A_ineq'], cert['b_ineq'])]
        rec['retained_certificate_rows_present_initially'] = present
        for row, b in zip(cert['E_eq'], cert['h_eq']):
            require((any(tuple(row) == e and b == h for e, h in zip(eqs, eqrhs))), 'export_initial_lp.py:107: verification check failed')
        if result['steps'] == 1:
            require((all(present)), lambda: (ident))
    return rec


def main():
    start = time.monotonic()
    OUT.mkdir(exist_ok=True)
    entries = [item for item in sources() if item[4] != 'trivial']
    counts = Counter()
    index = []
    target = OUT/'initial_lp_matrices.jsonl.gz'
    with gzip.open(target, 'wt', encoding='utf-8', compresslevel=6) as out:
        for line_number, item in enumerate(entries, 1):
            rec = build_initial(item)
            out.write(json.dumps(encode(rec), separators=(',', ':'))+'\n')
            row = {'id': rec['id'], 'line_number': line_number,
                   'columns': len(rec['column_order']),
                   'inequality_rows': len(rec['A_ineq']),
                   'equality_rows': len(rec['E_eq'])}
            index.append(row)
            counts['problems'] += 1
            counts['inequality_rows'] += row['inequality_rows']
            counts['equality_rows'] += row['equality_rows']
            counts['saved_certificate_rows_from_later_cuts'] += sum(
                not x for x in rec.get('retained_certificate_rows_present_initially', []))
            if line_number % 100 == 0:
                print('Reconstructed', line_number, 'of', len(entries), 'initial problems', flush=True)
    require((len(entries) == 964), 'export_initial_lp.py:136: verification check failed')
    with (OUT/'initial_lp_index.csv').open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=list(index[0]))
        writer.writeheader()
        writer.writerows(index)
    summary = {'statistics': dict(counts),
               'compressed_bytes': target.stat().st_size,
               'seconds': round(time.monotonic()-start, 3),
               'provenance': 'Reconstruction from the original builder rules, not a solver-state log.'}
    (OUT/'initial_lp_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
