"""Check every reconstructed initial row and optionally the original builders.

The default check uses only the standard library. --compare-original also
compares selected complete problems with check_cyclic.py/check_code.py and
therefore requires NumPy, SciPy and SymPy, as the original search did.
"""

from validation import require
from validation import save_summary
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse, csv, gzip, json, math, time
from export_all_matrices import sources, parameter_layout, decode_certificate
from noncyclic import code_generated

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'matrices'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compare-original', action='store_true')
    args = parser.parse_args()
    start = time.monotonic()
    original = {item[0]: item for item in sources() if item[4] != 'trivial'}
    selected = {
        'cyclic_n8_0001', 'cyclic_n8_0016', 'cyclic_n9_0001',
        'cyclic_n9_0291', 'noncyclic_n8_0000', 'noncyclic_n9_0004',
        'noncyclic_n9_0008', 'extension_n8_0000', 'extension_n8_0007',
        'extension_n8_0010', 'extension_n8_0646',
    }
    extension_items = [v for v in original.values() if v[4] == 'extension']
    selected.add(max(extension_items, key=lambda v: v[-1]['steps'])[0])
    selected.add(max(extension_items, key=lambda v: len(parameter_layout(v[3], v[5], v[6], v[7])))[0])
    if args.compare_original:
        from check_cyclic import build_problem
        from check_code import build_code_problem
    with (OUT/'initial_lp_index.csv').open(newline='') as f:
        csv_rows = list(csv.DictReader(f))
    require((len(csv_rows) == len(original)), 'verify_initial_lp_export.py:38: verification check failed')
    seen_ids, counts, compared = set(), Counter(), []
    with gzip.open(OUT/'initial_lp_matrices.jsonl.gz', 'rt', encoding='utf-8') as f:
        for line_no, line in enumerate(f, 1):
            rec = json.loads(line)
            ident = rec['id']
            require((ident not in seen_ids), 'verify_initial_lp_export.py:44: verification check failed')
            seen_ids.add(ident)
            _, name, key, n, phase, mod, gens, extra, result = original[ident]
            require(((rec['source_file'], rec['source_entry'], rec['rank'], rec['phase'],
                    rec['lattice_modulus'], rec['lattice_generator_rows'], rec['extra_vector']) == (
                    name, key, n, phase, mod, gens, extra)), 'verify_initial_lp_export.py:47: verification check failed')
            layout = parameter_layout(n, mod, gens, extra)
            require((rec['gram_parameters'] == layout), 'verify_initial_lp_export.py:51: verification check failed')
            require((rec['column_order'] == [p['name'] for p in layout]+['delta']), 'verify_initial_lp_export.py:52: verification check failed')
            m = len(layout)
            pairpos = [(i-1, j-1, k) for k, p in enumerate(layout)
                       for i, j in p['coordinate_pairs']]
            code = set(code_generated(gens, mod))

            def features(b):
                fs = [0]*m
                for i, j, k in pairpos:
                    fs[k] += 2*b[i]*b[j]
                return fs

            A, bs = rec['A_ineq'], rec['b_ineq']
            nums, dens = rec['inequality_vector_numerators'], rec['inequality_vector_denominators']
            require((len(A) == len(bs) == len(nums) == len(dens)), 'verify_initial_lp_export.py:66: verification check failed')
            unique_rows = set()
            for j, (row, rhs, v, den) in enumerate(zip(A, bs, nums, dens)):
                require((len(row) == m+1), 'verify_initial_lp_export.py:69: verification check failed')
                if v is None:
                    require((den is None and j == len(A)-1), 'verify_initial_lp_export.py:71: verification check failed')
                    require((row == [0]*m+[1] and rhs == 1), 'verify_initial_lp_export.py:72: verification check failed')
                else:
                    require((len(v) == n and den > 0 and math.gcd(den, *v) == 1 and any(v)), 'verify_initial_lp_export.py:74: verification check failed')
                    require((all(x*mod % den == 0 for x in v)), 'verify_initial_lp_export.py:75: verification check failed')
                    require((tuple((x*mod//den) % mod for x in v) in code), 'verify_initial_lp_export.py:76: verification check failed')
                    norm = sum(x*x for x in v)
                    require((not (norm == den*den and sum(bool(x) for x in v) == 1)), 'verify_initial_lp_export.py:78: verification check failed')
                    if extra is not None:
                        require((any(F(x, den) != F(y) for x, y in zip(v, extra))), 'verify_initial_lp_export.py:80: verification check failed')
                        require((any(F(x, den) != -F(y) for x, y in zip(v, extra))), 'verify_initial_lp_export.py:81: verification check failed')
                    require((row == [-x for x in features(v)]+[den*den]), 'verify_initial_lp_export.py:82: verification check failed')
                    require((rhs == norm-den*den), 'verify_initial_lp_export.py:83: verification check failed')
                    require(((tuple(row), rhs) not in unique_rows), 'verify_initial_lp_export.py:84: verification check failed')
                unique_rows.add((tuple(row), rhs))
            E, hs = rec['E_eq'], rec['h_eq']
            if extra is None:
                require((E == hs == []), 'verify_initial_lp_export.py:88: verification check failed')
            else:
                x = list(map(F, extra))
                den = math.lcm(*(v.denominator for v in x))
                v = [int(den*y) for y in x]
                require((E == [features(v)+[0]]), 'verify_initial_lp_export.py:93: verification check failed')
                require((hs == [den*den-sum(y*y for y in v)]), 'verify_initial_lp_export.py:94: verification check failed')
            if result['status'] == 'excluded':
                c = decode_certificate(result, phase)
                flags = [(tuple(a), b) in unique_rows for a, b in zip(c['A_ineq'], c['b_ineq'])]
                require((rec['retained_certificate_rows_present_initially'] == flags), 'verify_initial_lp_export.py:98: verification check failed')
                if result['steps'] == 1:
                    require((all(flags)), 'verify_initial_lp_export.py:100: verification check failed')
                counts['saved_certificate_rows_from_later_cuts'] += sum(not flag for flag in flags)
            index_row = csv_rows[line_no-1]
            require(((index_row['id'], int(index_row['line_number']), int(index_row['columns']),
                    int(index_row['inequality_rows']), int(index_row['equality_rows'])) == (
                    ident, line_no, m+1, len(A), len(E))), 'verify_initial_lp_export.py:103: verification check failed')
            if args.compare_original and ident in selected:
                if phase == 'cyclic':
                    labels, old_A, old_b, old_vectors, _ = build_problem(n, mod, gens[0])
                    old_E, old_h = [], []
                else:
                    _, _, labels, _, old_A, old_b, old_vectors, old_E, old_h, _, _ = build_code_problem(n, mod, gens, extra)
                require((len(labels) == m), 'verify_initial_lp_export.py:112: verification check failed')
                require((A == [list(row) for row in old_A] and bs == old_b), 'verify_initial_lp_export.py:113: verification check failed')
                require((E == [list(row) for row in old_E] and hs == old_h), 'verify_initial_lp_export.py:114: verification check failed')
                vectors = [None if v is None else [str(F(x, d)) for x in v] for v, d in zip(nums, dens)]
                require((vectors == old_vectors), 'verify_initial_lp_export.py:116: verification check failed')
                compared.append(ident)
                print('Matches original builder:', ident, flush=True)
            counts['problems'] += 1
            counts['inequality_rows'] += len(A)
            counts['equality_rows'] += len(E)
    require((seen_ids == set(original)), 'verify_initial_lp_export.py:122: verification check failed')
    if args.compare_original:
        require((set(compared) == selected), 'verify_initial_lp_export.py:124: verification check failed')
    saved_summary = json.loads((OUT/'initial_lp_summary.json').read_text())
    require((dict(counts) == saved_summary['statistics']), 'verify_initial_lp_export.py:126: verification check failed')
    summary = {'status': 'passed', 'checks': dict(counts),
               'whole_problems_compared_with_original_builders': compared,
               'seconds': round(time.monotonic()-start, 3),
               'scope': 'Every exported initial constraint checked exactly; optional whole-problem comparisons check original row order, deduplication, and representative vectors on the listed examples.'}
    save_summary(OUT/'initial_lp_verification_summary.json', json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
