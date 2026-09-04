#!/usr/bin/env python3
"""
Recompute BUSCO completeness for the homology arms restricted to the substrate
that the ab initio arm analysed, and write additional file 1.

Why this exists
---------------
Helixer was run through its web interface on the sequences of 25 kb or more (244
scaffolds, 79.6 % of the assembly), while the three homology arms used the whole
assembly. Comparing their BUSCO completeness therefore compares two different
substrates. This script removes that asymmetry without re-running BUSCO: it keeps,
for every BUSCO group, only the matches that fall on the scaffolds the ab initio
arm analysed, and reclassifies the group with BUSCO's own rule.

Why not re-run BUSCO
--------------------
In protein mode the HMMER score and the length criterion are per sequence, so
removing proteins does not change the classification of the ones that remain.
Run with --check, the unrestricted recomputation it prints reproduces the deposited
short summaries exactly; that is the validity control for the restricted figures.

Known limit
-----------
The full table records only the matches BUSCO reported: for a group classified as
complete, additional fragmented candidates are not listed. When the complete match
is removed by the restriction, this script therefore counts the group as missing
where a re-run might report it as fragmented. That affects the split between
fragmented and missing, never the percentage of complete BUSCOs, which is exact.

Inputs
------
--busco-dir     directory holding full_table.<arm>.tsv[.gz] for every arm
--gff3-dir      directory holding <arm>.gff3[.gz] with the mRNA features
--substrate     GFF3 of the ab initio arm; its ##sequence-region lines define the
                substrate (any GFF3 with those lines works)
--total         number of BUSCO groups in the lineage (12594 for artiodactyla_odb12)

Output: one CSV row per arm with the counts in both conditions.
"""
import argparse, csv, gzip, os, sys
from collections import Counter, defaultdict


DISPLAY = {'liftoff': 'Liftoff', 'miniprot': 'miniprot', 'lifton': 'LiftOn',
           'helixer': 'Helixer'}


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')


def find(directory, *stems):
    for stem in stems:
        for ext in ('.gz', ''):
            p = os.path.join(directory, stem + ext)
            if os.path.exists(p):
                return p
    sys.exit(f'not found in {directory}: {stems[0]}')


def substrate(path):
    """Scaffold -> length, from the ##sequence-region lines of a GFF3."""
    seqs = {}
    with opener(path) as fh:
        for line in fh:
            if line.startswith('##sequence-region'):
                f = line.split()
                seqs[f[1]] = int(f[3])
    if not seqs:
        sys.exit(f'no ##sequence-region lines in {path}')
    return seqs


def attr(a, key):
    i = a.find(key)
    if i < 0:
        return ''
    s = a[i + len(key):]
    j = s.find(';')
    return (s[:j] if j >= 0 else s).strip()


def transcript_scaffold(path, feature='mRNA'):
    """Transcript ID -> scaffold, from the mRNA features of a GFF3."""
    m = {}
    with opener(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != feature:
                continue
            tid = attr(f[8], 'ID=')
            if tid:
                m[tid] = f[0]
    return m


def full_table(path):
    """BUSCO id -> [(status, sequence), ...] from a BUSCO full table."""
    rows = defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) >= 2:
                rows[f[0]].append((f[1], f[2] if len(f) > 2 else ''))
    return rows


def classify(rows, total, keep=None):
    """BUSCO's rule: >=2 complete -> duplicated, 1 -> single, else fragmented, else missing."""
    c = Counter()
    for group, matches in rows.items():
        complete = [s for st, s in matches
                    if st in ('Complete', 'Duplicated') and (keep is None or keep(s))]
        fragmented = [s for st, s in matches
                      if st == 'Fragmented' and (keep is None or keep(s))]
        if len(complete) >= 2:
            c['D'] += 1
        elif len(complete) == 1:
            c['S'] += 1
        elif fragmented:
            c['F'] += 1
        else:
            c['M'] += 1
    c['M'] += total - len(rows)          # groups absent from the table are missing
    c['C'] = c['S'] + c['D']
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--busco-dir', required=True)
    ap.add_argument('--gff3-dir', required=True)
    ap.add_argument('--substrate', required=True,
                    help='GFF3 of the ab initio arm; its ##sequence-region lines are the substrate')
    ap.add_argument('--arms', default='liftoff,miniprot,lifton',
                    help='comma-separated homology arms (default: %(default)s)')
    ap.add_argument('--ab-initio-arm', default='helixer',
                    help='name of the ab initio arm, added as a reference row')
    ap.add_argument('--gff3-prefix', default='llama_',
                    help='prefix of the GFF3 file names (default: %(default)s)')
    ap.add_argument('--total', type=int, default=12594,
                    help='BUSCO groups in the lineage (default: %(default)s)')
    ap.add_argument('--check', action='store_true',
                    help='also print the unrestricted recomputation, which must equal '
                         'the deposited short summaries')
    ap.add_argument('--out', default='additional_file_1.csv')
    args = ap.parse_args()

    scaffolds = substrate(args.substrate)
    print(f'substrate: {len(scaffolds)} scaffolds, shortest {min(scaffolds.values()):,} bp, '
          f'{sum(scaffolds.values()):,} bp in total', file=sys.stderr)
    wl = set(scaffolds)

    rows = []
    for arm in args.arms.split(','):
        arm = arm.strip()
        tx = transcript_scaffold(find(args.gff3_dir, f'{args.gff3_prefix}{arm}.gff3'))
        table = full_table(find(args.busco_dir, f'full_table.{arm}.tsv'))
        unknown = {s for ms in table.values() for _, s in ms if s and s not in tx}
        if unknown:
            sys.exit(f'{arm}: {len(unknown)} BUSCO sequences absent from the GFF3, '
                     f'e.g. {sorted(unknown)[:3]}')
        outside = sum(1 for v in tx.values() if v not in wl)
        full = classify(table, args.total)
        rest = classify(table, args.total, keep=lambda s: tx.get(s) in wl)
        if args.check:
            print(f'{arm:9s} unrestricted  C:{100*full["C"]/args.total:.1f}% '
                  f'[S:{100*full["S"]/args.total:.1f}%,D:{100*full["D"]/args.total:.1f}%],'
                  f'F:{100*full["F"]/args.total:.1f}%,M:{100*full["M"]/args.total:.1f}%',
                  file=sys.stderr)
        rows.append(dict(
            arm=DISPLAY.get(arm, arm), mRNA_total=len(tx), mRNA_on_sequences_under_25kb=outside,
            percentage=round(100 * outside / len(tx), 2),
            C_full=full['C'], S_full=full['S'], D_full=full['D'], F_full=full['F'], M_full=full['M'],
            C_restricted=rest['C'], S_restricted=rest['S'], D_restricted=rest['D'],
            F_restricted=rest['F'], M_restricted=rest['M'],
            C_pct_full=round(100 * full['C'] / args.total, 1),
            C_pct_restricted=round(100 * rest['C'] / args.total, 1),
            delta_pp=round(round(100 * full['C'] / args.total, 1)
                           - round(100 * rest['C'] / args.total, 1), 1)))

    if args.ab_initio_arm:
        arm = args.ab_initio_arm
        tx = transcript_scaffold(find(args.gff3_dir, f'{args.gff3_prefix}{arm}.gff3'))
        table = full_table(find(args.busco_dir, f'full_table.{arm}.tsv'))
        c = classify(table, args.total)
        pct = round(100 * c['C'] / args.total, 1)
        rows.append(dict(
            arm=DISPLAY.get(arm, arm), mRNA_total=len(tx), mRNA_on_sequences_under_25kb=0, percentage=0.0,
            C_full=c['C'], S_full=c['S'], D_full=c['D'], F_full=c['F'], M_full=c['M'],
            C_restricted=c['C'], S_restricted=c['S'], D_restricted=c['D'],
            F_restricted=c['F'], M_restricted=c['M'],
            C_pct_full=pct, C_pct_restricted=pct, delta_pp=0.0))

    with open(args.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f'{r["arm"]:9s} {r["C_pct_full"]:5.1f}% -> {r["C_pct_restricted"]:5.1f}% '
              f'({r["delta_pp"]:+.1f})', file=sys.stderr)
    print(f'written: {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
