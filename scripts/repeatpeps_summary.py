#!/usr/bin/env python3
"""Resumen del alineamiento contra proteinas de elementos transponibles.

Lee las salidas de run_repeatpeps.sh y produce:

  repeatpeps/repeatpeps_por_conjunto.tsv
      fraccion de secuencias con hit (e <= 1e-10) en cada uno de los siete
      conjuntos, con la sensibilidad a la cobertura de consulta (>= 30, 50, 70 %).
  repeatpeps/repeatpeps_novel_loci.tsv
      los 790 loci posicionalmente noveles de Helixer (helixer_novel_loci.tsv)
      con su categoria de la Tabla 6 y su mejor hit contra RepeatPeps, si lo hay.
  repeatpeps/additional_file_9.csv
      mejor hit de cada secuencia con hit, en los siete conjuntos, con familia,
      identidad, e-value, bitscore y ambas coberturas.

QUE DEBE VERSE
  La fraccion con hit es la misma en los cuatro brazos, en el combinado y en los
  dos proteomas de referencia de RefSeq (en torno al 1,3 %): es una propiedad de
  los proteomas de mamifero (genes derivados de elementos transponibles, hAT-Ac,
  PiggyBac, Gypsy...), no del metodo ni del enmascaramiento. De los 790 loci
  noveles, 12 tienen hit, y de los 225 sin ortologo, ninguno: los ORF de
  elementos transponibles no explican ese grupo.

Escrito el 3 de septiembre de 2026.
"""
import csv, gzip, os, re
from collections import Counter, OrderedDict

BASE = os.environ.get("BASE", os.path.expanduser("~/llama_annotation_pipeline"))
WORK = os.path.join(BASE, "repeatpeps")
NOVEL = os.environ.get("NOVEL", os.path.join(BASE, "blastp", "helixer_novel_loci.tsv"))
SPROT = os.environ.get("SPROT", os.path.join(BASE, "blastp", "novel_vs_sprot.tsv"))

SETS = OrderedDict([("helixer", "Helixer"), ("liftoff", "Liftoff"), ("miniprot", "miniprot"),
                    ("lifton", "LiftOn"), ("combined", "Combined candidate proteome"),
                    ("vpac_ref", "V. pacos RefSeq (reference)"), ("cdro_ref", "C. dromedarius RefSeq (external)")])

def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")

def n_seqs(path):
    with open(path) as fh:
        return sum(1 for l in fh if l.startswith(">"))

def load_hits(s):
    d = {}
    with open(os.path.join(WORK, f"hits_{s}.tsv")) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            qs, qe, ss, se, ql, sl = (int(x) for x in c[6:12])
            m = re.search(r"#(\S+)", c[1])
            d[c[0]] = dict(subject=c[1].split()[0], family=(m.group(1) if m else ""),
                           pident=float(c[2]), length=int(c[3]), evalue=float(c[4]), bitscore=float(c[5]),
                           qcov=100.0 * (qe - qs + 1) / ql, scov=100.0 * (abs(se - ss) + 1) / sl, qlen=ql)
    return d

# ------------------------------------------------------------------ por conjunto
rows_set = []
hits = {}
for s, label in SETS.items():
    n = n_seqs(os.path.join(WORK, "q", f"{s}.faa"))
    h = load_hits(s); hits[s] = h
    q30 = sum(1 for v in h.values() if v["qcov"] >= 30)
    q50 = sum(1 for v in h.values() if v["qcov"] >= 50)
    q70 = sum(1 for v in h.values() if v["qcov"] >= 70)
    rows_set.append([label, n, len(h), f"{100*len(h)/n:.2f}", q30, q50, q70])
    print(f"{label:34s} n={n:6d}  hit={len(h):4d} ({100*len(h)/n:.2f} %)  qcov>=30/50/70: {q30}/{q50}/{q70}")
with open(os.path.join(WORK, "repeatpeps_por_conjunto.tsv"), "w") as out:
    out.write("set\tsequences\thits_e1e-10\tpct\thits_qcov30\thits_qcov50\thits_qcov70\n")
    for r in rows_set: out.write("\t".join(map(str, r)) + "\n")

# ------------------------------------------------------------------ los 790 loci noveles
novel = [l.rstrip("\n").split("\t") for l in opener(NOVEL)][1:]
sprot = set(l.split("\t")[0] for l in opener(SPROT))
def categoria(r):
    dro, alp = r[6] == "si", r[7] == "si"
    if dro and alp: return "orthologue in both references"
    if dro: return "orthologue in C. dromedarius only"
    if alp: return "orthologue in V. pacos only"
    return "no camelid orthologue, Swiss-Prot hit" if r[0] in sprot else "no orthologue at all"
hx = hits["helixer"]
tot, hit = Counter(), Counter()
with open(os.path.join(WORK, "repeatpeps_novel_loci.tsv"), "w") as out:
    out.write("mrna_id\tlocus_length\tcategory\trepeat_hit\tsubject\tfamily\tpident\tevalue\tqcov\n")
    for r in novel:
        k = categoria(r); tot[k] += 1
        h = hx.get(r[0])
        if h: hit[k] += 1
        out.write("\t".join([r[0], r[4], k, "yes" if h else "no",
                             h["subject"] if h else "", h["family"] if h else "",
                             f"{h['pident']:.1f}" if h else "", f"{h['evalue']:.2e}" if h else "",
                             f"{h['qcov']:.0f}" if h else ""]) + "\n")
print("\nLoci noveles de Helixer (Tabla 6):")
for k in ["orthologue in both references", "orthologue in C. dromedarius only", "orthologue in V. pacos only",
          "no camelid orthologue, Swiss-Prot hit", "no orthologue at all"]:
    print(f"  {k:40s} n={tot[k]:4d}  repeat hit={hit[k]:3d}")
print(f"  {'TOTAL':40s} n={sum(tot.values()):4d}  repeat hit={sum(hit.values()):3d}")

# ------------------------------------------------------------------ additional file 9
with open(os.path.join(WORK, "additional_file_9.csv"), "w", newline="") as out:
    w = csv.writer(out)
    w.writerow(["set", "query_id", "subject_id", "repeat_family", "pident", "aln_length", "evalue", "bitscore", "query_cov_pct", "subject_cov_pct"])
    for s, label in SETS.items():
        for q, v in sorted(hits[s].items()):
            w.writerow([label, q, v["subject"], v["family"], f"{v['pident']:.1f}", v["length"], f"{v['evalue']:.2e}",
                        f"{v['bitscore']:.1f}", f"{v['qcov']:.1f}", f"{v['scov']:.1f}"])
print("\nescrito:", os.path.join(WORK, "additional_file_9.csv"))
