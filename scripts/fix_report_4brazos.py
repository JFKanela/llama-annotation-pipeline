#!/usr/bin/env python3
"""Recomputa la tabla comparativa CORRECTA desde los ficheros ya generados.

No recalcula nada. Arregla los bugs del informe:
  1. Conteo de genes: se lee del GFF directamente (no del texto de AGAT).
  2. gffcompare: se derivan exones de las CDS de miniprot antes de comparar.
  3. LiftOn: usa el fichero _v2 (nativo NCBI) y su BUSCO con nombre .specific.
  4. Incluye Helixer como cuarto brazo (ab initio).

Uso:  python fix_report.py     (desde ~/llama_annotation_pipeline)
"""
import os, re, subprocess, shutil, glob

R = "results"
# GFF por metodo
ANNOT = {
    "liftoff":  f"{R}/liftoff/llama_liftoff.gff3",
    "miniprot": f"{R}/miniprot/llama_miniprot.gff3",
    "lifton":   f"{R}/lifton/llama_lifton_v2.gff3",
    "helixer":  f"{R}/helixer/llama_helixer.gff3",
}
# proteoma (.faa) por metodo -> nº real de proteinas del BUSCO
PROT = {
    "liftoff":  f"{R}/proteomes/llama_liftoff.faa",
    "miniprot": f"{R}/proteomes/llama_miniprot.faa",
    "lifton":   f"{R}/proteomes/llama_lifton_v2.faa",
    "helixer":  f"{R}/proteomes/llama_helixer.faa",
}

def count_feat(gff, feat):
    if not os.path.isfile(gff): return None
    n = 0
    with open(gff) as f:
        for ln in f:
            if ln.startswith("#"): continue
            c = ln.split("\t")
            if len(c) > 2 and c[2] == feat: n += 1
    return n

def count_prot(faa):
    if not os.path.isfile(faa): return None
    return sum(1 for l in open(faa) if l.startswith(">"))

def find_summary(name):
    """Busca el short_summary en sus posibles ubicaciones/nombres."""
    cands = [f"{R}/busco/{name}/short_summary.txt",
             *glob.glob(f"{R}/busco/{name}/short_summary.specific.*.txt"),
             f"{R}/busco/{name}/run_artiodactyla_odb12/short_summary.txt"]
    for p in cands:
        if os.path.isfile(p): return p
    return None

def busco(name):
    p = find_summary(name)
    if not p: return {}
    m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)", open(p).read())
    return dict(zip("C S D F M n".split(), m.groups())) if m else {}

# ---- 1. Conteos reales ----
print("="*70)
print("CONTEO REAL (genes desde GFF, proteinas desde .faa)")
print("="*70)
prot_counts = {}
for meth in ANNOT:
    g  = count_feat(ANNOT[meth], "gene")
    mr = count_feat(ANNOT[meth], "mRNA")
    pr = count_prot(PROT[meth]); prot_counts[meth] = pr
    print(f"{meth:10s} | genes={g} | mRNA={mr} | proteinas(.faa)={pr}")

# ---- 2. Fix gffcompare (exones de miniprot) + concordancia ----
def add_exons(src, dst):
    with open(src) as fi, open(dst, "w") as fo:
        for ln in fi:
            if ln.startswith("#"): fo.write(ln); continue
            c = ln.rstrip("\n").split("\t")
            fo.write(ln)
            if len(c) > 2 and c[2] == "CDS":
                c2 = c[:]; c2[2] = "exon"; fo.write("\t".join(c2) + "\n")

gffcmp = {}
if shutil.which("gffcompare") and os.path.isfile(ANNOT["miniprot"]) and os.path.isfile(ANNOT["liftoff"]):
    os.makedirs(f"{R}/gffcompare", exist_ok=True)
    mp_fix = f"{R}/gffcompare/miniprot.withexon.gff3"
    add_exons(ANNOT["miniprot"], mp_fix)
    subprocess.run(["gffcompare","-r",ANNOT["liftoff"],"-o",f"{R}/gffcompare/cmp_fixed",mp_fix],
                   capture_output=True, text=True)
    st = f"{R}/gffcompare/cmp_fixed.stats"
    if os.path.isfile(st):
        for ln in open(st):
            m = re.search(r"(Transcript|Locus) level:\s+([\d.]+)\s+\|\s+([\d.]+)", ln)
            if m: gffcmp[m.group(1).lower()] = (m.group(2), m.group(3))

# ---- 3. unmapped Liftoff ----
unm = f"{R}/liftoff/unmapped.txt"
n_unm = sum(1 for l in open(unm) if l.strip() and not l.startswith("#")) if os.path.isfile(unm) else "NA"

# ---- 4. Tabla final ----
print("\n" + "="*70)
print("TABLA COMPARATIVA (4 brazos + techo alpaca + baseline chaku)")
print("="*70)
print(f"{'proteoma':11s} {'C%':7s} {'Unico%':8s} {'Dup%':6s} {'Frag%':6s} {'Miss%':6s} {'proteinas':10s}")
order = ["liftoff","miniprot","lifton","helixer","alpaca_ref","chaku_v1"]
for name in order:
    b = busco(name)
    if not b: 
        print(f"{name:11s}  (sin BUSCO)")
        continue
    loci = prot_counts.get(name, "-")
    print(f"{name:11s} {b['C']:7s} {b['S']:8s} {b['D']:6s} {b['F']:6s} {b['M']:6s} {str(loci):10s}")

if gffcmp:
    print("\nConcordancia A(Liftoff) vs B(miniprot), con exones derivados:")
    for lvl,(sn,pr) in gffcmp.items():
        print(f"  {lvl:11s} sensibilidad={sn}%  precision={pr}%")
print(f"\nGenes de referencia NO transferidos por Liftoff (unmapped): {n_unm}")
print("\nNOTA: helixer (ab initio) corrio sobre los scaffolds >=10 kb: 3.640 scaffolds")
print("      (0,34 % del total) que reunen 1.915.763.599 bp, el 81,46 % de la secuencia")
print("      ensamblada. Los brazos de homologia, sobre el ensamblado completo.")
