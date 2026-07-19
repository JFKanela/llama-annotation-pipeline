
# =====================================================================
# run_helixer.py - Helixer (ab initio) branch, external GPU execution
#
# Run on Kaggle (NVIDIA Tesla P100) inside an isolated `uv` virtual environment;
# NOT automated in Snakemake. helixerlite 25.5.27, model
# `vertebrate_v0.3_m_0080.h5` (DOI 10.5281/zenodo.10836346).
#
# IMPORTANT - read before comparing this branch with the homology branches:
# this run used scaffolds >= 10 kb only, NO subsequence overlap, and loaded the
# model by path rather than by lineage. Those choices make the reported
# completeness a lower bound, not a benchmark of the method. See the
# "Helixer run conditions (important)" section of README.md and the
# "Helixer run parameters" section of REPRODUCIBILITY.md.
#
# Citation: Nature Methods (2025), DOI 10.1038/s41592-025-02939-1.
# =====================================================================
import os, sys, time, shutil
os.chdir("/kaggle/working")
import pyfastx
from helixerlite.__main__ import fasta2hdf5
from helixerlite.hybrid_model import HybridModel
from helixerlite.utilities import preds2gff3
from gfftk.gff import gff2dict, dict2gff3
FASTA, GFF = sys.argv[1], sys.argv[2]
CFG = {"url":"https://zenodo.org/records/10836346/files/vertebrate_v0.3_m_0080.h5?download=1",
       "name":"vertebrate_v0.3_m_0080.h5","length":213840}
BATCH = 16; CHUNK_BP = 300_000_000
MODEL = "/kaggle/working/"+CFG["name"]
CKDIR = "/kaggle/working/hxchunks"; TMP="/tmp/hxbig"
os.makedirs(CKDIR, exist_ok=True); os.makedirs(TMP, exist_ok=True)
if not os.path.isfile(MODEL):
    os.system('curl -L -C - --retry 10 -o "%s" "%s"' % (MODEL, CFG["url"]))
print("Modelo:", round(os.path.getsize(MODEL)/1e6,1), "MB", flush=True)

# trocear por bp acumulados
chunks=[]; cur=[]; curbp=0; idx=0
for name, seq in pyfastx.Fasta(FASTA, build_index=False):
    cur.append((name,seq)); curbp+=len(seq)
    if curbp>=CHUNK_BP:
        p=f"{TMP}/chunk{idx}.fasta"; open(p,"w").write("".join(f">{n}\n{s}\n" for n,s in cur))
        chunks.append((idx,p)); cur=[]; curbp=0; idx+=1
if cur:
    p=f"{TMP}/chunk{idx}.fasta"; open(p,"w").write("".join(f">{n}\n{s}\n" for n,s in cur)); chunks.append((idx,p))
print(f">> {len(chunks)} chunks de ~{CHUNK_BP/1e6:.0f} Mb", flush=True)

t0=time.time(); gffs=[]
for i,cf in chunks:
    gff_i=f"{CKDIR}/chunk{i}.gff3"
    if os.path.isfile(gff_i) and os.path.getsize(gff_i)>0:
        print(f">> chunk {i+1}/{len(chunks)}: ya hecho, salto", flush=True); gffs.append(gff_i); continue
    print(f">> chunk {i+1}/{len(chunks)}: prediciendo (batch={BATCH}, N/M abajo)", flush=True)
    h5=f"{TMP}/c{i}.h5"; pred=f"{TMP}/c{i}.pred.h5"
    fasta2hdf5(cf, h5, species=f"c{i}", subseqlen=CFG["length"])
    HybridModel(["--load-model-path",MODEL,"--test-data",h5,"--val-test-batch-size",str(BATCH),
                 "-v","--prediction-output-path",pred,"--cpus","4"]).run()
    preds2gff3(h5, pred, gff_i, peak_threshold=0.8, min_coding_length=60)
    for f in (h5,pred,cf):
        if os.path.isfile(f): os.remove(f)
    gffs.append(gff_i)
    print(f">> chunk {i+1}/{len(chunks)} OK ({(time.time()-t0)/60:.1f} min acum) | disco /tmp libre: {shutil.disk_usage('/tmp').free/1e9:.1f} GB", flush=True)

# concatenar y limpiar
raw="/tmp/all_helixer.gff3"
with open(raw,"w") as o:
    o.write("##gff-version 3\n")
    for g in gffs:
        for ln in open(g):
            if not ln.startswith("#"): o.write(ln)
Genes=gff2dict(raw, FASTA); dict2gff3(Genes, output=GFF)
print("Tiempo total: %.1f min | genes: %d" % ((time.time()-t0)/60, len(Genes)), flush=True)
assert os.path.getsize(GFF)>0
print("GFF3 OK:", round(os.path.getsize(GFF)/1e6,2), "MB", flush=True)
