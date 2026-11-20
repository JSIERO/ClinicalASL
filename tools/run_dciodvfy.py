"""Run dciodvfy (per file) and dcentvfy (all files) on ClinicalASL DICOM output; summarise per series.
Values inside <...> or after '=' are redacted in the summary (no patient tags in chat). Full logs go to <logdir>.
Usage: run_dciodvfy.py <output dir> <logdir>"""
import sys, os, re, subprocess, collections, glob
import pydicom
TOOLS = r"I:\ProtonDrive\My files\Documents\studies - IMAGR\dicom3tools_clunie"
out, logdir = sys.argv[1], sys.argv[2]
os.makedirs(logdir, exist_ok=True)
files = [p for p in glob.glob(os.path.join(out, "**", "*"), recursive=True)
         if os.path.isfile(p) and not p.lower().endswith((".log", ".json", ".txt", ".png", ".nii", ".gz"))]
def red(s):
    # keep attribute paths <...> and Module=/Element= names; hide any value after '=' (could be patient data)
    s = re.sub(r"(?<!Module)(?<!Element)(=\s*)<[^>]*>", r"\1..", s)
    s = re.sub(r"(Value dubious for this VR \[\w+\]\s*=\s*).*$", r"\1..", s)
    s = re.sub(r"(for file\s*).*$", r"\1..", s)
    return s.strip()
per_series = collections.defaultdict(lambda: {"n": 0, "err": collections.Counter(), "warn": collections.Counter()})
full = []
for p in sorted(files):
    try:
        ds = pydicom.dcmread(p, stop_before_pixels=True)
    except Exception:
        continue
    key = (int(ds.get("SeriesNumber", -1)), str(ds.get("SeriesDescription", "")), str(ds.get("SOPClassUID", "")))
    r = subprocess.run([os.path.join(TOOLS, "dciodvfy.exe"), "-new", p], capture_output=True, text=True, errors="replace")
    txt = (r.stdout or "") + (r.stderr or "")
    full.append(f"===== {os.path.relpath(p, out)}\n{txt}")
    s = per_series[key]; s["n"] += 1
    for line in txt.splitlines():
        line = line.strip()
        if line.startswith("Error"): s["err"][red(line)] += 1
        elif line.startswith("Warning"): s["warn"][red(line)] += 1
open(os.path.join(logdir, "dciodvfy_full.txt"), "w", encoding="utf-8").write("\n".join(full))
r = subprocess.run([os.path.join(TOOLS, "dcentvfy.exe")] + sorted(files), capture_output=True, text=True, errors="replace")
ent = (r.stdout or "") + (r.stderr or "")
open(os.path.join(logdir, "dcentvfy_full.txt"), "w", encoding="utf-8").write(ent)

tot_e = tot_w = 0
for key in sorted(per_series):
    s = per_series[key]; e, w = sum(s["err"].values()), sum(s["warn"].values()); tot_e += e; tot_w += w
    print(f"\n## series {key[0]} '{key[1]}' SOP {key[2]}  files {s['n']}  errors {e}  warnings {w}")
    for m, c in s["err"].most_common(): print(f"   E x{c:<4} {m[:230]}")
    for m, c in s["warn"].most_common(): print(f"   W x{c:<4} {m[:230]}")
print(f"\nTOTAL files {len(files)} errors {tot_e} warnings {tot_w}")
print("\n## dcentvfy (entity consistency across all files):")
ents = collections.Counter(red(l) for l in ent.splitlines() if l.strip())
for m, c in ents.most_common(40): print(f"   x{c:<4} {m[:230]}")
