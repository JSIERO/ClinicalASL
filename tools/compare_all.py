"""Where do two runs start to differ? All work NIfTIs (non-QASL) + DICOM series counts. Usage: compare_all.py A B"""
import sys, os, glob, collections
import numpy as np, nibabel as nib, pydicom
A, B = sys.argv[1], sys.argv[2]
fa = {os.path.relpath(p, A + "/work"): p for p in glob.glob(A + "/work/**/*.nii*", recursive=True)}
fb = {os.path.relpath(p, B + "/work"): p for p in glob.glob(B + "/work/**/*.nii*", recursive=True)}
print("nii only in A:", sorted(set(fa) - set(fb))[:10], " only in B:", sorted(set(fb) - set(fa))[:10])
for r in sorted(set(fa) & set(fb)):
    if "QASL" in r: continue
    a = np.asanyarray(nib.load(fa[r]).dataobj).astype(float); b = np.asanyarray(nib.load(fb[r]).dataobj).astype(float)
    if a.shape != b.shape: print(f"{r:<55} SHAPE {a.shape} {b.shape}"); continue
    d = np.abs(np.nan_to_num(a) - np.nan_to_num(b))
    print(f"{r:<55} {'IDENTICAL' if d.max() == 0 else f'max|d| {d.max():.4g} mean|d| {d.mean():.4g} sum|a| {np.abs(np.nan_to_num(a)).mean():.4g}'}")
# QASL fit dirs: compare first-level parameter maps for a sense of fit noise
for r in sorted(set(fa) & set(fb)):
    if "QASL" in r and os.path.basename(r).startswith(("perfusion", "arrival", "mask")) and r.count(os.sep) <= 4:
        a = nib.load(fa[r]).get_fdata(); b = nib.load(fb[r]).get_fdata()
        if a.shape == b.shape:
            d = np.abs(a - b); print(f"{r:<70} {'IDENTICAL' if d.max() == 0 else f'max|d| {d.max():.4g} mean|d| {d.mean():.4g}'}")

def series(root):
    c = collections.Counter()
    for p in glob.glob(root + "/output/**/*", recursive=True):
        if os.path.isfile(p) and not p.endswith((".nii", ".nii.gz", ".png", ".log", ".json", ".txt")):
            try:
                ds = pydicom.dcmread(p, stop_before_pixels=True, specific_tags=["SeriesNumber", "SeriesDescription", "NumberOfFrames"])
                c[(int(ds.get("SeriesNumber", -1)), str(ds.get("SeriesDescription", "")))] += 1
            except Exception:
                pass
    return c
sa, sb = series(A), series(B)
print(f"DICOM files: A {sum(sa.values())} in {len(sa)} series, B {sum(sb.values())} in {len(sb)} series, same list: {sa == sb}")
for k in sorted(set(sa) | set(sb)):
    print(f"  {k[0]:>5} {k[1]:<50} A {sa.get(k, 0):>4} B {sb.get(k, 0):>4}")
print("output files A:", len(glob.glob(A + "/output/**/*", recursive=True)), "B:", len(glob.glob(B + "/output/**/*", recursive=True)))
