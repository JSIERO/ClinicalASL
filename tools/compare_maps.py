"""Compare CBF/CVR (and AAT/ATA) NIfTI maps between two ClinicalASL runs. Usage: compare_maps.py <dirA> <dirB>"""
import sys, os, glob, re
import numpy as np, nibabel as nib

a_root, b_root = sys.argv[1], sys.argv[2]
pat = re.compile(r"(CBF|CVR|AAT|ATA)", re.I)

def rel_files(root):
    out = {}
    for p in glob.glob(os.path.join(root, "**", "*.nii*"), recursive=True):
        r = os.path.relpath(p, root)
        if pat.search(os.path.basename(r)) and "QASL" not in r:
            out[r] = p
    return out

A, B = rel_files(a_root), rel_files(b_root)
print(f"maps: A {len(A)}, B {len(B)}, common {len(set(A) & set(B))}")
for r in sorted(set(A) ^ set(B)):
    print("  only in", "A" if r in A else "B", r)
print(f"{'file':<60} {'shape':>14} {'voxels':>8} {'identical%':>10} {'max|d|':>10} {'mean|d|':>10} {'meanA':>9} {'meanB':>9} {'r':>8}")
for r in sorted(set(A) & set(B)):
    a = nib.load(A[r]).get_fdata(); b = nib.load(B[r]).get_fdata()
    if a.shape != b.shape:
        print(f"{r:<60} SHAPE {a.shape} vs {b.shape}"); continue
    m = np.isfinite(a) & np.isfinite(b) & ((a != 0) | (b != 0))
    n = int(m.sum())
    if n == 0:
        print(f"{r:<60} {str(a.shape):>14} {0:>8} (empty)"); continue
    d = np.abs(a[m] - b[m])
    ident = 100.0 * np.mean(d == 0)
    rr = np.corrcoef(a[m], b[m])[0, 1] if np.std(a[m]) > 0 and np.std(b[m]) > 0 else float("nan")
    print(f"{r:<60} {str(a.shape):>14} {n:>8} {ident:>10.2f} {d.max():>10.4g} {d.mean():>10.4g} {a[m].mean():>9.4g} {b[m].mean():>9.4g} {rr:>8.5f}")
