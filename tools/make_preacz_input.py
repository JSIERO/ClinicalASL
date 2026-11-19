"""Copy ASLtest to a pre-ACZ-only input: skip every series whose SeriesDescription contains 'postACZ'.
Prints only series number/description counts (no patient tags)."""
import os, shutil, sys, collections, pydicom
src, dst = sys.argv[1], sys.argv[2]
kept, skipped = collections.Counter(), collections.Counter()
for root, _, files in os.walk(src):
    for f in files:
        p = os.path.join(root, f)
        try:
            ds = pydicom.dcmread(p, stop_before_pixels=True, force=True)
            key = f"{getattr(ds,'SeriesNumber','?')} {getattr(ds,'SeriesDescription','')}"
        except Exception:
            key = "unreadable"
        if 'postacz' in key.lower():
            skipped[key] += 1; continue
        out = os.path.join(dst, os.path.relpath(p, src))
        os.makedirs(os.path.dirname(out), exist_ok=True)
        shutil.copy2(p, out); kept[key] += 1
print("KEPT:", dict(kept)); print("SKIPPED:", dict(skipped))
