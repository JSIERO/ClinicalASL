"""Compare the DICOM output of two ClinicalASL runs file by file (matched by file name): pixel data and every header
element except UIDs, content/instance dates and times. Prints only tag keywords and counts, never values (no
patient data). Usage: compare_dicom.py <outputA> <outputB>"""
import sys, os, glob, collections
import pydicom

IGNORE = {"SOPInstanceUID", "SeriesInstanceUID", "MediaStorageSOPInstanceUID", "ContentDate", "ContentTime",
          "InstanceCreationDate", "InstanceCreationTime", "SoftwareVersions", "Manufacturer"}
a_root, b_root = sys.argv[1], sys.argv[2]
names = lambda r: {os.path.basename(p) for p in glob.glob(os.path.join(r, "*.dcm"))}
A, B = names(a_root), names(b_root)
print(f"DICOM files: A {len(A)}, B {len(B)}, common {len(A & B)}, only A {len(A - B)}, only B {len(B - A)}")
for n in sorted(A ^ B)[:20]:
    print("  only in", "A" if n in A else "B", n)


def flat(ds, prefix=""):
    # {path: value} for all elements, recursing into sequences; UIDs/times in IGNORE skipped
    out = {}
    for el in ds:
        key = prefix + (el.keyword or str(el.tag))
        if el.keyword in IGNORE or el.tag == 0x7FE00010:
            continue
        if el.VR == "SQ":
            for i, item in enumerate(el.value):
                out.update(flat(item, f"{key}[{i}]."))
        else:
            out[key] = el.value
    return out


diff_tags, pixel_diff, n_same = collections.Counter(), 0, 0
for n in sorted(A & B):
    da, db = pydicom.dcmread(os.path.join(a_root, n)), pydicom.dcmread(os.path.join(b_root, n))
    fa, fb = flat(da), flat(db)
    fa.update({"meta." + k: v for k, v in flat(da.file_meta).items()})
    fb.update({"meta." + k: v for k, v in flat(db.file_meta).items()})
    tags = {k for k in fa.keys() | fb.keys() if fa.get(k) != fb.get(k)}
    pix = da.PixelData != db.PixelData
    pixel_diff += pix
    diff_tags.update(t.split("[")[0] if "[" in t else t for t in tags)
    n_same += not tags and not pix
print(f"identical (ignoring UIDs/times): {n_same}/{len(A & B)} files; pixel data differs in {pixel_diff}")
for t, c in diff_tags.most_common(40):
    print(f"  differs x{c:<4} {t}")
