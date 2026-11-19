---
name: dicom-export
description: ClinicalASL step 12 DICOM/PNG output specialist. Use when changing or checking the grayscale, PALETTE COLOR or PNG Secondary Capture DICOMs (tags, series numbers, scaling, RWVM, LUTs), and to run the fast export-only test plus dciodvfy/dcentvfy validation.
tools: Read, Grep, Glob, Edit, Bash, PowerShell
model: inherit
---

Files: `utils/dicom_common.py` (scaling, UIDs + patient/study defaults, common metadata, source reference,
multiframe frame selection, single-frame template selection, timing cleanup, series numbers, file names, one log
line per series), `utils/save_data_dicom_grayscale.py` and `utils/save_data_dicom_color.py` (only what differs:
MONOCHROME2 + Rescale/W-L/RWVM vs PALETTE COLOR LUT), `utils/save_png_to_dicom.py`, `utils/save_figure_to_png.py`;
orchestration in `asl_save_results_cbfaatcvr.py` (`save_nifti_and_dicom`, `save_png`). Export failures raise.

## Data rules
DICOM headers hold real patient data. Never print PatientName/ID/BirthDate/AccessionNumber or dump a header.
Inspect with pydicom on named, non-identifying tags only. `tools/run_dciodvfy.py` redacts values;
`tools/compare_dicom.py` prints only tag keywords. Don't log patient/study identifiers from exporters.

## What step 12 writes (outputdir)
| Series | SeriesNumber | SOP class (single-frame / multiframe) | File name |
|---|---|---|---|
| Map grayscale CBF/AAT/ATA per ctx | source SN + 10 + i (i = 1 CBF, 2 AAT, 3 ATA) | MR Image / Enhanced MR | `ASL_CBF_preACZ_<SN>_<inst>.dcm` |
| Map colour | grayscale + 100 | SC (1.2.840.10008.5.1.4.1.1.7) / Multi-frame True Color SC (...7.4) | `ASL_CBF_preACZ_COLOR_<SN>_<inst>.dcm` |
| CVR grayscale / colour | 888 / 988 (hard-coded) | as above | `ASL_CVR_888_<inst>.dcm` |
| PNG montages | 999, one SeriesInstanceUID | Secondary Capture | `ASL_CBF_preACZ_999_PNG.dcm`, `ASL_CVR_999_PNG.dcm` |
Single mode: no study tag in names (`ASL CBF`), no CVR. Expected on ASLtest (16 slices, single-frame PACS):
CVR mode 231 files / 15 series; single mode 99 files / 7 series.

## Encoding rules
- Grayscale: uint16 scaled by 65535/nanmax (CVR: int16, 32767/abs max, PixelRepresentation 1). Per-image
  RescaleSlope = 1/scale, intercept 0; W/L from the display range; RealWorldValueMappingSequence rebuilt from the
  map's own scaling (First/Last value mapped follow PixelRepresentation, LUTLabel = type, units UCUM text).
  Inherited Philips RWVM, TriggerTime and VelocityEncodingDirection are removed.
- Colour: PALETTE COLOR, 65536-entry 16-bit LUT built from the same colormap and range as the PNG (viridis CBF/ATA,
  devon AAT, vik CVR from `colormaps/*.mat`); index 0 black = outside the mask. CVR is offset to uint16 (+32768).
  No Rescale, no W/L, no RWVM on colour. ConversionType WSD.
- PNG SC: 8-bit RGB, PatientOrientation L\P, Laterality copied, map name in ImageComments, one SeriesDescription
  "ASL results (ClinicalASL, NOT CE MARKED)". Instance 1 = CVR (cvr mode), then CBF, AAT, ATA x contexts.
- Every series: SeriesDescription `<name> (ClinicalASL, NOT CE MARKED)`, ImageType DERIVED\SECONDARY\QUANTITATIVE,
  DerivationCodeSequence DCM 126302, ReferencedSeriesSequence to the source, UIDs under root 1.3.6.1.4.1.54321.1.1
  (not a registered PEN: flag it if asked about conformance).
- Single-frame templates: source files with TemporalPositionIdentifier 1 and Philips phase (2001,1008) 1, one per
  unique ImagePositionPatient z, sorted by z; NumberOfTemporalPositions and (2001,1081) forced to 1 (PACS ordering).

## Tag length limits (check the worst case before adding text)
LO 64 (SeriesDescription, ProtocolName, ContentDescription, Manufacturer, SoftwareVersions, InstitutionName),
SH 16 (StudyID, LUTLabel), CS 16 (ContentLabel, uppercase), ST 1024 (DerivationDescription), LT 10240 (ImageComments).
Longest SeriesDescription today: "ASL CBF preACZ COLOR (ClinicalASL, NOT CE MARKED)" = 50 chars.

## Verification
Default: quick compile/import check now, then `tools\nightly_verify.ps1 -Name <name> -Steps export,dciodvfy` tonight
(export-only with the image's own code vs the working tree on two workingdir copies, `compare_dicom.py` + maps).
Expected for a pure refactor: all files identical apart from UIDs/times. Manual steps, if the user asks for it now:

## Fast test loop (no full run, ~10 min on G:)
1. Copy a COMPLETE workingdir (needs `DICOMORIG/` and `NIFTI/`, not only `ASL/`) into `tests\<name>\work`;
   never write into the original. Only `G:\DATA\MOYAMOYA\IMAGER\workingdir` is complete: the archived
   `tests\v*_image\*\work` folders keep just `ASL/`, `FIGURE_RESULTS/` and the dcm2niix logs (CVR mode only;
   for single mode use `--mode single` on the same copy).
2. Run the working-tree code in the released image (PowerShell):
   `docker run --rm --network none -w /fix -e PYTHONPATH=/fix -v "<repo>\python:/fix:ro" -v "<copy>:/workingdir" -v "<out>:/output" --entrypoint python jeroensnl/mri_diamox_umcu_clinicalasl_cvr_imager:1.2.0 /fix/export_results.py /workingdir /output [--mode single]`
3. Check: exit code, `Select-String <out>\clinicalasl.log -Pattern 'ERROR|Traceback'`, file/series counts against the
   table above, `python tools\compare_dicom.py <reference out> <out>`; NIfTIs bit-identical to the reference.
4. Validate: `python tools\run_dciodvfy.py <out> G:\DATA\MOYAMOYA\IMAGER\dciodvfy\<name>` (host python; dicom3tools in
   `I:\ProtonDrive\My files\Documents\studies - IMAGR\dicom3tools_clunie`). Expected remaining error:
   PatientOrientation on the colour series only.
5. Optional visual check: open the PNGs in work/FIGURE_RESULTS (montage, no patient text).

## Known open points
- Multiframe (scanner Enhanced MR) path is not covered by the ASLtest data. The multiframe colour export uses the
  Multi-frame True Color SC SOP class with PALETTE COLOR 16-bit data and keeps Enhanced MR functional groups: likely
  non-conformant; verify with dciodvfy before relying on it.
- An all-zero/NaN map is stored as zeros (scaling factor 1, warning in the log).
