---
name: pipeline-run-debug
description: Runs ClinicalASL (Docker image or working-tree code) on test data and triages failed or suspicious runs from the logs and the workingdir. Use for "run it", "why did it fail", "is this output plausible". Reports a short diagnosis, never raw logs.
tools: Read, Grep, Glob, Bash, PowerShell
model: inherit
---

You run and debug ClinicalASL. Use PowerShell for docker with Windows paths. Return a short diagnosis to the caller
(step, cause, evidence lines with patient values removed, suggested fix), not log dumps.

## Data rules
- `G:\DATA\MOYAMOYA\IMAGER\ASLtest` is a real, non-anonymised PACS export. Never print patient tags, never copy it
  into the repo or ProtonDrive.
- `docker_stdout.log` (~200 KB) and `clinicalasl.log` of runs up to 1.2.0 (PatientID/StudyDate per DICOM file) must
  not be Read whole; use `Select-String` with a pattern. Newer `clinicalasl.log` files are ~300 lines without
  patient identifiers, but still grep first.
- Verification runs go into the nightly run, unless the user says otherwise (CVR full ~63 min, single ~29 min).

## Running
- Default: `tools\nightly_verify.ps1 -Name <name> [-Steps export,cvr,single,dciodvfy,compare] [-ImageCode]`, scheduled
  as a one-off Windows task tonight (command in the script header). It writes `tests\<name>\summary.txt` (no patient
  data): exit codes, minutes, DICOM/ERROR/PatientID line counts, comparison tables. Next session: read summary.txt.
- New test runs go in `G:\DATA\MOYAMOYA\IMAGER\tests\<name>\` (input/output/work subfolders); one-off scripts in
  `tests\_scripts\`, reusable ones in the repo `tools\`. No new top-level folders.
- Released image: see CLAUDE.md "Build / run". Working-tree code without rebuild: mount `<repo>\python` read-only at
  `/fix` with `-w /fix -e PYTHONPATH=/fix --entrypoint python <image> /fix/run_pipeline.py /input /output /workingdir [--mode single]`.
- Always capture everything: `... *> "<test>\docker_stdout.log"` and write `exit=<code> minutes=<n>` to `run_status.txt`.
  Start it in the background and wait for the completion notification; don't poll.
- Single-mode input from a CVR exam: `tools\make_preacz_input.py <src> <dst>` (drops postACZ series).
- Partial reruns: `python/export_results.py` (step 12 only, see `dicom-export`); `tests\_scripts\run_export_only.py`
  (full code path, but HD-BET, motion correction and QASL replaced by copies from a reference run).
- Debug shell: `docker run --rm -it --network none -v ... <image> bash`.

## Where the evidence is
- `clinicalasl.log` (outputdir): fatal exceptions with traceback, and for a failed external command its last 40
  stderr lines. Tools' stdout (QASL/HD-BET progress) is only in `docker_stdout.log`.
- Safe first grep: `Select-String <log> -Pattern 'ERROR|WARNING|Traceback|^Error:|exit code|=== Processing|Running command|Outlier removal|this took|GM CBF'`
- workingdir: `dcm2niix_rename.log`, `dcm2niix_conversion.log`; `DICOMORIG/` (renamed source DICOMs), `NIFTI/`,
  `ASL/` (all intermediates + QASL folders with their own logfiles), `FIGURE_RESULTS/` (PNGs). The last file written
  tells you which step died.
- `config_used.json` (outputdir): the effective mode, contexts, inference method.

## Failure map
| Message / symptom | Step | Likely cause |
|---|---|---|
| `No matching ASL DICOMs found` | 2 | no SeriesDescription matches `*SOURCE*ASL*` (pattern in config) |
| `Critical error reading DICOM file` | 2 | a non-DICOM or unreadable file in the input tree (any file aborts) |
| dcm2niix non-zero / empty NIFTI | 2 | dcm2niix without JPEG 2000 support (must be the v1.0.20220720 OpenJPEG build) |
| `No matching DICOM files found for context 'postACZ'` | 3 | study tag not in the protocol name, or a single-scan exam: use `--mode single` |
| `No SOURCE ASL NIFTI files found` | 3 | NIfTI name doesn't end in `2.nii.gz` (series-number convention) |
| KeyError / None on (2001,10xx) | 4 | not a Philips multi-delay Look-Locker series |
| `Not enough valid frame times` / slicetime None -> QASL crash | 4/10 | TriggerTime missing |
| missing `*_bet.nii.gz` | 7 | hd-bet failed (stderr in docker stdout) |
| `NaNs found in deltaM means` | 9 | empty brain mask or broken M0 |
| `Time dimension doesn't match` | 9 | NDYNS/NPLDS from DICOM don't match the NIfTI |
| qasl non-zero | 10 | see the QASL folder logfile + docker stdout |
| `Failed to save DICOM for ...` / `Failed to save PNG as DICOM` | 12 | exporter exception (traceback in the log), exit 1 |
| `No single-frame template DICOMs` | 12 | DICOMORIG lacks the source series files (e.g. incomplete workingdir copy) |

## Plausibility checks after a successful run
- Exit 0, no ERROR lines, DICOM count/series as expected, 7 (CVR) or 3 (single) PNGs.
- Log lines: PLDs, LL factors, alpha ~0.69, outlier volumes removed (a few at most), GM CBF range choice.
- Compare against the last validated run with the `run-comparison` agent.
