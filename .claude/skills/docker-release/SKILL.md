---
name: docker-release
description: Build, check and release the ClinicalASL Docker image for IMAGR (jeroensnl/mri_diamox_umcu_clinicalasl_cvr_imager). Use for Dockerfile/dependency changes, version bumps and the release checklist (version sync, build, nightly validation, changelog, tag, push).
---

Read the Dockerfile and AGENTS.md "Build / run" first. Use PowerShell for docker with Windows paths.

## Dockerfile structure (lean multi-stage, like FetalSVR)
- Stage `dcm2niix`: builds dcm2niix v1.0.20220720 with OpenJPEG + CharLS.
- Stage `pyenv`: micromamba env `/opt/conda/envs/qi`; CPU-only torch FIRST (so no nvidia-*/triton CUDA wheels, ~3.4 GB);
  QASL packages from bitbucket with the key as BuildKit secret; HD-BET; requirements.txt; strips test/dev packages,
  caches, headers.
- Final: ubuntu:22.04 + runtime libs (libgomp1, libopenjp2-7, libcharls2) + env + HD-BET params + `python/`.
- `docker-entrypoint.sh`: `bash|python|qasl|hd-bet|dcm2niix|...` pass through, otherwise `run_pipeline.py "$@"`.
  Kept LF via `.gitattributes` (core.autocrlf=true on this machine); the build also strips CR.
- Debug shell: `docker run --rm -it ... <image> bash`. On Linux hosts the output folder must be writable by the
  container user.

## Hard rules
- The QASL SSH key (`G:\DATA\MOYAMOYA\IMAGER\qi-key-jeroen`) is only ever passed as a BuildKit secret
  (`--secret id=qasl_key,src=...`). Never COPY it, echo it, put it in a log, the repo or chat.
- The old `:latest` image (26 GB) has the key in a layer: never push it or retag it.
- Building, pushing and tagging are outward-facing or long: ask before each (build time is long; >30 min rule).
  Never `docker push` without an explicit yes for that tag.
- Image tag = app version. The build overwrites `__version__.py` with `APP_VERSION`.

## Image facts (1.2.0)
6.6 GB, Ubuntu 22.04, micromamba env `qi` (python 3.11) copied from the build stage; CPU-only torch 2.13.0 installed
first (no nvidia-*/triton); PyPI pinned by `docker/qi-constraints.txt` (PIP_CONSTRAINT); QASL packages `--no-deps`
at the commits in `docker/qasl-commits.txt` (QASL 1.4.2, asserted); HD-BET v2.0.1 with params in `/opt/hd-bet_params`
(patched `paths.py`); toblerone numba module pre-compiled (no compiler at runtime); dcm2niix v1.0.20220720 built with
OpenJPEG + CharLS. The FINAL CHECK RUN fails the build on missing imports, CUDA torch, qasl/hd-bet CLI, dcm2niix libs.
Entrypoint: `bash|python|qasl|hd-bet|dcm2niix|...` pass through, else `run_pipeline.py "$@"`; LF endings enforced.

## Release checklist
1. Version sync (all equal): `python/clinical_asl_pipeline/__version__.py`, `config/config_default.json` "version",
   README "Version:" line + config example, `utils/load_parameters.py` fallback, Dockerfile `ARG APP_VERSION` default
   + header build comment, AGENTS.md and INSTALL.md build/run commands. `Grep` for the old version string to be sure.
2. CHANGELOG.md: new section at the top, short bullets, user-visible changes only.
3. Build (ask first):
   `docker build --secret id=qasl_key,src=G:\DATA\MOYAMOYA\IMAGER\qi-key-jeroen --build-arg APP_VERSION=<v> -t jeroensnl/mri_diamox_umcu_clinicalasl_cvr_imager:<v> .`
   Run in the background with output to a log; check for "FINAL OK" and "QASL + HD-BET OK"; `docker images` size ~6.6 GB.
4. Smoke: `docker run --rm <image:v> python run_pipeline.py --version`.
5. Validation (~1.5 h): schedule it tonight with `tools\nightly_verify.ps1 -Image <image:v> -Name v<v>_image`
   (Windows Task Scheduler, see the script header): CVR + single full runs, export-only, dciodvfy into
   `IMAGER\dciodvfy\v<v>_image`, comparison with the previous release. Read `summary.txt` the next session.
6. Record the validation numbers in CHANGELOG (duration, DICOM count/series, map differences, dciodvfy result).
7. Commit (author JSIERO, only when asked), then offer git tag `v<v>` and docker push as separate questions.

## Dependency changes
- Change a PyPI pin in `docker/qi-constraints.txt` (it constrains every pip install) and/or `requirements.txt`.
- Never let antspyx move off 0.5.4 (0.6.x registration is worse). Never add a CUDA torch.
- QASL upgrade = new commits in `docker/qasl-commits.txt` + `QASL_VERSION` arg; qasl_setup HEAD (1.6.x) needs the
  ANTs CLI, which the image does not have. A QASL change needs a full validation run and a map comparison.
