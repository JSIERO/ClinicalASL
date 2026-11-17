#!/usr/bin/env python3
"""
ClinicalASL - re-run only the results/export step (step 12) from an existing working directory.

Uses the QASL maps, brain masks and stimulus->baseline registered maps already in <workingdir>/ASL, and the source
DICOM/NIfTI in <workingdir>/DICOMORIG and NIFTI as templates, to regenerate the result NIfTIs, PNGs and DICOMs.
For testing PNG/DICOM export changes in about a minute instead of a full pipeline run (~1 h).
Overwrites the result NIfTIs/PNGs in <workingdir> (work on a copy to keep the original); DICOMs go to <outputdir>.

Usage: python export_results.py <workingdir> <outputdir> [--mode cvr|single] [--inference-method ssvb|vaby] [--config cfg.json]
"""
import argparse
import glob
import logging
import os
import pydicom
from clinical_asl_pipeline import main_pipeline as mp
from clinical_asl_pipeline.asl_extract_params_dicom import asl_extract_params_dicom
from clinical_asl_pipeline.asl_save_results_cbfaatcvr import asl_save_results_cbfaatcvr
from clinical_asl_pipeline.utils.load_parameters import load_parameters
from clinical_asl_pipeline.utils.setup_logging import setup_logging
from clinical_asl_pipeline.__version__ import __version__

MULTIFRAME_MR_UID = "1.2.840.10008.5.1.4.1.1.4.1"

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("workingdir", help="working directory of an earlier full run (ASL, NIFTI, DICOMORIG)")
    p.add_argument("outputdir", help="output directory for the DICOMs and log")
    p.add_argument("--mode", choices=list(mp.ANALYSIS_MODES), default=None, help="as in run_pipeline.py (default: config, 'cvr')")
    p.add_argument("--inference-method", choices=["ssvb", "vaby"], default=None)
    p.add_argument("--config", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                     "clinical_asl_pipeline", "config", "config_default.json"))
    a = p.parse_args()

    params = load_parameters(config_path=a.config)
    if a.inference_method:
        params["inference_method"] = a.inference_method
    params = mp.apply_analysis_mode(params, a.mode)
    os.makedirs(a.outputdir, exist_ok=True)
    setup_logging(a.outputdir)
    logging.info(f"ClinicalASL v{__version__}: export only (step 12) from {a.workingdir}")

    subject = mp.prepare_subject_paths(dict(params), None, a.outputdir, a.workingdir)
    subject = mp.prepare_input_output_paths(subject)
    orig = sorted(glob.glob(os.path.join(subject["DICOMorigdir"], "*")))[:20]
    subject["is_multiframe"] = any(
        pydicom.dcmread(f, stop_before_pixels=True, force=True).get("SOPClassUID", "") == MULTIFRAME_MR_UID for f in orig)
    subject["is_singleframe"] = not subject["is_multiframe"]
    for context, study_tag in zip(subject["ASL_CONTEXT"], subject["context_study_tags"]):
        subject = mp.get_latest_source_data(subject, study_tag, context_tag=context)
        subject = asl_extract_params_dicom(subject, context_tag=context)

    asl_save_results_cbfaatcvr(subject)
    logging.info("Export finished.")

if __name__ == "__main__":
    main()
