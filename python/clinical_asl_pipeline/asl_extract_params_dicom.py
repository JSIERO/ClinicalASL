"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

DICOM parameter extraction module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Function for extracting ASL MRI parameters from DICOM Philips multiframe metadata.

License: BSD 3-Clause License
"""

import os
import logging
import warnings
import pydicom
import numpy as np
from pydicom.misc import is_dicom

def asl_extract_params_dicom(subject, context_tag):
    context_data = subject[context_tag]
    dicom_dir = subject['DICOMsubjectdir']
    dicom_filename = os.path.basename(context_data['sourceDCM_path'])
    dicom_path = os.path.join(dicom_dir, dicom_filename)

    logging.info(f"Reading DICOM info... : {dicom_path} for context: {context_tag}")   

    if not os.path.exists(dicom_path):
        raise FileNotFoundError(f"DICOM file not found: {dicom_path}")

    # Load the first DICOM file
    info = pydicom.dcmread(dicom_path, stop_before_pixels=True)

    # Save the original header
    context_data[f"{dicom_filename}_DCMinfo"] = info

    # ========== MULTIFRAME FORMAT ==========
    if hasattr(info, 'PerFrameFunctionalGroupsSequence'):
        logging.info("Detected multiframe Philips DICOM.")
        logging.info(f"Extracting scan parameters from multiframe DICOM series with basename: {dicom_filename}")

        private_seq = info.PerFrameFunctionalGroupsSequence[0].get((0x2005, 0x140f))
        echo_time = float(private_seq[0].get((0x0018, 0x0081)).value) 
        flipangle = float(private_seq[0].get((0x0018, 0x1314)).value)
        age_patient = getattr(info, 'PatientAge', None)
        if age_patient:
            age_number = int(age_patient[:3])  # Takes the first 3 characters and converts to int
        else:
            age_number = None

        nplds = info.get((0x2001, 0x1017)).value
        ndyns = int(private_seq[0].get((0x0020, 0x0105)).value)
        nslices = info.get((0x2001, 0x1018)).value

        voxel_spacing_xy = private_seq[0].PixelSpacing
        slice_thickness = private_seq[0].SliceThickness
        voxelsize = np.array(list(voxel_spacing_xy) + [slice_thickness], dtype=np.float32)

        # Extract frametimes
        frametimes = []
        num_frames = nslices * nplds * ndyns
        for t in range(num_frames):
            try:
                frame = info.PerFrameFunctionalGroupsSequence[t]
                private_elem = frame.get((0x2005, 0x140f))
                item = private_elem.value[0] if private_elem else None
                trigger_elem = item.get((0x0018, 0x1060)) if item else None
                frametimes.append(float(trigger_elem.value) if trigger_elem else None)
            except Exception:
                frametimes.append(None)
                
        plds = np.array([ft / 1e3 for ft in frametimes[:nplds] if ft is not None], dtype=np.float32)       

    # ========== SINGLEFRAME FORMAT ==========
    else:
        logging.info("Detected singleframe PACS-exported DICOM series.")

        # Load all DICOMs (filter for files only, ie no folders) in same directory with matching prefix
        all_files = [f for f in os.listdir(dicom_dir) if os.path.isfile(os.path.join(dicom_dir, f)) and is_dicom(os.path.join(dicom_dir, f))]       
        series_prefix = "_".join(dicom_filename.split("_")[:-1])
        matching_files = [f for f in all_files if f.startswith(series_prefix)]
        matching_files = sorted(matching_files, key=lambda x: int(x.split("_")[-1]))        
        logging.info(f"Extracting scan parameters from singleframe DICOM series with basename: {series_prefix}")

        if not matching_files:
            raise FileNotFoundError(f"No matching DICOMs found for prefix {series_prefix} in {dicom_dir}")

        # Read all headers
        headers = [pydicom.dcmread(os.path.join(dicom_dir, f), stop_before_pixels=True) for f in matching_files]

        # Use first header to extract common fields
        hdr0 = headers[0]
        echo_time = float(getattr(hdr0, 'EchoTime', 0.0))
        flipangle = float(getattr(hdr0, 'FlipAngle', 0.0))
        age_patient = getattr(hdr0, 'PatientAge', None)
        if age_patient:
            age_number = int(age_patient[:3])  # Takes the first 3 characters and converts to int
        else:
            age_number = None

        voxelsize = np.array([
            float(hdr0.PixelSpacing[0]),
            float(hdr0.PixelSpacing[1]),
            float(getattr(hdr0, 'SliceThickness', 0.0))
        ], dtype=np.float32)

        # Infer NSLICES
        nslices = info.get((0x2001, 0x1018)).value
        nplds = info.get((0x2001, 0x1017)).value
        ndyns = info.get((0x2001, 0x1081)).value
        frametimes = [float(getattr(h, 'TriggerTime', 0)) for h in headers]

        plds = np.array([ft / 1e3 for ft in frametimes[:nplds] if ft is not None], dtype=np.float32)

    # Estimate slice timing if enough unique values in frametimes
    clean_frametimes = [ft for ft in frametimes if ft is not None]
    unique_sorted_frametimes = np.unique(np.sort(clean_frametimes))
    if len(unique_sorted_frametimes) >= nslices:
        diffs = np.diff(unique_sorted_frametimes[:nslices])
        slicetime = round(float(np.mean(diffs)), 1)
    else:
        slicetime = None
        warnings.warn("Not enough valid frame times to estimate slice timing.")

    # ==== Common fields for both formats ====
    context_data['tau'] = subject['tau']
    context_data['N_BS'] = subject['N_BS']
    context_data['labeleff'] = subject['labeleff']
    context_data['TE'] = echo_time
    context_data['slicetime'] = slicetime
    context_data['PLDS'] = plds
    context_data['NPLDS'] = nplds
    context_data['NDYNS'] = ndyns
    context_data['NREPEATS'] = ndyns - 1
    context_data['VOXELSIZE'] = voxelsize
    context_data['NSLICES'] = nslices
    context_data['FLIPANGLE'] = flipangle
    context_data['TIS'] = context_data['PLDS'] + context_data['tau']
    context_data['TR_M0'] = context_data['TIS'][0]
    context_data['age'] = age_number  

    # Labeling efficiencies
    context_data['alpha_inv'] = context_data['labeleff']
    context_data['alpha_BS'] = 0.95 ** context_data['N_BS']
    context_data['alpha'] = round(context_data['alpha_inv'] * context_data['alpha_BS'], 2)

    logging.info(f"ECHOTIME: {round(echo_time,2)} ms")
    logging.info(f"SLICETIME: {slicetime} ms")
    logging.info(f"FLIPANGLE: {flipangle} degrees")
    logging.info(f"VOXEL SIZE (XYZ): {voxelsize} mm")
    logging.info(f"NSLICES (number of slices): {nslices}")
    logging.info(f"label duration: {context_data['tau']} s")
    logging.info(f"PLDs: {plds} s")
    logging.info(f"ALPHA (effective - incl. background suppression): {context_data['alpha']}")
    logging.info(f"NPLDS (number of post-label delays): {nplds}")
    logging.info(f"NDYNS (number of dynamics): {ndyns}")
    logging.info(f"NREPEATS (control/label pair repeats): {ndyns - 1}")
    logging.info(f"Patient age: {age_patient} (numeric: {age_number})")

    return subject
