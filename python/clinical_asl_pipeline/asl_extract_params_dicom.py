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
import time
import numpy as np

def asl_extract_params_dicom(subject, filename, context_tag):
    #
    # Extracts ASL (Arterial Spin Labeling) parameters from a DICOM file and updates the subject dictionary.
    # This function reads DICOM metadata for ASL data, extracts relevant imaging parameters such as TR, TE, 
    # number of post-label delays (PLDs), slice timing, voxel size, and other acquisition details, and saves 
    # them into the provided subject dictionary under the specified context tag.
    # Parameters
    # ----------
    #subject : dict
    #    Dictionary containing at least the keys 'DICOMsubjectdir', 'tau', 'labeleff', and 'N_BS' - these are NOT read from the DICOM (action for Philips)
    #    The function will update this dictionary with extracted DICOM parameters.
    # filename : str
    #    Name of the DICOM file to read (relative to 'DICOMsubjectdir').
    # context_tag : str
    #     Key under which to store the extracted parameters in the subject dictionary.
    # Returns
    # -------
    # subject : dict
    #     The updated subject dictionary with extracted DICOM parameters and metadata.
    # Notes
    # -----
    # - The function expects certain private DICOM tags and sequences, which may be specific to the scanner/vendor.
    # - Some parameters (e.g., 'tau', 'labeleff', 'N_BS') are not stored in DICOM and must be provided by the user.
    # - The function also estimates slice timing and computes effective labeling efficiency including background suppression.
    #- Full DICOM metadata is stored in the subject dictionary for reference.

    start_time = time.time()
    logging.info(f"  Reading DICOM info... : {filename} for context: {context_tag}")

    # ? Safely join paths to avoid missing slashes
    dicom_path = os.path.join(subject['DICOMsubjectdir'], filename)
    if not os.path.exists(dicom_path):
        raise FileNotFoundError(f"DICOM file not found: {dicom_path}")

    # Load DICOM metadata
    info = pydicom.dcmread(dicom_path, stop_before_pixels=True)

    elapsed_time = round(time.time() - start_time)
    logging.info(f'  DICOM read completed in: {elapsed_time} s')   

    # Access the tags (nested)
    # link to the private tag of multiframe Philips DICOM 
    private_seq = info.PerFrameFunctionalGroupsSequence[0].get((0x2005, 0x140f))
    
    TR = info.CardiacRRIntervalSpecified
    logging.info(f"TR: {TR/1e3} s")

    echo_time = float(private_seq[0].get((0x0018, 0x0081)).value) 
    logging.info(f"ECHOTIME: {round(echo_time,2)} ms")
    
    flipangle = float(private_seq[0].get((0x0018, 0x1314)).value)
    logging.info(f"FLIPANGLE: {flipangle} degrees")

    nplds = info.get((0x2001, 0x1017)).value
    logging.info(f"NPLDS (number of post-label delays): {nplds}")

    ndyns = int(private_seq[0].get((0x0020, 0x0105)).value)
    logging.info(f"NDYNS (number of dynamics): {ndyns}")

    nrepeats = ndyns - 1
    logging.info(f"NREPEATS (control/label pair repeats): {nrepeats}")

    # Extract values
    voxel_spacing_xy = private_seq[0].PixelSpacing  # MultiValue of two elements
    slice_thickness = private_seq[0].SliceThickness  # single float

    # Combine and convert to NumPy array of floats
    voxelsize = np.array(list(voxel_spacing_xy) + [slice_thickness], dtype=np.float32)
    logging.info(f"VOXEL SIZE (XYZ): {voxelsize} mm")

    nslices = info.get((0x2001, 0x1018)).value
    logging.info(f"NSLICES (number of slices): {nslices}")

    # --- 1.  PLD and slice timing extraction ---      
    # Read all frametimes from DICOM
    frametimes = []
    
    num_frames = nslices * nplds * ndyns
    
    for t in range(num_frames):
        try:
            frame = info.PerFrameFunctionalGroupsSequence[t]
            private_elem = frame.get((0x2005, 0x140f))
            private_seq = private_elem.value if private_elem else None
    
            if private_seq and len(private_seq) > 0:
                item = private_seq[0]
                trigger_elem = item.get((0x0018, 0x1060))  # TriggerTime
                if trigger_elem:
                    frametimes.append(float(trigger_elem.value))
                else:
                    frametimes.append(None)
            else:
                frametimes.append(None)
        except IndexError:
            frametimes.append(None)

    # Save values to subject context dict
    subject[context_tag]['tau'] = subject['tau'] # labelling duration in  not stored in DICOM,  Philips issue, but provided by user
    subject[context_tag]['N_BS'] = subject['N_BS'] # not stored in DICOM, but provided by user      
    subject[context_tag]['labeleff'] = subject['labeleff'] # not stored in DICOM, but provided by user
    subject[context_tag]['TR'] = TR/1e3 # in s
    subject[context_tag]['TE'] = echo_time# in ms
    subject[context_tag]['NPLDS'] = nplds
    subject[context_tag]['NDYNS'] = ndyns
    subject[context_tag]['NREPEATS'] =subject[context_tag]['NDYNS'] - 1
    subject[context_tag]['VOXELSIZE'] = voxelsize # in mm
    subject[context_tag]['NSLICES'] = nslices
    subject[context_tag]['FLIPANGLE'] = flipangle
    
    # Save raw frametimes
    subject[context_tag]['frametimes'] = frametimes
    
    # Step 2: Extract PLDs (first NPLDS frames), convert to seconds
    subject[context_tag]['PLDS'] = np.array([ft / 1e3 for ft in frametimes[:subject[context_tag]['NPLDS']] if ft is not None], dtype=np.float32)
    logging.info(f"PLDs: {subject[context_tag]['PLDS']} s")
    
    # Step 3: Compute TIs = PLD + tau
    subject[context_tag]['TIS'] = subject[context_tag]['PLDS'] + subject[context_tag]['tau']

    # TR fo the M0 at the first dynamic
    subject[context_tag]['TR_M0'] = subject[context_tag]['tau'] + subject[context_tag]['PLDS']

    # Step 4: Estimate slice timing
    # Unique + sorted frametimes for first NSLICES frames
    clean_frametimes = [ft for ft in frametimes if ft is not None]
    unique_sorted_frametimes = np.unique(np.sort(clean_frametimes))
    
    if len(unique_sorted_frametimes) >= subject[context_tag]['NSLICES']:
        diffs = np.diff(unique_sorted_frametimes[:subject[context_tag]['NSLICES']])
        subject[context_tag]['slicetime'] = round(float(np.mean(diffs)), 1)  # ms
        logging.info(f"SLICETIMING: {subject[context_tag]['slicetime']} ms")
    else:
        subject[context_tag]['slicetime'] = None
        warnings.warn("Not enough valid frame times to estimate slice timing.")

    # Labeling (inversion) efficiencies
    subject[context_tag]['alpha_inv'] = subject[context_tag]['labeleff']
    subject[context_tag]['alpha_BS'] = 0.95 ** subject[context_tag]['N_BS']
    subject[context_tag]['alpha'] = round(subject[context_tag]['alpha_inv'] * subject[context_tag]['alpha_BS'],2)
    logging.info(f"ALPHA (effective - incl. background suppresion): {subject[context_tag]['alpha']}")

    # Store full DICOM info
    subject[context_tag][f"{filename}_DCMinfo"] = info

    return subject
