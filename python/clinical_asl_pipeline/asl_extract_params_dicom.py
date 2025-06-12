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
    #   Extract parameters from DICOM metadata for ASL data.
    #   This function reads DICOM metadata, extracts relevant parameters such as TR, TE, PLDs, slice timing, and voxel size,
    #   and saves them into the provided subject dictionary.
    #   Parameters:
    #   subject: dict containing at least 'DICOMsubjectdir', 'tau', 'labeleff', and 'N_BS'
    #   filename: string, name of the DICOM file    
    #   Returns: updated subject dictionary with extracted parameters and DICOM info    
    #

    start_time = time.time()
    logging.info(f"  Reading DICOM info... : {filename} for context: {context_tag}")

    # ? Safely join paths to avoid missing slashes
    dicom_path = os.path.join(subject['DICOMsubjectdir'], filename)
    if not os.path.exists(dicom_path):
        raise FileNotFoundError(f"DICOM file not found: {dicom_path}")

    # Load DICOM metadata
    info = pydicom.dcmread(dicom_path, stop_before_pixels=True)
    elapsed_time = round(time.time() - start_time)
    logging.info(f'..DICOM read completed in: {elapsed_time} s')   

    # Access the tags (nested)
    TR = info.CardiacRRIntervalSpecified
    logging.info(f"TR: {TR/1e3} s")

    nplds = info.get((0x2001, 0x1017)).value
    logging.info(f"NPLDS (Number of PLDs): {nplds}")  # Number of PLDS 

    private_seq = info.PerFrameFunctionalGroupsSequence[0].get((0x2005, 0x140f))

    echo_time = float(private_seq[0].get((0x0018, 0x0081)).value)  # EchoTime tag
    logging.info(f"EchoTime: {echo_time} ms")

    ndyns = int(private_seq[0].get((0x0020, 0x0105)).value)  # Number of Temporal Positions
    logging.info(f"NDYNS (Number of Dynamics): {ndyns}")

    # Extract values
    voxel_spacing_xy = private_seq[0].PixelSpacing  # MultiValue of two elements
    slice_thickness = private_seq[0].SliceThickness  # single float

    # Combine and convert to NumPy array of floats
    voxelsize = np.array(list(voxel_spacing_xy) + [slice_thickness], dtype=np.float32)
    logging.info(f"Voxel Size (XYZ in mm): {voxelsize}")

    nslices = info.get((0x2001, 0x1018)).value
    logging.info(f"NSLICES (Number of Slices): {nslices}")

    flipangle = float(private_seq[0].get((0x0018, 0x1314)).value)
    logging.info(f"FlipAngle: {flipangle} degrees")

    # Save values to subject
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

    # --- 1.  PLD and slice timing extraction ---      
    # Read all frametimes from DICOM
    frametimes = []
    
    num_frames =subject[context_tag]['NSLICES'] * subject[context_tag]['NPLDS'] * subject[context_tag]['NDYNS'] * 2
    
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
    
    # Save raw frametimes
    subject[context_tag]['frametimes'] = frametimes
    
    # Step 2: Extract PLDs (first NPLDS frames), convert to seconds
    subject[context_tag]['PLDS'] = np.array([
        ft / 1e3 for ft in frametimes[:subject[context_tag]['NPLDS']] if ft is not None
    ], dtype=np.float32)
    
    # Step 3: Compute TIs = PLD + tau
    subject[context_tag]['TIS'] = subject[context_tag]['PLDS'] + subject[context_tag]['tau']
    
    # Step 4: Estimate slice timing
    # Unique + sorted frametimes for first NSLICES frames
    clean_frametimes = [ft for ft in frametimes if ft is not None]
    unique_sorted_frametimes = np.unique(np.sort(clean_frametimes))
    
    if len(unique_sorted_frametimes) >= subject[context_tag]['NSLICES']:
        diffs = np.diff(unique_sorted_frametimes[:subject[context_tag]['NSLICES']])
        subject[context_tag]['slicetime'] = round(float(np.mean(diffs)), 1)  # ms
        logging.info(f"Estimated slice timing (ms): {subject[context_tag]['slicetime']}")
    else:
        subject[context_tag]['slicetime'] = None
        warnings.warn("Not enough valid frame times to estimate slice timing.")

    # Labeling efficiencies
    subject[context_tag]['alpha_inv'] = subject[context_tag]['labeleff']
    subject[context_tag]['alpha_BS'] = 0.95 ** subject[context_tag]['N_BS']
    subject[context_tag]['alpha'] = subject[context_tag]['alpha_inv'] * subject[context_tag]['alpha_BS']

    # M0 timing
    subject[context_tag]['TR_M0'] = subject[context_tag]['tau'] + subject[context_tag]['PLDS']

    # Store full DICOM info
    subject[context_tag][f"{filename}_DCMinfo"] = info

    return subject
