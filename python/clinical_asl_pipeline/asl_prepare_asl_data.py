"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

ASL data preparation module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Functions for preparing ASL MRI data, including interleaving control/label images,
    Look-Locker correction, and saving processed NIfTI files.

License: BSD 3-Clause License
"""

import os
import logging
import numpy as np
import nibabel as nib
from clinical_asl_pipeline.asl_interleave_control_label import asl_interleave_control_label
from clinical_asl_pipeline.utils.save_data_nifti import save_data_nifti

def asl_prepare_asl_data(subject, filename, phase_tag):
    # Prepare (multidelay ASL data for analysis by interleaving control and label files per PLD, M0, and performing Look-Locker Correction.
    # Parameters:
    # subject: dict containing subject information including paths and parameters
    # filename: string, name of the NIfTI file containing ASL data
    # phase_tag: string, e.g. preAZ or postACZ phase_tag for the keys in the subject dictionary to store results
    # Returns:
    # subject: updated subject dictionary with processed ASL data     
   
    # Use a shorter alias for subject[phase_tag]
    phase_data = subject[phase_tag]

    # Load data
    data_path = os.path.join(subject['NIFTIdir'], filename)

    # remove scaling, as nibabel nib.load consumes the slope and intercept automatically
    data = nib.load(data_path).get_fdata()/nib.load(data_path).dataobj.slope 
    
    # Log NIFTI template path  
    logging.info(f"Template NIFTI path: {phase_data['templateNII_path']}")
    
    dims = data.shape

    phase_data['M0ASL_allPLD'] = np.zeros((dims[0], dims[1], dims[2], phase_data['NDYNS'], phase_data['NPLDS'], 2))
    phase_data['ASL_label1label2_allPLD'] = np.zeros((dims[0], dims[1], dims[2], phase_data['NREPEATS'] * 2, phase_data['NPLDS']))

    # Split label/control, apply Look Locker correction
    for i in range(phase_data['NPLDS']):
        idx_label = slice(i, phase_data['NPLDS'] * phase_data['NDYNS'] * 2, 2 * phase_data['NPLDS'])
        idx_control = slice(i + phase_data['NPLDS'], phase_data['NPLDS'] * phase_data['NDYNS'] * 2, 2 * phase_data['NPLDS'])

        phase_data['M0ASL_allPLD'][:, :, :, :phase_data['NDYNS'], i, 0] = data[:, :, :, idx_label] / phase_data['LookLocker_correction_factor_perPLD'][i]
        phase_data['M0ASL_allPLD'][:, :, :, :phase_data['NDYNS'], i, 1] = data[:, :, :, idx_control] / phase_data['LookLocker_correction_factor_perPLD'][i]

    # Interleave control/label, save per-PLD
    for i in range(phase_data['NPLDS']):
        interleaved = asl_interleave_control_label(np.squeeze(phase_data['M0ASL_allPLD'][:, :, :, 1:, i, 0]), np.squeeze(phase_data['M0ASL_allPLD'][:, :, :, 1:, i, 1]))
        phase_data['ASL_label1label2_allPLD'][:, :, :, :, i] = interleaved

    logging.info("Saving ASL data interleaved label control: all PLDs for AAT")
    # Swap axes to interleave PLDs and time correctly
    reordered = np.transpose(phase_data['ASL_label1label2_allPLD'], (0, 1, 2, 4, 3))  # (x, y, z, PLD, time)    
    # Now reshape so time dimension becomes interleaved PLDs
    PLDall = reordered.reshape(dims[0], dims[1], dims[2], phase_data['NPLDS'] * phase_data['NREPEATS'] * 2)
    save_data_nifti(PLDall, phase_data['PLDall_labelcontrol_path'], phase_data['templateNII_path'], 1, None, phase_data['TR'])

    logging.info("Saving ASL data interleaved label control: 2-to-last PLDs for CBF")
    # Extract PLDs 2 to end → shape: (x, y, z, time, NPLDS-1)
    PLD2tolast = phase_data['ASL_label1label2_allPLD'][:, :, :, :, 1:]
    # Reorder to (x, y, z, PLD, time), then reshape
    PLD2tolast = np.transpose(PLD2tolast, (0, 1, 2, 4, 3)).reshape(dims[0], dims[1], dims[2], (phase_data['NPLDS'] - 1) * phase_data['NREPEATS'] * 2)
    save_data_nifti(PLD2tolast, phase_data['PLD2tolast_labelcontrol_path'], phase_data['templateNII_path'], 1, None, phase_data['TR'])
    
    logging.info("Saving ASL data interleaved label control: 1-to-2 PLDs for ATA")
    # Extract first 2 PLDs → shape: (x, y, z, time, 2)
    PLD1to2 = phase_data['ASL_label1label2_allPLD'][:, :, :, :, 0:2]    
    # Reorder to (x, y, z, PLD, time), then reshape
    PLD1to2 = np.transpose(PLD1to2, (0, 1, 2, 4, 3)).reshape(dims[0], dims[1], dims[2], 2 * phase_data['NREPEATS'] * 2)
    save_data_nifti(PLD1to2, phase_data['PLD1to2_labelcontrol_path'], phase_data['templateNII_path'], 1, None, phase_data['TR'])
    
    # M0 image construction
    phase_data['M0_allPLD'] = np.mean(phase_data['M0ASL_allPLD'][:, :, :, 0, :, :], axis=4)
    phase_data['M0'] = phase_data['M0_allPLD'][:, :, :, 0]

    logging.info("Saving M0 image")
    save_data_nifti(phase_data['M0'], phase_data['M0_path'] , phase_data['templateNII_path'], 1, None, phase_data['TR'])

    return subject
