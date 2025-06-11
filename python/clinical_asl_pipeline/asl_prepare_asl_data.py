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

def asl_prepare_asl_data(subject, filename, context_tag):
    # Prepare (multidelay ASL data for analysis by interleaving control and label files per PLD, M0, and performing Look-Locker Correction.
    # Parameters:
    # subject: dict containing subject information including paths and parameters
    # filename: string, name of the NIfTI file containing ASL data
    # context_tag: string, e.g. 'baseline' or 'stimulus' context_tag for the keys in the subject dictionary to store results
    # Returns:
    # subject: updated subject dictionary with processed ASL data     

    # Use a shorter alias for subject[context_tag]
    context_data = subject[context_tag]

    # Load data
    data_path = os.path.join(subject['NIFTIdir'], filename)

    # remove scaling, as nibabel nib.load consumes the slope and intercept automatically
    data = nib.load(data_path).get_fdata()/nib.load(data_path).dataobj.slope 
    
    # Log NIFTI template path  
    logging.info(f"Template NIFTI path: {context_data['templateNII_path']}")
    
    dims = data.shape

    context_data['M0ASL_allPLD'] = np.zeros((dims[0], dims[1], dims[2], context_data['NDYNS'], context_data['NPLDS'], 2))
    context_data['ASL_label1label2_allPLD'] = np.zeros((dims[0], dims[1], dims[2], context_data['NREPEATS'] * 2, context_data['NPLDS']))

    # Split label/control, apply Look Locker correction
    for i in range(context_data['NPLDS']):
        idx_label = slice(i, context_data['NPLDS'] * context_data['NDYNS'] * 2, 2 * context_data['NPLDS'])
        idx_control = slice(i + context_data['NPLDS'], context_data['NPLDS'] * context_data['NDYNS'] * 2, 2 * context_data['NPLDS'])

        context_data['M0ASL_allPLD'][:, :, :, :context_data['NDYNS'], i, 0] = data[:, :, :, idx_label] / context_data['LookLocker_correction_factor_perPLD'][i]
        context_data['M0ASL_allPLD'][:, :, :, :context_data['NDYNS'], i, 1] = data[:, :, :, idx_control] / context_data['LookLocker_correction_factor_perPLD'][i]

    # Interleave control/label, save per-PLD
    for i in range(context_data['NPLDS']):
        interleaved = asl_interleave_control_label(np.squeeze(context_data['M0ASL_allPLD'][:, :, :, 1:, i, 0]), np.squeeze(context_data['M0ASL_allPLD'][:, :, :, 1:, i, 1]))
        context_data['ASL_label1label2_allPLD'][:, :, :, :, i] = interleaved

    logging.info("Saving ASL data interleaved label control: all PLDs for AAT")
    # Swap axes to interleave PLDs and time correctly
    reordered = np.transpose(context_data['ASL_label1label2_allPLD'], (0, 1, 2, 4, 3))  # (x, y, z, PLD, time)    
    # Now reshape so time dimension becomes interleaved PLDs
    PLDall = reordered.reshape(dims[0], dims[1], dims[2], context_data['NPLDS'] * context_data['NREPEATS'] * 2)
    save_data_nifti(PLDall, context_data['PLDall_labelcontrol_path'], context_data['templateNII_path'], 1, None, context_data['TR'])

    logging.info("Saving ASL data interleaved label control: 2-to-last PLDs for CBF")
    # Extract PLDs 2 to end → shape: (x, y, z, time, NPLDS-1)
    PLD2tolast = context_data['ASL_label1label2_allPLD'][:, :, :, :, 1:]
    # Reorder to (x, y, z, PLD, time), then reshape
    PLD2tolast = np.transpose(PLD2tolast, (0, 1, 2, 4, 3)).reshape(dims[0], dims[1], dims[2], (context_data['NPLDS'] - 1) * context_data['NREPEATS'] * 2)
    save_data_nifti(PLD2tolast, context_data['PLD2tolast_labelcontrol_path'], context_data['templateNII_path'], 1, None, context_data['TR'])
    
    logging.info("Saving ASL data interleaved label control: 1-to-2 PLDs for ATA")
    # Extract first 2 PLDs → shape: (x, y, z, time, 2)
    PLD1to2 = context_data['ASL_label1label2_allPLD'][:, :, :, :, 0:2]    
    # Reorder to (x, y, z, PLD, time), then reshape
    PLD1to2 = np.transpose(PLD1to2, (0, 1, 2, 4, 3)).reshape(dims[0], dims[1], dims[2], 2 * context_data['NREPEATS'] * 2)
    save_data_nifti(PLD1to2, context_data['PLD1to2_labelcontrol_path'], context_data['templateNII_path'], 1, None, context_data['TR'])
    
    # M0 image construction
    context_data['M0_allPLD'] = np.mean(context_data['M0ASL_allPLD'][:, :, :, 0, :, :], axis=4)
    context_data['M0'] = context_data['M0_allPLD'][:, :, :, 0]

    logging.info("Saving M0 image")
    save_data_nifti(context_data['M0'], context_data['M0_path'] , context_data['templateNII_path'], 1, None, context_data['TR'])

    return subject
