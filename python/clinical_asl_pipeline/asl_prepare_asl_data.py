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
from clinical_asl_pipeline.utils.append_mc import append_mc
from clinical_asl_pipeline.asl_motioncorrection_ants import asl_motioncorrection_ants

def asl_prepare_asl_data(subject, filename, context_tag, motion_correction=True):
    # Prepare (multidelay ASL data for analysis by interleaving control and label files per PLD, M0, and performing Look-Locker Correction.
    # Parameters:
    # subject: dict containing subject information including paths and parameters
    # filename: string, name of the NIfTI file containing ASL data
    # motioncorrection: True/Fals, perform motion correction on PLD ordered data (makes sense for cbf/aat fit) using ANTs with DenseRigid
    # context_tag: string, e.g. 'baseline' or 'stimulus' context_tag for the keys in the subject dictionary to store results
    # Returns:
    # subject: updated subject dictionary with processed ASL data     

    # Use a shorter alias for subject[context_tag]
    context_data = subject[context_tag]

    # Load source multidelay ASL data: dims [x, y, z, timepoints = control/label 2 x NPLDS x NDYNS], in this order
    data_path = os.path.join(subject['NIFTIdir'], filename)

    # remove scaling, as nibabel nib.load consumes the slope and intercept automatically
    M0ASL_allPLD = nib.load(data_path).get_fdata()/nib.load(data_path).dataobj.slope 
    dims = M0ASL_allPLD.shape

    # multidelay LookLocker ASL data, M0ASL_allPLD
    NREPEATS = context_data['NREPEATS'] # excluding the M0 (first volume)
    NDYNS = context_data['NDYNS']
    NPLDS = context_data['NPLDS']
    
    # M0ASL_allPLD is 6D numpy array (x, y, z,  NDYNS, NPLDS, control/label) -> to compute deltaM for outlier removal
    context_data['M0ASL_allPLD'] = np.zeros((dims[0], dims[1], dims[2], NDYNS, NPLDS, 2))
    # ASL_controllabel_allPLD is 5D numpy array (x, y, z,  NREPEATS x 2, NPLDS) with interleaved control label volumes -> fed to QASL analysis
    context_data['ASL_controllabel_allPLD'] = np.zeros((dims[0], dims[1], dims[2], NREPEATS * 2, NPLDS))

    # Split control/label, store in array, apply Look Locker correction
    for i in range(NPLDS):
        # slice object to index array, to extract control and label volumes sorted per PLD and DYNAMIC
        idx_control = slice(i, NPLDS * NDYNS * 2, 2 * NPLDS)
        idx_label = slice(i + NPLDS, NPLDS * NDYNS * 2, 2 * NPLDS)
        #  extract control and label volumes sorted and store in M0ASL_allPLD [x, y , z, NDYNS, NPLDS, control/label], apply Look Locker correction for each PLD
        context_data['M0ASL_allPLD'][:, :, :, :NDYNS, i, 0] = M0ASL_allPLD[:, :, :, idx_control] / context_data['LookLocker_correction_factor_perPLD'][i]
        context_data['M0ASL_allPLD'][:, :, :, :NDYNS, i, 1] = M0ASL_allPLD[:, :, :, idx_label] / context_data['LookLocker_correction_factor_perPLD'][i]

    # M0 image construction
    context_data['M0_allPLD'] = np.mean(context_data['M0ASL_allPLD'][:, :, :, 0, :, :], axis=4)
    context_data['M0'] = context_data['M0_allPLD'][:, :, :, 0] # take M0 from first PLD as calibration M0 for quantification

    # Log NIFTI template path  
    logging.info(f"Template NIFTI path: {context_data['templateNII_path']}")
    logging.info("Saving M0 image")

    save_data_nifti(context_data['M0'], context_data['M0_path'] , context_data['templateNII_path'], 1, None, context_data['TR'])

    # Interleave control/label, save per-PLD
    # resulting ASL_controllabel_allPLD is 5D numpy array (x, y, z,  NREPEATS x 2, NPLDS) with interleaved control label volumes -> QASL analysis
    for i in range(NPLDS):
        interleaved = asl_interleave_control_label(np.squeeze(context_data['M0ASL_allPLD'][:, :, :, 1:, i, 0]), np.squeeze(context_data['M0ASL_allPLD'][:, :, :, 1:, i, 1]))
        context_data['ASL_controllabel_allPLD'][:, :, :, :, i] = interleaved # 4D numpy array (x, y, z,  NREPEATS *2 ) with interleaved control label volumes

    # Swap axes to interleave PLDs and time correctly
    reordered = np.transpose(context_data['ASL_controllabel_allPLD'], (0, 1, 2, 4, 3))  # (x, y, z, PLD, time)
    # Now reshape so time dimension becomes interleaved PLDs
    PLDall = reordered.reshape(dims[0], dims[1], dims[2], NPLDS * NREPEATS * 2)
    save_data_nifti(PLDall, context_data['PLDall_controllabel_path'], context_data['templateNII_path'], 1, None, context_data['TR'])


    if motion_correction==True:
            # perform motion correction routine
            asl_motioncorrection_ants(context_data['PLDall_controllabel_path'], context_data['M0_path'], append_mc(context_data['PLDall_controllabel_path']))
            
            logging.info("Saving ASL motion-corrected data interleaved label control: all PLDs for AAT")
            PLDall_motioncorrected = nib.load(append_mc(context_data['PLDall_controllabel_path'])).get_fdata()
            context_data['PLDall_controllabel_path'] = append_mc(context_data['PLDall_controllabel_path']) # update path to motion corrected data

            logging.info("Saving ASL motion-corrected data interleaved label control: 2-to-last PLDs for CBF")
            PLD2tolast = PLDall_motioncorrected[:, :, :, NREPEATS*2: ]
            save_data_nifti(PLD2tolast, append_mc(context_data['PLD2tolast_controllabel_path']), context_data['templateNII_path'], 1, None, context_data['TR'])
            context_data['PLD2tolast_controllabel_path'] =  append_mc(context_data['PLD2tolast_controllabel_path']) # update path to motion corrected data

            logging.info("Saving ASL motion-corrected data interleaved label control: 1-to-2 PLDs for ATA")
            PLD1to2 = PLDall_motioncorrected[:, :, :, 0:NREPEATS*2*2]
            save_data_nifti(PLD1to2, append_mc(context_data['PLD1to2_controllabel_path']), context_data['templateNII_path'], 1, None, context_data['TR'])
            context_data['PLD1to2_controllabel_path'] =  append_mc(context_data['PLD1to2_controllabel_path']) # update path to motion corrected data
    elif motion_correction==False:
            logging.info("Saving ASL data interleaved label control: all PLDs for AAT")
            PLDall_motioncorrected = nib.load(append_mc(context_data['PLDall_controllabel_path'])).get_fdata()
            context_data['PLDall_controllabel_path'] = append_mc(context_data['PLDall_controllabel_path']) # update path to motion corrected data

            logging.info("Saving ASL data interleaved label control: 2-to-last PLDs for CBF")
            PLD2tolast = PLDall_motioncorrected[:, :, :, NREPEATS*2: ]
            save_data_nifti(PLD2tolast, append_mc(context_data['PLD2tolast_controllabel_path']), context_data['templateNII_path'], 1, None, context_data['TR'])
            context_data['PLD2tolast_controllabel_path'] =  append_mc(context_data['PLD2tolast_controllabel_path']) # update path to motion corrected data

            logging.info("Saving ASL data interleaved label control: 1-to-2 PLDs for ATA")
            PLD1to2 = PLDall_motioncorrected[:, :, :, 0:NREPEATS*2*2]
            save_data_nifti(PLD1to2, append_mc(context_data['PLD1to2_controllabel_path']), context_data['templateNII_path'], 1, None, context_data['TR'])
            context_data['PLD1to2_controllabel_path'] =  append_mc(context_data['PLD1to2_controllabel_path']) # update path to motion corrected data 
            
    return subject
