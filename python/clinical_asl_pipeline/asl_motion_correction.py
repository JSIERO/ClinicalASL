"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

ASL data preparation module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Perform motion correction on 4D datasets using a reference image.
    Method is by ANTs

License: BSD 3-Clause License
"""

import os
import ants 
import logging
import nibabel as nib
from clinical_asl_pipeline.utils.append_filename import append_mc
from clinical_asl_pipeline.utils.save_data_nifti import save_data_nifti


def asl_motion_correction(subject, context_tag):
    # perform motion correction on PLD ordered data (makes sense for cbf/aat fit) using ANTs with DenseRigid

    context_data = subject[context_tag]
    NREPEATS = context_data['NREPEATS']

    # perform motion correction routine, appenc file name with '_mc' prefix
    asl_motioncorrection_ants(context_data['PLDall_controllabel_path'], context_data['M0_path'], append_mc(context_data['PLDall_controllabel_path']))
    
    PLDall_motioncorrected = nib.load(append_mc(context_data['PLDall_controllabel_path'])).get_fdata()
    PLD2tolast = PLDall_motioncorrected[:, :, :, NREPEATS*2: ]
    PLD1to2 = PLDall_motioncorrected[:, :, :, 0:NREPEATS*2*2]
        
    logging.info("Saving ASL motion-corrected data interleaved label control: all PLDs for AAT")
    context_data['PLDall_controllabel_path'] = append_mc(context_data['PLDall_controllabel_path']) # update path to motion corrected data

    logging.info("Saving ASL motion-corrected data interleaved label control: 2-to-last PLDs for CBF")
    save_data_nifti(PLD2tolast, append_mc(context_data['PLD2tolast_controllabel_path']), context_data['templateNII_path'], 1, None, context_data['TR'])
    context_data['PLD2tolast_controllabel_path'] =  append_mc(context_data['PLD2tolast_controllabel_path']) # update path to motion corrected data

    logging.info("Saving ASL motion-corrected data interleaved label control: 1-to-2 PLDs for ATA")
    save_data_nifti(PLD1to2, append_mc(context_data['PLD1to2_controllabel_path']), context_data['templateNII_path'], 1, None, context_data['TR'])
    context_data['PLD1to2_controllabel_path'] =  append_mc(context_data['PLD1to2_controllabel_path']) # update path to motion corrected data
    
    return subject

def asl_motioncorrection_ants(inputdata, refdata, outputdata):
    # Perform motion correction using ANTs
    # inputdata: path to the input ASL data in NIfTI formatd
    # refdata: path to the reference image for motion correction
    # outputdata: path to save the motion-corrected output data in NIfTI format
    logging.info(f"Perform motion correction (using ANTs): input: {os.path.basename(inputdata)}  reference: {os.path.basename(refdata)} ")
    results_dict = ants.motion_correction(
        ants.image_read(inputdata),
        ants.image_read(refdata),
        type_of_transform="DenseRigid",
        aff_metric="mattes",
        smoothing_in_mm=True,
        singleprecision=True,
        verbose=True
    )
    ants.image_write(results_dict["motion_corrected"], outputdata)
