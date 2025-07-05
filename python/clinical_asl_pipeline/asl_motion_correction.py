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
    # perform motion correction on PLD ordered data (makes sense for CBF/AAT fit)
    # method: ANTs with DenseRigid, 6DOF, MattesMutualInformation as cost function metric
    #
    # Parameters:
    #   subject: dict containing subject information including paths and parameters
    #   context_tag: string, e.g. 'baseline' or 'stimulus' context_tag for the keys in the subject dictionary to store results
    # Returns:
    #   motion-corrected PLD ordered NIFTIs, filename appended with '_mc' before '.nii.gz'

    # Use a shorter alias for subject[context_tag]
    context_data = subject[context_tag]

    inputdata_path = context_data['PLDall_controllabel_path']
    refdata_path = context_data['M0_path']
    outputdata_path = append_mc(inputdata_path)
    nifti_template_path = context_data['sourceNIFTI_path']
    NREPEATS = context_data['NREPEATS']
    
    # update path to motion corrected data, appeding '_mc' to filename using append_mc
    context_data['PLDall_controllabel_path'] =  append_mc(context_data['PLDall_controllabel_path'])
    context_data['PLD2tolast_controllabel_path'] =  append_mc(context_data['PLD2tolast_controllabel_path'])
    context_data['PLD1to2_controllabel_path'] =  append_mc(context_data['PLD1to2_controllabel_path'])

    # perform motion correction routine, append output file name with '_mc' prefix
    asl_motioncorrection_ants(inputdata_path, refdata_path, outputdata_path)    

    PLDall_motioncorrected = nib.load(outputdata_path).get_fdata()
    PLD2tolast = PLDall_motioncorrected[:, :, :, NREPEATS*2: ]
    PLD1to2 = PLDall_motioncorrected[:, :, :, 0:NREPEATS*2*2]        

    logging.info("Saving ASL motion-corrected data interleaved label control: all PLDs for AAT")
    logging.info("Saving ASL motion-corrected data interleaved label control: 2-to-last PLDs for CBF")
    logging.info("Saving ASL motion-corrected data interleaved label control: 1-to-2 PLDs for ATA")

    save_data_nifti(PLD2tolast, context_data['PLD2tolast_controllabel_path'], nifti_template_path, 1, None, None)
    save_data_nifti(PLD1to2, context_data['PLD1to2_controllabel_path'], nifti_template_path, 1, None, None)
    
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
