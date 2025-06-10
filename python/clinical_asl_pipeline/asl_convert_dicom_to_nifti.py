"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

DICOM to NIfTI conversion module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Function for converting ASL MRI DICOM files to NIfTI format using dcm2niix.

License: BSD 3-Clause License
"""

import os
import logging
import shutil
import subprocess
from glob import glob

def asl_convert_dicom_to_nifti(subject):

    # Convert DICOM files to NIfTI format using dcm2niix.
    # This function performs the following steps:
    # 1. Copies all DICOM files from the input directory to the subject-specific DICOM directory.
    # 2. Creates an 'ORIG' subdirectory within the subject DICOM directory for original DICOM files.
    # 3. Runs dcm2niix to rename DICOM files for readability.
    # 4. Removes files matching specific patterns from the subject DICOM directory.
    # 5. Moves files matching certain patterns to the 'ORIG' directory.
    # 6. Renames files in the subject DICOM directory to remove unwanted characters.
    # 7. Runs dcm2niix again to convert DICOM files to compressed NIfTI format in the output directory.
    # 8. Performs final cleanup by renaming files in the NIfTI output directory and removing any raw files.
    # Parameters:
    #    - subject: Dictionary containing paths for DICOM and NIfTI directories.
    # Returns:
    #    - subject: Updated subject dictionary with processed paths.

    nifti_output_dir = subject['NIFTIdir'] # Directory where NIfTI files will be saved
    dicom_input_dir = subject['DICOMinputdir']  # Directory containing the original DICOM files
    dicom_subject_dir = subject['DICOMsubjectdir'] # Directory for original DICOM files copied to subject directory

    def run_command(cmd):
        logging.info(f"Running command: {cmd}")
        subprocess.run(cmd, shell=True, check=False)

    def rename_files(directory):
        for file in os.listdir(directory):
            old_path = os.path.join(directory, file)
            if os.path.isfile(old_path):
                new_name = file.replace("-_", "").replace("-", "")
                new_path = os.path.join(directory, new_name)
                if old_path != new_path:
                    os.rename(old_path, new_path)

    def remove_files_by_pattern(directory, patterns):
        for pattern in patterns:
            for file_path in glob(os.path.join(directory, pattern)):
                try:
                    os.remove(file_path)
                except FileNotFoundError:
                    pass

    def move_files_by_pattern(src_dir, dst_dir, patterns):
        for pattern in patterns:
            for file_path in glob(os.path.join(src_dir, pattern)):
                try:
                    shutil.move(file_path, dst_dir)
                except (FileNotFoundError, shutil.Error):
                    pass

    # Copy all DICOMs dicom_input_dir into dicom_subject_dir
    run_command(f'cp {dicom_input_dir}/* {dicom_subject_dir}')
        
    # --- Prepare directories ---
    orig_dir = os.path.join(dicom_subject_dir, 'ORIG')
    os.makedirs(orig_dir, exist_ok=True)

    logging.info('Converting DICOMs to NIFTI using dcm2niix (initial rename)')
    # Run dcm2niix to rename DICOM files
    # This step is to ensure that the DICOM files are renamed in a readable format before further processing.
    run_command(f'dcm2niix -w 0 -r y -f %p_%s {dicom_subject_dir} > {os.path.join(nifti_output_dir, "dcm2niix_rename.log")} 2>&1')

    # Clean and organize files
    remove_files_by_pattern(dicom_subject_dir, ['*._Raw', '*_PS'])
    move_files_by_pattern(dicom_subject_dir, orig_dir, ['IM_*', 'XX_*', 'PS_*'])
    rename_files(dicom_subject_dir)

    # --- Main conversion ---
    logging.info('Converting DICOMs to NIFTI using dcm2niix (final conversion)')
    # Run dcm2niix to convert DICOM files to NIfTI format
    # This step compresses the NIfTI files and organizes them in the specified output directory.
    run_command(f'dcm2niix -w 1 -z y -b y -f %p_%s -o {nifti_output_dir} {dicom_subject_dir} > {os.path.join(nifti_output_dir, "dcm2niix.log")} 2>&1')
    logging.info('DICOMs converted to NIFTI')

    # Final cleanup
    rename_files(nifti_output_dir)
    remove_files_by_pattern(nifti_output_dir, ['*Raw*'])
    return subject  # subject dictionary with updated paths


