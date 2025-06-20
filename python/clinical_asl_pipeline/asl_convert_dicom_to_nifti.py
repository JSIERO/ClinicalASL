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
from clinical_asl_pipeline.utils.run_command_with_logging import run_command_with_logging
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
    dcmniixlog_dir = subject['SUBJECTdir'] # # Directory for log files dcm2niix for rename and dicom conversiopn 

    def run_command(cmd):
        logging.info(f"Running command: {cmd}")
        run_command_with_logging(cmd)

    def copy_dicom_series(dicom_input_dir, dicom_subject_dir):
        # Check read access to input directory
        if not os.access(dicom_input_dir, os.R_OK):
            logging.error(f"No read access to input directory: {dicom_input_dir}")
            raise PermissionError(f"Cannot read from {dicom_input_dir}")

        # Create destination directory if it doesn't exist
        os.makedirs(dicom_subject_dir, exist_ok=True)

        # Check write access to destination directory
        if not os.access(dicom_subject_dir, os.W_OK):
            logging.error(f"No write access to output directory: {dicom_subject_dir}")
            raise PermissionError(f"Cannot write to {dicom_subject_dir}")

        logging.info(f"Copying DICOM files from {dicom_input_dir} to {dicom_subject_dir}")

        # Copy files
        copied_count = 0
        for fname in os.listdir(dicom_input_dir):
            src = os.path.join(dicom_input_dir, fname)
            dst = os.path.join(dicom_subject_dir, fname)
            if os.path.isfile(src):
                try:
                    shutil.copy2(src, dst)
                    logging.info(f"Copied: {src} → {dst}")
                    copied_count += 1
                except Exception as e:
                    logging.error(f"Failed to copy {src}: {e}")
            else:
                logging.warning(f"Skipped (not a file): {src}")

        logging.info(f"Done. Total files copied: {copied_count}")

    def rename_files(directory):
        for file in os.listdir(directory):
            old_path = os.path.join(directory, file)
            if os.path.isfile(old_path):
                new_name = file.replace("-_", "").replace("-", "")
                new_path = os.path.join(directory, new_name)
                if old_path != new_path:
                    os.replace(old_path, new_path)

    def remove_files_by_pattern(directory, patterns):
        for pattern in patterns:
            for file_path in glob(os.path.join(directory, pattern)):
                try:
                    os.remove(file_path)
                except FileNotFoundError:
                    pass

    def move_files_by_pattern(src_dir, dst_dir, patterns):
        os.makedirs(dst_dir, exist_ok=True)  # Ensure target dir exists
        for pattern in patterns:
            for file_path in glob(os.path.join(src_dir, pattern)):
                filename = os.path.basename(file_path)
                dst_path = os.path.join(dst_dir, filename)
                try:
                    os.replace(file_path, dst_path)
                except FileNotFoundError:
                    logging.warning(f"File not found (skipped): {file_path}")
                except PermissionError:
                    logging.error(f"Permission error moving {file_path} -> {dst_path}")

    # Copy all DICOMs from dicom_input_dir into dicom_subject_dir
    copy_dicom_series(dicom_input_dir, dicom_subject_dir)
        
    # --- Prepare directories ---
    orig_dir = os.path.join(dicom_subject_dir, 'ORIG')
    os.makedirs(orig_dir, exist_ok=True)

    logging.info('Converting DICOMs to NIFTI using dcm2niix v1.0.20220720 (initial rename)')
    # Run dcm2niix to rename DICOM files
    # This step is to ensure that the DICOM files are renamed in a readable format before further processing.
    run_command(f'dcm2niix -v 1 -w 0 -r y -f %p_%s {dicom_subject_dir} > {os.path.join(dcmniixlog_dir, "dcm2niix_rename.log")} 2>&1')

    # Clean and organize files
    remove_files_by_pattern(dicom_subject_dir, ['*_Raw', '*_PS'])
    move_files_by_pattern(dicom_subject_dir, orig_dir, ['IM_*', 'XX_*', 'PS_*'])
    rename_files(dicom_subject_dir)

    # --- Main conversion ---
    logging.info('Converting DICOMs to NIFTI using dcm2niix v1.0.20220720 (final conversion)')
    # Run dcm2niix to convert DICOM files to NIfTI format
    # This step compresses the NIfTI files and organizes them in the specified output directory.
    run_command(f'dcm2niix -w 1 -z y -b y -f %p_%s -o {nifti_output_dir} {dicom_subject_dir} > {os.path.join(dcmniixlog_dir, "dcm2niix.log")} 2>&1')
    logging.info('DICOMs converted to NIFTI')

    # Final cleanup
    rename_files(nifti_output_dir)
    return subject  # subject dictionary with updated paths


