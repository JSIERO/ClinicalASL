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
import pydicom
from clinical_asl_pipeline.utils.run_command_with_logging import run_command_with_logging
from glob import glob

def asl_convert_dicom_to_nifti(subject):
# Convert DICOM files to NIfTI format using dcm2niix.
# This function performs the following steps:
# 1. Recursively walks through all files and subdirectories in the input directory.
# 2. Identifies and copies only DICOM files whose Series Description matches user-defined patterns (default: SOURCE*ASL*, SWIP*ASL*).
# 3. Logs skipped files with unmatched Series Descriptions.
# 4. Creates an 'ORIG' subdirectory within the subject DICOM directory for organizing non-ASL and matched ASL DICOM files.
# 5. Uses dcm2niix to rename DICOM files for readability and debugging.
# 6. Removes files matching specific patterns from the subject DICOM directory.
# 7. Moves files matching certain patterns to the 'ORIG' directory.
# 8. Moves matching filtered ASL DICOMs to the 'ORIG' directory.
# 9. Renames files in the subject DICOM directory to remove problematic characters.
# 10. Runs dcm2niix again to convert DICOM files to compressed NIfTI format in the output directory.
# 11. Performs final cleanup by renaming files in the NIfTI output directory.
# Parameters:
#    - subject: Dictionary containing paths for DICOM and NIfTI directories and optional matching patterns.
# Returns:
#    - subject: Updated subject dictionary with processed paths.

    nifti_output_dir = subject['NIFTIdir']
    dicom_input_dir = subject['DICOMinputdir']
    dicom_subject_dir = subject['DICOMsubjectdir']
    dcmniixlog_dir = subject['SUBJECTdir']
    # check for input DICOM series description patterns, default ['SOURCE*ASL*', 'SWIP*ASL*']), case insensitive
    patterns = subject.get('include_dicomseries_description_patterns', ['SOURCE*ASL*', 'SWIP*ASL*'])

    def run_command(cmd):
        logging.info(f"Running command: {cmd}")
        run_command_with_logging(cmd)

    def matches_series_description(desc, patterns):
        import fnmatch
        desc_upper = desc.upper()
        for pattern in patterns:
            if fnmatch.fnmatchcase(desc_upper, pattern.upper()):
                return True
        return False

    def copy_filtered_dicoms(src_dir, dst_dir, patterns):
        if not os.access(src_dir, os.R_OK):
            logging.error(f"No read access to input directory: {src_dir}")
            raise PermissionError(f"Cannot read from {src_dir}")

        os.makedirs(dst_dir, exist_ok=True)

        if not os.access(dst_dir, os.W_OK):
            logging.error(f"No write access to output directory: {dst_dir}")
            raise PermissionError(f"Cannot write to {dst_dir}")

        logging.info(f"Walking input directory recursively: {src_dir}")
        copied_count = 0
        skipped_files = []
        matched_files = []

        for root, _, files in os.walk(src_dir):
            for file in files:
                src_path = os.path.join(root, file)
                try:
                    ds = pydicom.dcmread(src_path, stop_before_pixels=True, force=True)
                    desc = str(ds.get('SeriesDescription', ''))
                    if matches_series_description(desc, patterns):
                        dst_path = os.path.join(dst_dir, os.path.basename(src_path))
                        shutil.copy2(src_path, dst_path)
                        logging.info(f"Copied: {src_path} → {dst_path} [SeriesDescription: {desc}]")
                        copied_count += 1
                        matched_files.append(dst_path)
                    else:
                        skipped_files.append((src_path, desc))
                except Exception as e:
                    logging.error(f"DICOM read error at {src_path}: {e}")
                    raise RuntimeError(f"Critical error reading DICOM file: {src_path}") from e

        for skipped_file, desc in skipped_files:
            logging.info(f"Skipped file (SeriesDescription did not match): {skipped_file} [SeriesDescription: {desc}]")

        if copied_count == 0:
            logging.error("No matching ASL DICOMs found.")
            raise RuntimeError("DICOM conversion aborted: no matching ASL series found.")
        else:
            logging.info(f"Done. Total files copied: {copied_count}")

        return matched_files

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
        os.makedirs(dst_dir, exist_ok=True)
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

    matched_files = copy_filtered_dicoms(dicom_input_dir, dicom_subject_dir, patterns)

    orig_dir = os.path.join(dicom_subject_dir, 'ORIG')
    os.makedirs(orig_dir, exist_ok=True)

    move_files_by_pattern(dicom_subject_dir, orig_dir, ['IM_*', 'XX_*', 'PS_*'])
    # move files for the matched DICOMs from PACS
    for file_path in matched_files:
        try:
            shutil.move(file_path, os.path.join(orig_dir, os.path.basename(file_path)))
        except Exception as e:
            logging.error(f"Failed to move matched file to ORIG: {file_path}: {e}")

    logging.info('Renaming Philips DICOMs using protocolname_seriesnumber filenaming - dcm2niix v1.0.20220720 (initial rename)')
    run_command(f'dcm2niix -v 1 -w 0 -r y -f %p_%s {dicom_subject_dir} > {os.path.join(dcmniixlog_dir, "dcm2niix_rename.log")} 2>&1')

    remove_files_by_pattern(dicom_subject_dir, ['*_Raw', '*_PS'])
    rename_files(dicom_subject_dir)

    logging.info('Converting DICOMs to NIFTI - dcm2niix v1.0.20220720 (final conversion)')
    run_command(f'dcm2niix -w 1 -z y -b y -f %p_%s -o {nifti_output_dir} {dicom_subject_dir} > {os.path.join(dcmniixlog_dir, "dcm2niix_conversion.log")} 2>&1')
    logging.info('DICOMs converted to NIFTI')

    rename_files(nifti_output_dir)
    return subject
