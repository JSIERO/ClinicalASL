import os
import shutil
import subprocess
from glob import glob

def asl_convert_dicom_to_nifti(dicom_input_dir, nifti_output_dir, vtr_only=None, imager=None, subject_dir=None):
    """
    Updated version to match MATLAB code.
    If imager == 'IMAGER', will copy DICOMs to SUBJECTdir/DICOMORIG and process from there.

    Args:
        dicom_input_dir: original input dir
        nifti_output_dir: output dir for NIfTI
        vtr_only: 'vTR_only' or None
        imager: 'IMAGER' or None
        subject_dir: SUBJECTdir full path (required if imager == 'IMAGER')
    Returns:
    - output: DICOMORIG directoy with renamed DICOMS for further processing
    """

    def run_command(cmd):
        print(f"Running command: {cmd}")
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

    # --- Handle IMAGER case ---
    if imager == 'IMAGER':
        if subject_dir is None:
            raise ValueError("subject_dir must be provided if imager == 'IMAGER'")

        dicom_orig_dir = os.path.join(subject_dir, 'DICOMORIG')

        # Make DICOMORIG if not exists
        os.makedirs(dicom_orig_dir, exist_ok=True)

        # Copy all DICOMs into DICOMORIG
        run_command(f'cp {dicom_input_dir}/* {dicom_orig_dir}')
        
    # Now common logic

    orig_dir = os.path.join(dicom_orig_dir, 'ORIG')
    os.makedirs(orig_dir, exist_ok=True)

    # Initial DICOM conversion (rename phase)
    run_command(f'dcm2niix -w 0 -r y -f %p_%s {dicom_orig_dir} > {os.path.join(nifti_output_dir, "dcm2niix_rename.log")} 2>&1')

    # Clean and organize files
    remove_files_by_pattern(dicom_orig_dir, ['*._Raw', '*_PS'])
    move_files_by_pattern(dicom_orig_dir, orig_dir, ['IM_*', 'XX_*', 'PS_*'])
    rename_files(dicom_orig_dir)

    # --- Main conversion ---
    if vtr_only is None:
        run_command(f'dcm2niix -w 1 -z y -b y -f %p_%s -o {nifti_output_dir} {dicom_orig_dir} > {os.path.join(nifti_output_dir, "dcm2niix.log")} 2>&1')
    elif vtr_only == 'vTR_only':
        # Handle vTR ASL files
        files_vtr = glob(os.path.join(dicom_input_dir, '*PRIDE*SOURCE*vTR*'))
        if files_vtr:
            remove_files_by_pattern(nifti_output_dir, ['PRIDE*SOURCE*vTR*.*'])
            for file_path in files_vtr:
                run_command(f'dcm2niix -w 2 -z y -b y -f %p_%s -s y -o {nifti_output_dir} {file_path} > {os.path.join(nifti_output_dir, "dcm2niix.log")} 2>&1')

        # Handle SOURCE M0 files
        files_m0 = glob(os.path.join(dicom_input_dir, '*SOURCE*M0*'))
        if files_m0:
            remove_files_by_pattern(nifti_output_dir, ['*SOURCE*M0*.*'])
            for file_path in files_m0:
                run_command(f'dcm2niix -w 2 -z y -b y -f %p_%s -s y -o {nifti_output_dir} {file_path} > {os.path.join(nifti_output_dir, "dcm2niix.log")} 2>&1')
    
    # Final cleanup
    rename_files(nifti_output_dir)
    remove_files_by_pattern(nifti_output_dir, ['*Raw*'])
    return dicom_orig_dir  # update DICOM input directory to DICOMORIG directory for further processing


