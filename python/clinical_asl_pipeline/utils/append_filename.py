"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

ASL Look-Locker correction module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Little utility functions for ClinicalASL, such as appending '_mc' to filenames for motion corrected files.
    Little utility functions for ClinicalASL, such as appending '_or' to filenames for outlier removed files.


License: BSD 3-Clause License
"""

def append_mc(filename):
    # Append '_mc' at the end of the filename .nii.gz in filename (for motion corrected files)."""
    if filename.endswith('.nii.gz'):
        return filename[:-7] + '_mc.nii.gz'
    return filename


def append_or(filename):
    # Append '_or' before .nii.gz in filename (for outlier removedd files)."""
    if filename.endswith('.nii.gz'):
        return filename[:-7] + '_or.nii.gz'
    return filename