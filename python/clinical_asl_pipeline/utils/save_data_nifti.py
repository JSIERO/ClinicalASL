"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

NIfTI saving utility module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Utility function to save data as NIfTI files using header information from a template file.

License: BSD 3-Clause License
"""

import numpy as np
import nibabel as nib

def save_data_nifti(data, output_filename, templateNII_filename, scaleslope, datarange=None, TR=None):
    # Save data to a NIfTI file using header information from a reference (dummy) file.
    #
    # Parameters:
    #   data           : numpy array
    #                    The image data to be saved.
    #   output_filename: str
    #                    Path to save the new NIfTI file.
    #   dummy_filename : str
    #                    Reference NIfTI file from which to copy the header and affine.
    #   scaleslope     : float or 'samescaling'
    #                    Value to set for the NIfTI scl_slope field, or 'samescaling' to keep original.
    #   datarange      : tuple (min, max), optional
    #                    Intensity range for display (sets cal_min and cal_max in header).
    #   TR             : float, optional
    #                    Repetition time in seconds (sets pixdim[4] in header).
    #
    # Returns:
    #   None. The function saves the data to the specified output_filename.
    #
    # Notes:
    #   - The function copies the header and affine from the dummy file.
    #   - If data is floating point, it is saved as float64.
    #   - If scaleslope is not 'samescaling', scl_slope is set to the provided value.
    #   - If the scaling intercept (scl_inter) is not zero, it is reset to zero with a warning.
    #   - The function updates header dimensions and calibration min/max as needed.
    #   - If TR is provided, it is set in the header.
    
    data_info = nib.load(templateNII_filename).header.copy()
    affine = nib.load(templateNII_filename).affine

    if np.issubdtype(data.dtype, np.floating):
        data = data.astype(np.float64)
        data_info.set_data_dtype(np.float64)
        data_info['bitpix'] = 64

    if scaleslope != 'samescaling':
        data_info['scl_slope'] = scaleslope
        
    scl_inter = data_info.get('scl_inter', 0)    
    if not np.isnan(scl_inter) and scl_inter != 0:
        print('WARNING: Scaling intercept is not 0, setting to 0 now')
        data_info['scl_inter'] = 0

    if data.ndim == 4:
        data_info['dim'][0] = 4
        data_info['dim'][4] = data.shape[3]
    else:
        data_info['dim'][0] = 3
        data_info['dim'][4] = 1

    if datarange is None or len(datarange) != 2:
        data_info['cal_min'] = np.nanmin(data)
        data_info['cal_max'] = np.nanmax(data)
    else:
        data_info['cal_min'] = datarange[0]
        data_info['cal_max'] = datarange[1]

    if TR is not None:
        data_info['pixdim'][4] = TR
        # Assume mm and s as default units in nibabel
        # No enforced string field for units in nibabel header

    nib.save(nib.Nifti1Image(data, affine, header=data_info), output_filename)
