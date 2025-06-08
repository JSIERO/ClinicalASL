import numpy as np
import nibabel as nib

def save_data_nifti(data, output_filename, dummy_filename, scaleslope, datarange=None, TR=None):
    """
    Save data to a NIfTI file using header information from a dummy file.

    Parameters:
    - data: numpy array
    - output_filename: path to save the new NIfTI file
    - dummy_filename: reference NIfTI file for header
    - scaleslope: float or 'samescaling'
    - datarange: optional (min, max) for display intensity range
    - TR: optional repetition time in seconds
    """
    data_info = nib.load(dummy_filename).header.copy()
    affine = nib.load(dummy_filename).affine

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
