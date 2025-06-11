
"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

ASL T1-from-M0 computation module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    
    This module provides a function to dilate a 3D boolean mask slice by slice

License: BSD 3-Clause License
"""

import numpy as np
from scipy.ndimage import binary_dilation, generate_binary_structure

def dilate_mask(mask_3d, mode='2D', iterations=1, conservative=False):
    # Dilate a 3D mask either slice by slice (2D mode) or as full 3D volume.
    #
    # Args:
    #     mask_3d (np.ndarray): 3D boolean mask.
    #     mode (str): '2D' or '3D'.
    #     iterations (int): Number of dilation iterations (default 1).
    #     conservative (bool): If True, use 4-connectivity; else 8-connectivity.
    #
    # Returns:
    #     dilated_mask (np.ndarray): 3D boolean dilated mask.

    # Select connectivity level
    connectivity = 1 if conservative else 2

    # Prepare output mask
    dilated_mask = np.zeros_like(mask_3d, dtype=bool)

    if mode == '2D':
        # 2D dilation: dilate each slice independently
        structure = generate_binary_structure(2, connectivity)

        for z in range(mask_3d.shape[2]):
            slice_mask = mask_3d[:, :, z]
            dilated_slice = binary_dilation(slice_mask, structure=structure, iterations=iterations)
            dilated_mask[:, :, z] = dilated_slice

    elif mode == '3D':
        # 3D dilation: dilate full volume
        structure = generate_binary_structure(3, connectivity)
        dilated_mask = binary_dilation(mask_3d, structure=structure, iterations=iterations)

    else:
        raise ValueError("mode must be '2D' or '3D'")

    return dilated_mask
