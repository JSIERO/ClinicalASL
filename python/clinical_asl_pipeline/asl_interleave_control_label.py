"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

ASL control/label interleaving module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Function for interleaving control and label ASL MRI data volumes.

License: BSD 3-Clause License
"""

import warnings
import numpy as np

def asl_interleave_control_label(ctrl, label=None):
    # Interleave control and label ASL data for analysis.
    # This function takes ASL data from control and label volumes and interleaves them
    # such that the output has alternating control and label volumes.
    # Interleave CTRL and LABEL ASL data:
    # CTRL1, LABEL1, CTRL2, LABEL2, ...
    # Parameters:
    # - ctrl: 5D numpy array (x, y, z, time, 2) if label is None
    #        or (x, y, z, time) if label is provided separately
    #- label:  Optional. Same shape as ctrl if provided.
    #
    # Returns:
    #- output: 4D numpy array (x, y, z, time*2) with interleaved volumes
    dims = ctrl.shape

    # Determine interleaving indices
    num_volumes = dims[3]
    odd_index = np.arange(0, 2 * num_volumes, 2)
    even_index = np.arange(1, 2 * num_volumes, 2)

    if label is None:
        data = ctrl
        if data.ndim < 5:
            warnings.warn("WARNING: probably single-slice data!")

        # Preallocate interleaved array
        output = np.zeros(data.shape[:3] + (2 * num_volumes,), dtype=data.dtype)
        output[:, :, :, odd_index] = data[:, :, :, :, 0]
        output[:, :, :, even_index] = data[:, :, :, :, 1]
    else:
        # Preallocate interleaved array
        output = np.zeros(ctrl.shape[:3] + (2 * num_volumes,), dtype=ctrl.dtype)
        output[:, :, :, odd_index] = ctrl
        output[:, :, :, even_index] = label

    return output
