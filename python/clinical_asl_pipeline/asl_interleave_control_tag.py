"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

ASL control/tag interleaving module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Function for interleaving control and tag ASL MRI data volumes.

License: BSD 3-Clause License
"""

import numpy as np

def asl_interleave_control_tag(ctrl, tag=None):
    # Interleave control and tag ASL data for analysis.
    # This function takes ASL data from control and tag volumes and interleaves them
    # such that the output has alternating control and tag volumes.
    # Interleave CTRL and TAG ASL data:
    # CTRL1, TAG1, CTRL2, TAG2, ...
    # Parameters:
    # - ctrl: 5D numpy array (x, y, z, time, 2) if tag is None
    #        or (x, y, z, time) if tag is provided separately
    #- tag:  Optional. Same shape as ctrl if provided.
    #
    # Returns:
    #- output: 4D numpy array (x, y, z, time*2) with interleaved volumes
   
    dims = ctrl.shape

    # Determine interleaving indices
    num_volumes = dims[3]
    odd_index = np.arange(0, 2 * num_volumes, 2)
    even_index = np.arange(1, 2 * num_volumes, 2)

    if tag is None:
        data = ctrl
        if data.ndim < 5:
            print("WARNING: probably single-slice data!")

        # Preallocate interleaved array
        output = np.zeros(data.shape[:3] + (2 * num_volumes,), dtype=data.dtype)
        output[:, :, :, odd_index] = data[:, :, :, :, 0]
        output[:, :, :, even_index] = data[:, :, :, :, 1]
    else:
        # Preallocate interleaved array
        output = np.zeros(ctrl.shape[:3] + (2 * num_volumes,), dtype=ctrl.dtype)
        output[:, :, :, odd_index] = ctrl
        output[:, :, :, even_index] = tag

    return output
