"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

ASL Look-Locker correction module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Function for computing Look-Locker correction factors for ASL MRI data.

License: BSD 3-Clause License
"""

import logging
import numpy as np

def asl_look_locker_correction(subject, phase_tag):
    # Compute Look-Locker scaling correction per PLD for mDelay PCASL, by JCW SIERO 2020.
    # This function computes the Look-Locker correction factor for each PLD in the ASL data.
    # Parameters:   
    # subject: dictionary containing ASL parameters and data
    # phase_tag: string indicating the phase/tag for the subject's data, e.g., 'preACZ', 'postACZ', etc.
    # Returns: subject dictionary with Look-Locker correction factor for each PLD, to be used in ALS quantification (QASL)
    # Ensure subject dictionary has required keys
    # required_keys = ['FLIPANGLE', 'PLDS', 'T1b']s, ie correction only depends on flip angle, PLDs timing, and T1b

    # Flip angle in degrees
    flip_angle = subject[phase_tag]['FLIPANGLE']

    # Assume M0 = 1
    M0 = 1

    # PLDs in milliseconds
    PLD = np.array(subject[phase_tag]['PLDS']) * 1000  # convert from s to ms

    # Time vector t
    t = np.arange(1, int(PLD[-1] * 1.1) + 1)  # MATLAB-style 1:PLD(end)*1.1

    # Blood T1 in ms
    T1_blood = subject['T1b'] * 1e3  # 1.65 * 1000

    # Average delta PLD
    delta_PLD = np.mean(np.diff(PLD))

    # Effective T1 during Look-Locker
    T1_eff_blood = 1 / (1 / T1_blood - np.log(np.cos(np.radians(flip_angle))) / delta_PLD)

    # Equilibrium magnetization during LL readout, # Brix et al MRI 1990
    Meq_LL_blood = M0 * (1 - np.exp(-delta_PLD / T1_blood)) / (1 - np.cos(np.radians(flip_angle)) * np.exp(-delta_PLD / T1_blood))

    # Magnetization without Look-Locker
    Mz_noLL_C = M0 * (1 - np.exp(-t / T1_blood)) # control data
    Mz_noLL_L = M0 * (1 - 2 * np.exp(-t / T1_blood)) # label data
    deltaM_noLL = Mz_noLL_C - Mz_noLL_L      

    # With Look-Locker
    Mz_LL_C = np.zeros_like(t, dtype=np.float64)
    Mz_LL_L = np.zeros_like(t, dtype=np.float64)

    t1 = np.arange(1, int(PLD[0]) + 1)
    Mz_LL_C[t1 - 1] = M0 * (1 - np.exp(-t1 / T1_blood))
    Mz_LL_L[t1 - 1] = M0 * (1 - 2 * np.exp(-t1 / T1_blood))

    t2 = np.arange(int(PLD[0]) + 1, int(t[-1]) + 1)
    offset = int(round(PLD[0]))
    decay = np.exp(-t[t2 - 1 - offset] / T1_eff_blood)
    recovery = Meq_LL_blood * (1 - decay)
    Mz_LL_C[t2 - 1] = Mz_LL_C[offset - 1] * decay + recovery
    Mz_LL_L[t2 - 1] = Mz_LL_L[offset - 1] * decay + recovery

    deltaM_LL = Mz_LL_C - Mz_LL_L

    with np.errstate(divide='ignore', invalid='ignore'):
        deltaM_ratio_LL_noLL = deltaM_LL / deltaM_noLL
        deltaM_ratio_LL_noLL[~np.isfinite(deltaM_ratio_LL_noLL)] = 0
  
    # Convert to Mxy using flip angle
    deltaMxy_ratio_LL_noLL = deltaM_ratio_LL_noLL * np.sin(np.radians(flip_angle))

    # Extract values at PLD time points (rounded to nearest integer indices)
    tpointsMxy = PLD.astype(int)
    LookLocker_correction_factor_perPLD = np.round(deltaMxy_ratio_LL_noLL[tpointsMxy - 1],3)  # -1 for 0-based index
    
    # Display summary
    logging.info(
        f"Look-Locker correction factor for flipangle = {flip_angle} (deg), "
        f"PLDs(ms): {PLD.tolist()}, and deltaPLD(ms) = {delta_PLD}: "
        f"{LookLocker_correction_factor_perPLD}, JCW SIERO 2020"
    )
    subject[phase_tag]['LookLocker_correction_factor_perPLD'] = LookLocker_correction_factor_perPLD
    return subject
