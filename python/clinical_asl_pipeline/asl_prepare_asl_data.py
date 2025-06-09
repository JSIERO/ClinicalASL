import numpy as np
import os
import nibabel as nib
from clinical_asl_pipeline.asl_interleave_control_tag import asl_interleave_control_tag
from clinical_asl_pipeline.utils.save_data_nifti import save_data_nifti

def asl_prepare_asl_data(subject, filename, prefix):
    # Prepare (multidelay ASL data for analysis by interleaving control and label files per PLD, M0, and performing Look-Locker Correction.
    # Parameters:
    # subject: dict containing subject information including paths and parameters
    # filename: string, name of the NIfTI file containing ASL data
    # prefix: string, prefix for the keys in the subject dictionary to store results
    # Returns:
    # subject: updated subject dictionary with processed ASL data   
  
   
    # Load data
    data_path = os.path.join(subject['NIFTIdir'], filename)
    data = nib.load(data_path).get_fdata()/nib.load(data_path).dataobj.slope # remove scaling, as nibabel nib.load consumes the slope and intcept automatically
    print(subject[f'{prefix}_templateNII_path'])
    
    dims = data.shape

    # Initialize arrays
    subject[prefix] = {}
    subject[prefix]['M0ASL_allPLD'] = np.zeros((dims[0], dims[1], dims[2], subject['NDYNS'], subject['NPLDS'], 2))
    subject[prefix]['ASL_label1label2_allPLD'] = np.zeros((dims[0], dims[1], dims[2], subject['NREPEATS'] * 2, subject['NPLDS']))

    # Split label/control, apply Look Locker correction
    for i in range(subject['NPLDS']):
        idx_label = slice(i, subject['NPLDS'] * subject['NDYNS'] * 2, 2 * subject['NPLDS'])
        idx_control = slice(i + subject['NPLDS'], subject['NPLDS'] * subject['NDYNS'] * 2, 2 * subject['NPLDS'])

        subject[prefix]['M0ASL_allPLD'][:, :, :, :subject['NDYNS'], i, 0] = data[:, :, :, idx_label] / subject['LookLocker_correction_factor_perPLD'][i]
        subject[prefix]['M0ASL_allPLD'][:, :, :, :subject['NDYNS'], i, 1] = data[:, :, :, idx_control] / subject['LookLocker_correction_factor_perPLD'][i]

    # Interleave control/label, save per-PLD
    for i in range(subject['NPLDS']):
        interleaved = asl_interleave_control_tag(np.squeeze(subject[prefix]['M0ASL_allPLD'][:, :, :, 1:, i, 0]), np.squeeze(subject[prefix]['M0ASL_allPLD'][:, :, :, 1:, i, 1]))
        subject[prefix]['ASL_label1label2_allPLD'][:, :, :, :, i] = interleaved

    print("Saving ASL data interleaved label control: all PLDs for AAT")
    # Swap axes to interleave PLDs and time correctly
    reordered = np.transpose(subject[prefix]['ASL_label1label2_allPLD'], (0, 1, 2, 4, 3))  # (x, y, z, PLD, time)    
    # Now reshape so time dimension becomes interleaved PLDs
    PLDall = reordered.reshape(dims[0], dims[1], dims[2], subject['NPLDS'] * subject['NREPEATS'] * 2)
    save_data_nifti(PLDall, subject[f'{prefix}_PLDall_labelcontrol_path'], subject[f'{prefix}_templateNII_path'], 1, None, subject['TR'])

    print("Saving ASL data interleaved label control: 2-to-last PLDs for CBF")
    # Extract PLDs 2 to end → shape: (x, y, z, time, NPLDS-1)
    PLD2tolast = subject[prefix]['ASL_label1label2_allPLD'][:, :, :, :, 1:]
    # Reorder to (x, y, z, PLD, time), then reshape
    PLD2tolast = np.transpose(PLD2tolast, (0, 1, 2, 4, 3)).reshape(dims[0], dims[1], dims[2], (subject['NPLDS'] - 1) * subject['NREPEATS'] * 2)
    save_data_nifti(PLD2tolast, subject[f'{prefix}_PLD2tolast_labelcontrol_path'], subject[f'{prefix}_templateNII_path'], 1, None, subject['TR'])
    
    print("Saving ASL data interleaved label control: 1-to-2 PLDs for ATA")
    # Extract first 2 PLDs → shape: (x, y, z, time, 2)
    PLD1to2 = subject[prefix]['ASL_label1label2_allPLD'][:, :, :, :, 0:2]    
    # Reorder to (x, y, z, PLD, time), then reshape
    PLD1to2 = np.transpose(PLD1to2, (0, 1, 2, 4, 3)).reshape(dims[0], dims[1], dims[2], 2 * subject['NREPEATS'] * 2)
    save_data_nifti(PLD1to2, subject[f'{prefix}_PLD1to2_labelcontrol_path'], subject[f'{prefix}_templateNII_path'], 1, None, subject['TR'])
    
    # M0 image construction
    subject[prefix]['M0_allPLD'] = np.mean(subject[prefix]['M0ASL_allPLD'][:, :, :, 0, :, :], axis=4)
    subject[prefix]['M0'] = subject[prefix]['M0_allPLD'][:, :, :, 0]

    print("Saving M0 image")
    save_data_nifti(subject[prefix]['M0'], subject[f'{prefix}_M0_path'] , subject[f'{prefix}_templateNII_path'], 1, None, subject['TR'])

    return subject
