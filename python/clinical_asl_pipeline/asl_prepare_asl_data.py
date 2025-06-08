import numpy as np
import os
import nibabel as nib
from clinical_asl_pipeline.asl_interleave_control_tag import  asl_interleave_control_tag
from clinical_asl_pipeline.save_data_nifti import  save_data_nifti

def asl_prepare_asl_data(subject, filename, prefix, fast):
    
    #Prepare (multidelay) ASL data: interleave control and label files per PLD, M0, and perform Look-Locker Correction.
   
    # Load data
    data_path = os.path.join(subject['NIFTIdir'], filename)
    data = nib.load(data_path).get_fdata()/nib.load(data_path).dataobj.slope # remove scaling, as nibabel nib.load consumes the slope and intcept automatically
    print(subject['dummyfilenameSaveNII'])
    
    dims = data.shape

    # Initialize arrays
    subject[prefix] = {}
    subject[prefix]['M0ASL_allPLD'] = np.zeros((dims[0], dims[1], dims[2], subject['NDYNS'], subject['NPLDS'], 2))
    subject[prefix]['ASL_label1label2_allPLD'] = np.zeros((dims[0], dims[1], dims[2], subject['NREPEATS'] * 2, subject['NPLDS']))

    # Split label/control, apply Look Locker correction, save individual PLD files
    for i in range(subject['NPLDS']):
        idx_label = slice(i, subject['NPLDS'] * subject['NDYNS'] * 2, 2 * subject['NPLDS'])
        idx_control = slice(i + subject['NPLDS'], subject['NPLDS'] * subject['NDYNS'] * 2, 2 * subject['NPLDS'])

        subject[prefix]['M0ASL_allPLD'][:, :, :, :subject['NDYNS'], i, 0] = data[:, :, :, idx_label] / subject['LookLocker_correction_factor_perPLD'][i]
        subject[prefix]['M0ASL_allPLD'][:, :, :, :subject['NDYNS'], i, 1] = data[:, :, :, idx_control] / subject['LookLocker_correction_factor_perPLD'][i]

        save_data_nifti(np.squeeze(subject[prefix]['M0ASL_allPLD'][:, :, :, :, i, 0]), os.path.join(subject['ASLdir'], f"{prefix}_PLD0{i+1}_label1.nii.gz"), subject['dummyfilenameSaveNII'], 1, None, subject['TR'])
        save_data_nifti(np.squeeze(subject[prefix]['M0ASL_allPLD'][:, :, :, :, i, 1]), os.path.join(subject['ASLdir'], f"{prefix}_PLD0{i+1}_label2.nii.gz"), subject['dummyfilenameSaveNII'], 1, None, subject['TR'])

    # Interleave control/label, save per-PLD
    for i in range(subject['NPLDS']):
        interleaved = asl_interleave_control_tag(np.squeeze(subject[prefix]['M0ASL_allPLD'][:, :, :, 1:, i, 0]), np.squeeze(subject[prefix]['M0ASL_allPLD'][:, :, :, 1:, i, 1]))
        subject[prefix]['ASL_label1label2_allPLD'][:, :, :, :, i] = interleaved
        save_data_nifti(interleaved, os.path.join(subject['ASLdir'], f"{prefix}_PLD0{i+1}_label1label2.nii.gz"), subject['dummyfilenameSaveNII'], 1, None, subject['TR'])

    print("Saving ASL data interleaved label control: all PLDs for AAT")
    # Swap axes to interleave PLDs and time correctly
    reordered = np.transpose(subject[prefix]['ASL_label1label2_allPLD'], (0, 1, 2, 4, 3))  # (x, y, z, PLD, time)
    # Now reshape so time dimension becomes interleaved PLDs
    data_allPLD = reordered.reshape(dims[0], dims[1], dims[2], subject['NPLDS'] * subject['NREPEATS'] * 2)
    save_data_nifti(data_allPLD, os.path.join(subject['ASLdir'], f"{prefix}_allPLD_label1label2.nii.gz"), subject['dummyfilenameSaveNII'], 1, None, subject['TR'])

    print("Saving ASL data interleaved label control: 2-to-last PLDs for CBF")
    # Extract PLDs 2 to end → shape: (x, y, z, time, NPLDS-1)
    data_2tolastPLD = subject[prefix]['ASL_label1label2_allPLD'][:, :, :, :, 1:]
    # Reorder to (x, y, z, PLD, time), then reshape
    data_2tolastPLD = np.transpose(data_2tolastPLD, (0, 1, 2, 4, 3)).reshape(dims[0], dims[1], dims[2], (subject['NPLDS'] - 1) * subject['NREPEATS'] * 2)
    save_data_nifti(data_2tolastPLD, os.path.join(subject['ASLdir'], f"{prefix}_2tolastPLD_label1label2.nii.gz"), subject['dummyfilenameSaveNII'], 1, None, subject['TR'])
    
    print("Saving ASL data interleaved label control: 1-to-2 PLDs for ATA")
    # Extract first 2 PLDs → shape: (x, y, z, time, 2)
    data_1to2PLD = subject[prefix]['ASL_label1label2_allPLD'][:, :, :, :, 0:2]    
    # Reorder to (x, y, z, PLD, time), then reshape
    data_1to2PLD = np.transpose(data_1to2PLD, (0, 1, 2, 4, 3)).reshape(dims[0], dims[1], dims[2], 2 * subject['NREPEATS'] * 2)
    save_data_nifti(data_1to2PLD, os.path.join(subject['ASLdir'], f"{prefix}_1to2PLD_label1label2.nii.gz"), subject['dummyfilenameSaveNII'], 1, None, subject['TR'])
    
    # M0 image construction
    subject[prefix]['M0_allPLD'] = np.mean(subject[prefix]['M0ASL_allPLD'][:, :, :, 0, :, :], axis=4)
    subject[prefix]['M0'] = subject[prefix]['M0_allPLD'][:, :, :, 0]

    if fast != 'fast':
        save_data_nifti(subject[prefix]['M0_allPLD'], os.path.join(subject['ASLdir'], f"{prefix}_M0_allPLD.nii.gz"), subject['dummyfilenameSaveNII'], 1, None, subject['TR'])

    print("Saving M0 image")
    save_data_nifti(subject[prefix]['M0'], os.path.join(subject['ASLdir'], f"{prefix}_M0.nii.gz"), subject['dummyfilenameSaveNII'], 1, None, subject['TR'])

    return subject
