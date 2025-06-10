import os
import logging
import shutil
import ants

def asl_registration_prepostACZ_ANTS(subject):
    # Register post-ACZ ASL data to pre-ACZ ASL data using ANTsPy
    # subject is a dictionary containing paths to the necessary files.
    # The dictionary should contain the following keys:
    # 'preACZ': {
    #     'T1fromM0_path', 'CBF_path', 'AAT_path', 'ATA_path', 'mask_path'
    # },
    # 'postACZ': {
    #     'T1fromM0_path', 'CBF_path', 'AAT_path', 'ATA_path', 'mask_path',
    #     'T1fromM0_2preACZ_path', 'CBF_2preACZ_path', 'AAT_2preACZ_path',
    #     'ATA_2preACZ_path', 'mask_2preACZ_path'
    # },
    # 'ASLdir'
    # resulting transform will be saved in 'ASLdir' as 'rigid_postACZ_to_preACZ.mat'

    # Load fixed and moving images for registration
    logging.info("Registration T1fromM0 postACZ to preACZ data (ANTsPy) *********************************************************************")

    fixed = ants.image_read(subject['preACZ']['T1fromM0_path'])
    moving = ants.image_read(subject['postACZ']['T1fromM0_path'])

    # Run registration
    reg = ants.registration(fixed=fixed, moving=moving, type_of_transform='Rigid', metric='Mattes', reg_iterations=(1000, 500, 250, 100)) 
    interpolator = 'bSpline'
    # Save transformed moving image
    ants.image_write(reg['warpedmovout'], subject['postACZ']['T1fromM0_2preACZ_path'])

    # Apply same transform to CBF, AAT, ATA, mask
    def apply_transform(moving_path, reference_path, output_path, transformlist, interpolation):
        moving_img = ants.image_read(moving_path)
        reference_img = ants.image_read(reference_path)
        warped = ants.apply_transforms(fixed=reference_img, moving=moving_img,
                                        transformlist=transformlist,
                                        interpolator=interpolation)
        ants.image_write(warped, output_path)

    # CBF
    logging.info("Registration CBF postACZ to preACZ (ANTsPy)")
    apply_transform(subject['postACZ']['QASL_CBF_path'], subject['preACZ']['QASL_CBF_path'], subject['postACZ']['CBF_2preACZ_path'], reg['fwdtransforms'], interpolator)

    # AAT
    logging.info("Registration AAT postACZ to preACZ (ANTsPy)")
    apply_transform(subject['postACZ']['QASL_AAT_path'], subject['preACZ']['QASL_AAT_path'], subject['postACZ']['AAT_2preACZ_path'], reg['fwdtransforms'], interpolator)

    # ATA
    logging.info("Registration ATA postACZ to preACZ (ANTsPy)")
    apply_transform(subject['postACZ']['QASL_ATA_path'], subject['preACZ']['QASL_ATA_path'], subject['postACZ']['ATA_2preACZ_path'], reg['fwdtransforms'], interpolator)

    # Mask (NearestNeighbor interpolation)
    logging.info("Registration mask postACZ to preACZ (ANTsPy)")
    apply_transform(subject['postACZ']['mask_path'], subject['preACZ']['mask_path'], subject['postACZ']['mask_2preACZ_path'], reg['fwdtransforms'], 'nearestNeighbor')

    transform_path = reg['fwdtransforms'][0]

    # Choose your save location for transform (ITK format)
    output_transform_path = os.path.join(subject['ASLdir'], 'rigid_postACZ_to_preACZ.mat')

    # Copy it to save location
    shutil.copy(transform_path, output_transform_path)