import ants # type: ignore

def asl_motioncorrection_ants(inputdata, refdata, outputdata):
    # Perform motion correction using ANTs
    # inputdata: path to the input ASL data in NIfTI format
    # refdata: path to the reference image for motion correction
    # outputdata: path to save the motion-corrected output data in NIfTI format
    print('Perform motion correction')
    results_dict = ants.motion_correction(
        ants.image_read(inputdata),
        ants.image_read(refdata),
        type_of_transform="DenseRigid",
        aff_metric="mattes",
        smoothing_in_mm=True,
        singleprecision=True,
        verbose=True
    )
    ants.image_write(results_dict["motion_corrected"], outputdata)
