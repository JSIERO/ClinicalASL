import ants # type: ignore

def asl_motioncorrection_ants(inputdata, refdata, outputdata):
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
