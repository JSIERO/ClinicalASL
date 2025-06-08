import subprocess
import time


def asl_basil_analysis(
    subject,
    location_asl_labelcontrol_pld_nifti,
    location_m0,
    location_mask,
    output_map,
    pld_list,
    location_basil_info=None,
    artoff=None,
    spatialoff=None
):

    # Generate comma-separated PLD string
    pld_string = ",".join([f"{pld:.5g}" for pld in pld_list])

    # Optional model options
    if location_basil_info and len(location_basil_info) > 0:
        location_basil_info_string = f" --model-options={location_basil_info}"
    else:
        location_basil_info_string = ""

    # Arterial component off (optional)
    artoff_string = " --artoff" if artoff == "artoff" else ""

    # Spatial regularization (optional)
    spatial_string = "off" if spatialoff == "spatialoff" else "on"

    # Extract parameter values and convert to string
    T1t = str(subject['T1t'])
    T1b = str(subject['T1b'])
    tau = str(subject['tau'])
    TR_M0 = str(subject['TR_M0'][0])
    alpha = str(subject['alpha'])
    slicetime = str(subject['slicetime'] / 1000)  # convert ms to seconds

    # Timing the execution
    start_time = time.time()

    # Build oxford_asl command
    cmd = (
        f"oxford_asl -i {location_asl_labelcontrol_pld_nifti} "
        f"-c {location_m0} -m {location_mask} -o {output_map} "
        f"{location_basil_info_string}{artoff_string} "
        f"--spatial={spatial_string} --bolus={tau} --slicedt={slicetime} "
        f"--t1={T1t} --t1b={T1b} --t1t={T1t} --plds={pld_string} "
        f"--tr={TR_M0} --alpha={alpha} "
        f"--iaf=ct --ibf=tis --casl --fixbolus --cmethod voxel --cgain 1.00"
    )

    # Run command
    print("Running BASIL analysis...")
    subprocess.run(cmd, shell=True, check=True)

    print("BASIL analysis finished")
    elapsed = round(time.time() - start_time, 2)
    print(f"..this took: {elapsed} s")
