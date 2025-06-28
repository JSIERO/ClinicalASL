"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

DICOM saving utility module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Utility function to save data as DICOM files using header information from a template file.

License: BSD 3-Clause License
"""
import os
import logging
import numpy as np
import pydicom
from pydicom.uid import generate_uid
from pydicom.sequence import Sequence
from pydicom.dataset import Dataset
from pydicom.dataelem import DataElement
from pydicom.tag import Tag

def save_data_dicom(image, template_dicom_path, output_dicom_dir, name, value_range, content_label):
# Save a 3D ASL-derived image as either a multiframe or single-frame DICOM series,
#     based on the structure of the provided Philips reference DICOM.

#     This function prepares and scales a quantitative ASL image (e.g., CBF, CVR, AAT, ATA),
#     injects it into a DICOM file or series, and updates relevant metadata fields for 
#     quantitative analysis and PACS compatibility.

#     The output format (multiframe Enhanced MR or single-frame classic DICOM) is automatically 
#     determined from the template DICOM's structure.

#     Parameters:
#         image (np.ndarray): 
#             3D ASL image array with shape (Height, Width, Slices), typically from a NIfTI file.
#         template_dicom_path (str): 
#             Path to a reference DICOM file to use as a DICOM metadata template.
#         output_dicom_dir (str): 
#             Folder path where the output DICOM file (or series folder) will be saved.
#         name (str): 
#             Series or protocol name to be set in the DICOM metadata.
#         value_range (tuple): 
#             Tuple (min, max) specifying the value range for the VOI LUT (windowing).
#         content_label (str): 
#             Label describing the quantitative map type, e.g., 'CBF', 'CVR', 'AAT', or 'ATA'.

#     Raises:
#         ValueError: If content_label is not provided or invalid.

#     Notes:
#         - The image is scaled into the 16-bit integer range for compatibility with DICOM.
#         - Metadata such as RescaleSlope, SeriesDescription, VOI LUT, and labeling details are inserted.
#         - All DICOM manipulation is done with `pydicom`.
#         - Orientation and layout are matched empirically to Philips DICOM expectations.
#         - Output format is decided based on whether the template DICOM is multiframe.

    if not content_label:
        raise ValueError("Please supply a content_label such as 'CBF', 'CVR', 'AAT', or 'ATA'")

    use_signed = content_label.upper() == 'CVR'
    unit_str = {
        'CBF': 'ml/100g/min',
        'CVR': 'ml/100g/min',
        'ATA': 'ml/100g/min',
        'AAT': 's'
    }.get(content_label.upper(), '')

    image = np.nan_to_num(image, nan=0.0)

    if use_signed:
        abs_max = np.max(np.abs(image))
        scalingfactor = (2**15 - 1) / abs_max
        image_scaled = (image * scalingfactor).clip(-2**15, 2**15 - 1).astype(np.int16)
    else:
        scalingfactor = (2**16 - 1) / np.nanmax(image)
        image_scaled = (image * scalingfactor).clip(0, 2**16 - 1).astype(np.uint16)

    ds = pydicom.dcmread(template_dicom_path)
    is_multiframe = hasattr(ds, 'PerFrameFunctionalGroupsSequence')
    is_singleframe = not is_multiframe

    if is_multiframe:
        logging.info("Detected Enhanced Multiframe DICOM template.")

        image_fordicom = np.flip(np.transpose(image_scaled, (2, 1, 0)), axis=1)
        ds.PixelData = image_fordicom.tobytes()

        ds.ImageType = ['DERIVED', 'SECONDARY', 'QUANTITATIVE']
        ds.ContentLabel = content_label.upper()
        ds.ContentDescription = f"ASL derived map: {content_label.upper()} [{unit_str}]"
        ds.SeriesDescription = f"{name} [{unit_str}]"
        ds.ProtocolName = name
        ds.RescaleSlope = 1.0 / scalingfactor
        ds.RescaleIntercept = 0.0
        ds.BitsAllocated = 16
        ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 1 if use_signed else 0
        ds.SOPClassUID = pydicom.uid.UID("1.2.840.10008.5.1.4.1.1.4.1")
        ds.SeriesInstanceUID = generate_uid()
        ds.NumberOfFrames = image.shape[2]

        smallest = int(np.min(image_scaled))
        largest = int(np.max(image_scaled))
        vr = 'SS' if use_signed else 'US'
        ds[Tag(0x0028, 0x0106)] = DataElement(Tag(0x0028, 0x0106), vr, smallest)
        ds[Tag(0x0028, 0x0107)] = DataElement(Tag(0x0028, 0x0107), vr, largest)

        ref_item = Dataset()
        ref_item.SeriesInstanceUID = ds.SeriesInstanceUID
        ref_item.ReferencedInstanceSequence = Sequence([])
        ds.ReferencedSeriesSequence = Sequence([ref_item])

        for frame in ds.PerFrameFunctionalGroupsSequence:
            if hasattr(frame, "PixelValueTransformationSequence"):
                frame.PixelValueTransformationSequence[0].RescaleSlope = 1.0 / scalingfactor
            if hasattr(frame, "FrameVOILUTSequence"):
                frame.FrameVOILUTSequence[0].WindowCenter = float(np.mean(value_range))
                frame.FrameVOILUTSequence[0].WindowWidth = float(np.ptp(value_range))

        output_filename = os.path.basename(template_dicom_path)
        output_path = os.path.join(output_dicom_dir, f"{output_filename}.dcm")
        ds.save_as(output_path)
        logging.info(f"Saved multiframe DICOM: {output_path}")

    elif is_singleframe: 
        # --- SINGLEFRAME EXPORT ---
        logging.info("Detected template directory. Saving as single-frame DICOM series.")
        # Extract directory and prefix
        template_dir = os.path.dirname(template_dicom_path)
        template_prefix = os.path.basename(template_dicom_path).rsplit('_', 1)[0] + '_'  # e.g. 'sWIP_ASL_CBF_preACZ_603_'

        template_files = [
            f for f in os.listdir(template_dir)
            if f.startswith(template_prefix) and os.path.isfile(os.path.join(template_dir, f))
        ]

        # Pair each file with its InstanceNumber
        file_instance_pairs = []
        for f in template_files:
            path = os.path.join(template_dir, f)
            try:
                ds = pydicom.dcmread(path, stop_before_pixels=True)
                instance_number = getattr(ds, 'InstanceNumber', None)
                if instance_number is not None:
                    file_instance_pairs.append((f, instance_number))
                else:
                    raise ValueError(f"Missing InstanceNumber in {f}")
            except Exception as e:
                raise RuntimeError(f"Error reading DICOM {f}: {e}")

        # Sort files by InstanceNumber, ie singeflarme dicoms will have correct anatomical order
        template_files_sorted = [f for f, _ in sorted(file_instance_pairs, key=lambda x: x[1])]

        if len(template_files_sorted) != image.shape[2]:
            raise ValueError("Mismatch between number of template DICOMs and image slices")
        
        series_instance_uid = generate_uid()

        for i, fname in enumerate(template_files_sorted):
            template_file = os.path.join(template_dir, fname)
            ds = pydicom.dcmread(template_file)
            ds.decompress()  # Ensure dataset is uncompressed before writing new PixelData
            ds.ImageType = ['DERIVED', 'SECONDARY', 'QUANTITATIVE']
            ds.ContentLabel = content_label.upper()
            ds.ContentDescription = f"ASL derived map: {content_label.upper()} [{unit_str}]"
            ds.SeriesDescription = f"{name} [{unit_str}]"
            ds.ProtocolName = name
            ds.RescaleSlope = 1.0 / scalingfactor
            ds.RescaleIntercept = 0.0
            ds.BitsAllocated = 16
            ds.BitsStored = 16
            ds.HighBit = 15
            ds.PixelRepresentation = 1 if use_signed else 0
            ds.SOPClassUID = pydicom.uid.UID("1.2.840.10008.5.1.4.1.1.4") # singleframe classic DICOM
            ds.SOPInstanceUID = generate_uid()
            ds.SeriesInstanceUID = series_instance_uid

            slice_img = np.flipud(image_scaled[:, :, i].T)
            ds.Rows, ds.Columns = slice_img.shape
            ds.PixelData = slice_img.tobytes()

            smallest = int(np.min(slice_img))
            largest = int(np.max(slice_img))

            output_path = os.path.join(output_dicom_dir, f"{fname}.dcm")
            ds.save_as(output_path)
            logging.info(f"Saved slice DICOM: {output_path}")
