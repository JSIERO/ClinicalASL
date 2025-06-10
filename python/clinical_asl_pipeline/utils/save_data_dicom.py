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
import subprocess

def save_data_dicom(image, template_dicom_path, output_dicom_path, name, value_range, content_label):    #
    # Save a 3D ASL-derived image as a multiframe DICOM file using a Philips reference DICOM and DCMTK tools.
    # This function prepares and scales a quantitative ASL image (e.g., CBF, CVR, AAT, ATA), injects it into a DICOM
    # Enhanced MR Image template, and updates relevant DICOM metadata fields for quantitative analysis and visualization.
    #
    # Parameters:
    #     image (np.ndarray): 
    #         3D ASL image array with shape (Height, Width, Slices), typically loaded from a NIfTI file.
    #     template_dicom_path (str): 
    #         Path to a reference DICOM file to use as a template for metadata and structure.
    #     output_dicom_path (str): 
    #         Path where the output DICOM file will be saved.
    #     name (str): 
    #         Series or protocol name to be set in the DICOM metadata.
    #     value_range (tuple): 
    #         Tuple (min, max) specifying the value range for the VOI LUT (windowing) in the DICOM.
    #     content_label (str): 
    #         Label describing the quantitative map type, e.g., 'CBF', 'CVR', 'AAT', or 'ATA'.
    #
    # Raises:
    #     ValueError: If content_label is not provided.
    #
    # Notes:
    #     - The function uses DCMTK command-line tools (`dcmodify`) to update robust DICOM tags and inject pixel data.
    #     - Pixel data is scaled to fit into 16-bit signed or unsigned integer range, depending on the map type.
    #     - Metadata fields such as SeriesDescription, ProtocolName, RescaleSlope, and VOI LUT are updated for quantitative interpretation.
    #     - Fragile or complex DICOM fields are updated using pydicom after DCMTK modifications.
    #     - The function expects the input image orientation and shape to match the reference DICOM after empirical adjustment.
    #     - Temporary files are created for pixel data injection and are cleaned up after use.

    if not content_label:
        raise ValueError("Please supply a content_label such as 'CBF', 'CVR', 'AAT', or 'ATA'")

    use_signed = content_label.upper() == 'CVR'
    unit_str = {
        'CBF': 'ml/100g/min',
        'CVR': 'ml/100g/min',
        'ATA': 'ml/100g/min',
        'AAT': 's'
    }.get(content_label.upper(), '')

    # Replace NaNs with 0 before scaling
    image = np.nan_to_num(image, nan=0.0)
    
    # Scale image to integer range    
    if use_signed:
        abs_max = np.max(np.abs(image))
        scalingfactor = (2**15 - 1) / abs_max
        image_scaled = (image * scalingfactor).clip(-2**15, 2**15 - 1).astype(np.int16)
    else:
        scalingfactor = (2**16 - 1) / np.nanmax(image)
        image_scaled = (image * scalingfactor).clip(0, 2**16 - 1).astype(np.uint16)

    # Prepare image for DICOM insert:
    # DICOM Enhanced MR expects (Frames, Rows, Columns) layout
    
    # Transpose + flip to match Philips DICOM pixel ordering - verified empirically against reference DICOM
    image_4d = np.flip(np.transpose(image_scaled, (2, 1, 0)), axis=1)   #(Slices, Columns flipped, Rows)

    # Save pixel data to raw file
    raw_path = os.path.join(os.getcwd(), "temp_pixeldata.raw")
    image_4d.tofile(raw_path)

    # Copy template to output
    subprocess.run(["cp", template_dicom_path, output_dicom_path], check=True)

    # Modify robust tags via dcmodify
    dcmodify_cmd = [
        "dcmodify", "-nb", "-g",
        "-i", f"(0008,103e)={name} [{unit_str}]",
        "-i", f"(0018,1030)={name}",
        "-i", f"(0028,1053)={1.0 / scalingfactor:.6f}",
        "-i", "(0028,1052)=0.000000",
        "-i", "(0028,0100)=16",
        "-i", "(0028,0101)=16",
        "-i", "(0028,0102)=15",
        "-i", f"(0028,0103)={'1' if use_signed else '0'}",
        output_dicom_path
    ]
    subprocess.run(dcmodify_cmd, check=True)

    # Inject PixelData
    insert_cmd = ["dcmodify", "-nb", "--insert-from-file", f"(7fe0,0010)={raw_path}", output_dicom_path]
    subprocess.run(insert_cmd, check=True)
    # Remove the temporary raw file
    os.remove(raw_path)

    # Modify fragile fields using pydicom
    ds = pydicom.dcmread(output_dicom_path)
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
    ds.SOPClassUID = pydicom.uid.UID("1.2.840.10008.5.1.4.1.1.4.1")  # Enhanced MR Image Storage
    ds.SeriesInstanceUID = generate_uid()
    ds.NumberOfFrames = image.shape[2]

    # -- Set explicit VR for ambiguous value representation
    smallest = int(np.min(image_scaled))
    largest = int(np.max(image_scaled))
    vr = 'SS' if use_signed else 'US'
    ds[Tag(0x0028, 0x0106)] = DataElement(Tag(0x0028, 0x0106), vr, smallest)
    ds[Tag(0x0028, 0x0107)] = DataElement(Tag(0x0028, 0x0107), vr, largest)

    # -- Add ReferencedSeriesSequence with empty ReferencedInstanceSequence
    ref_item = Dataset()
    ref_item.SeriesInstanceUID = ds.SeriesInstanceUID
    ref_item.ReferencedInstanceSequence = Sequence([])  # empty sequence
    ds.ReferencedSeriesSequence = Sequence([ref_item])

    # -- Update PerFrameFunctionalGroupsSequence
    if hasattr(ds, "PerFrameFunctionalGroupsSequence"):
        for frame in ds.PerFrameFunctionalGroupsSequence:
            if hasattr(frame, "PixelValueTransformationSequence"):
                frame.PixelValueTransformationSequence[0].RescaleSlope = 1.0 / scalingfactor
            if hasattr(frame, "Private_2005_140f"):
                frame.Private_2005_140f[0].RescaleSlope = 1.0 / scalingfactor
            if hasattr(frame, "FrameVOILUTSequence"):
                frame.FrameVOILUTSequence[0].WindowCenter = float(np.mean(value_range))
                frame.FrameVOILUTSequence[0].WindowWidth = float(np.ptp(value_range))

    ds.save_as(output_dicom_path)

    logging.info(f"DICOM written to: {output_dicom_path}")
    logging.info(f"  Content: {content_label.upper()} ({unit_str}), Slope: {1 / scalingfactor:.6f}")
