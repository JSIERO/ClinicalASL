import os
import numpy as np
import pydicom
from pydicom.uid import generate_uid
from pydicom.sequence import Sequence
from pydicom.dataset import Dataset
from pydicom.dataelem import DataElement
from pydicom.tag import Tag
import subprocess

def save_data_dicom(image, template_dicom_path, output_dicom_path, name, value_range, content_label):
    """
    Save ASL-derived image as multiframe DICOM using Philips reference and DCMTK.

    Parameters:
        image (np.ndarray): 3D ASL image (H x W x S), loaded from NIfTI
        template_dicom_path (str): path to reference DICOM file
        output_dicom_path (str): path to output DICOM
        name (str): series/protocol name
        value_range (tuple): (min, max) for VOI LUT
        content_label (str): 'CBF', 'CVR', 'AAT', 'ATA'
    """
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

    print(f"DICOM written to: {output_dicom_path}")
    print(f"  Content: {content_label.upper()} ({unit_str}), Slope: {1 / scalingfactor:.6f}")
