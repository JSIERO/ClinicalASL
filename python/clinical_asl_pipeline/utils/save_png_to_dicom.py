"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

DICOM to NIfTI conversion module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
This module provides functionality to convert PNG images to DICOM format.   
License: BSD 3-Clause License
"""

import pydicom
import numpy as np
import datetime
import logging
from pydicom.dataset import Dataset, FileDataset
from pydicom.uid import generate_uid, ExplicitVRLittleEndian, SecondaryCaptureImageStorage
from pydicom.uid import ExplicitVRLittleEndian

from PIL import Image

def save_png_to_dicom(png_path, output_path, series_description, instancenumber=1, template_dcm_path=None):
    #   Convert a PNG image to DICOM format.
    #   This function reads a PNG image, converts it to a DICOM dataset, and saves it to the specified output path.
    #   Parameters:
    #      png_path (str): Path to the input PNG image. 
    #      output_path (str): Path where the output DICOM file will be saved.   
    #      series_description (str): Title for the DICOM Series Description. 
    #      instancenumber (int): Instance number for the DICOM image. Default is 1.
    #      template_dcm_path (str): Optional path to a template DICOM file for metadata. If not provided, defaults will be used.
    #   #   Returns:
    #      None: The function saves the DICOM file to the specified output path.        

    # Load PNG image
    img = Image.open(png_path)
    img = img.convert('L' if img.mode == 'L' else 'RGB')  # Grayscale or RGB
    pixel_array = np.asarray(img)

    # Load template DICOM if provided
    template_ds = None
    if template_dcm_path:
        template_ds = pydicom.dcmread(template_dcm_path, stop_before_pixels=True)

    # File Meta
    file_meta = pydicom.Dataset()
    file_meta.MediaStorageSOPClassUID = SecondaryCaptureImageStorage
    file_meta.MediaStorageSOPInstanceUID = generate_uid()
    file_meta.ImplementationClassUID = generate_uid()
    file_meta.TransferSyntaxUID = ExplicitVRLittleEndian

    # Create new DICOM dataset
    ds = FileDataset(output_path, {}, file_meta=file_meta, preamble=b"\0" * 128)

    # Copy metadata from template or generate defaults
    ds.Modality = template_ds.Modality if template_ds and 'Modality' in template_ds else 'OT'
    ds.StudyInstanceUID = template_ds.StudyInstanceUID if template_ds and 'StudyInstanceUID' in template_ds else generate_uid()
    ds.SeriesInstanceUID = generate_uid()  # new series for SC images
    ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
    ds.SOPClassUID = file_meta.MediaStorageSOPClassUID
    ds.PatientName = template_ds.PatientName if template_ds and 'PatientName' in template_ds else "Anonymous^Patient"
    ds.PatientID = template_ds.PatientID if template_ds and 'PatientID' in template_ds else "000000"
    ds.PatientBirthDate = template_ds.PatientBirthDate
    ds.PatientSex = template_ds.PatientSex
    ds.StudyDate = template_ds.StudyDate
    ds.StudyTime = template_ds.StudyTime
    ds.AccessionNumber = template_ds.AccessionNumber
    ds.SeriesDescription = f"{series_description} - derived (ClinicalASL)"
    ds.InstanceNumber = instancenumber if instancenumber else 1 # Default to 1 if not provided
    ds.ImageType = "DERIVED\\SECONDARY\\QUANTITATIVE"
    ds.StudyID = template_ds.StudyID
    ds.SeriesNumber = template_ds.SeriesNumber + 100 + instancenumber  # Increment series number for new series
    ds.StudyDescription = template_ds.StudyDescription
    ds.SoftwareVersion = 'ClinicalASL https://github.com/JSIERO/ClinicalASL'
    
    # Add content date/time for this image
    now = datetime.datetime.now()
    ds.ContentDate = now.strftime('%Y%m%d')
    ds.ContentTime = now.strftime('%H%M%S.%f')

    # Image size and type
    ds.Rows, ds.Columns = pixel_array.shape[:2]
    if len(pixel_array.shape) == 2:  # Grayscale
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
    else:  # RGB
        ds.SamplesPerPixel = 3
        ds.PhotometricInterpretation = "RGB"
        ds.PlanarConfiguration = 0

    ds.BitsStored = 8
    ds.BitsAllocated = 8
    ds.HighBit = 7
    ds.PixelRepresentation = 0

    # Add pixel data
    ds.PixelData = pixel_array.tobytes()

    # Encoding
    ds.is_little_endian = True
    ds.is_implicit_VR = False

    # Save
    ds.save_as(output_path)
    logging.info(f"PNG-derived DICOM saved to {output_path}")