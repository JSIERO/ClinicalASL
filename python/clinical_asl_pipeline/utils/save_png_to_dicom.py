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
import os
from clinical_asl_pipeline.__version__ import __version__ as TOOL_VERSION
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import generate_uid, ExplicitVRLittleEndian, SecondaryCaptureImageStorage
from PIL import Image

def save_png_to_dicom(png_path, output_path, series_description, series_instance_uid,  instance_number=1, source_dicom_path=None):
    #   Convert a PNG image to DICOM format.
    #   This function reads a PNG image, converts it to a DICOM dataset, and saves it to the specified output path.
    #   Parameters:
    #      png_path (str): Path to the input PNG image. 
    #      output_path (str): Path where the output DICOM file will be saved.   
    #      series_description (str): Title for the DICOM Series Description. 
    #      SeriesInstanceUID (str): Unique identifier for the DICOM Series. This is used to group related images.
    #      instancenumber (int): Instance number for the DICOM image. Default is 1.
    #      source_dicom_path (str):  path to a source DICOM file for referencing and template for metadata. 
    #   #   Returns:
    #      None: The function saves the DICOM file to the specified output path. 

    now = datetime.datetime.now()
    IMPLEMENTATION_UID_ROOT = "1.3.6.1.4.1.54321.1.1" # Example root UID for ClinicalASL, fake PEN

    # Load PNG image
    img = Image.open(png_path)
    img = img.convert('L' if img.mode == 'L' else 'RGB')  # Grayscale or RGB
    pixel_array = np.asarray(img)

    # Load template DICOM if provided
    template_ds = None
    if source_dicom_path:
        template_ds = pydicom.dcmread(source_dicom_path, stop_before_pixels=True, force=True)

    # File Meta
    # Generate consistent UID
    sop_instance_uid  = generate_uid(prefix=IMPLEMENTATION_UID_ROOT + '.')
    file_meta = FileMetaDataset()
    file_meta.FileMetaInformationVersion = b'\x00\x01'
    file_meta.MediaStorageSOPClassUID = SecondaryCaptureImageStorage
    file_meta.MediaStorageSOPInstanceUID = sop_instance_uid 
    file_meta.ImplementationClassUID = IMPLEMENTATION_UID_ROOT
    file_meta.TransferSyntaxUID = ExplicitVRLittleEndian

    # Create new DICOM dataset
    ds = FileDataset(output_path, {}, file_meta=file_meta, preamble=b"\0" * 128)
    ds.file_meta = file_meta
    ds.SOPInstanceUID = sop_instance_uid 
    ds.SOPClassUID = file_meta.MediaStorageSOPClassUID

    # Copy metadata from template or generate defaults
    ds.Modality = getattr(template_ds, 'Modality', 'OT')
    ds.ImageType = "DERIVED\\SECONDARY\\QUANTITATIVE"

    ds.StudyInstanceUID = getattr(template_ds, 'StudyInstanceUID', None)
    ds.SeriesInstanceUID = series_instance_uid  # use provided SeriesInstanceUID
    ds.InstanceCreatorUID = IMPLEMENTATION_UID_ROOT
    ds.SeriesNumber = 999  # Secondary Capture images typically have a fixed series number, here 999 is used as a placeholder
    ds.InstanceNumber = instance_number
    ds.SeriesDescription = f"{series_description} - (ClinicalASL-Siero)"
    
    # 2) Copy StudyID exactly as your site uses it; enforce SH (<=12 ASCII)
    raw_sid = getattr(template_ds, 'StudyID', None)
    if raw_sid is None or str(raw_sid).strip() == '':
        # last-resort: derive something stable but short; better: pass study_id explicitly
        raw_sid = getattr(template_ds, 'AccessionNumber', '') or 'STUDY'
        
    # Keep ASCII-only and truncate to 12
    sid_ascii = ''.join(ch for ch in str(raw_sid) if ord(ch) < 128)[:12]
    ds.StudyID = sid_ascii

    ds.PerformedProcedureStepID = getattr(template_ds, 'PerformedProcedureStepID', 'UnknownProcedure')
    ds.RequestedProcedureID = getattr(template_ds, 'RequestedProcedureID', 'UnknownRequest')
    ds.StudyDescription = getattr(template_ds, 'StudyDescription', 'Unknown Study')
    ds.PatientName = getattr(template_ds, 'PatientName', 'Anonymous^Patient')
    ds.PatientID = getattr(template_ds, 'PatientID', '000000')
    ds.PatientBirthDate = getattr(template_ds, 'PatientBirthDate', '')
    ds.PatientAge = getattr(template_ds, 'PatientAge', '')
    ds.PatientWeight = getattr(template_ds, 'PatientWeight', '')
    ds.PatientSex = getattr(template_ds, 'PatientSex', '')       
    ds.StudyDate = getattr(template_ds, 'StudyDate', now.strftime('%Y%m%d'))
    ds.StudyTime = getattr(template_ds, 'StudyTime', now.strftime('%H%M%S'))
    ds.AccessionNumber = getattr(template_ds, 'AccessionNumber', '') 
    ds.BodyPartExamined = getattr(template_ds, 'BodyPartExamined', 'BRAIN')
    ds.PatientOrientation = ['L', 'P']
    ds.SoftwareVersions = f'ClinicalASL v{TOOL_VERSION}, https://github.com/JSIERO/ClinicalASL'
    ds.InstitutionName = getattr(template_ds, 'InstitutionName', 'University Medical Center Utrecht')
    ds.Manufacturer = f"ClinicalASL v{TOOL_VERSION}, https://github.com/JSIERO/ClinicalASL"
    ds.ReferringPhysicianName = getattr(template_ds, 'ReferringPhysicianName', 'Unknown^Referring Physician')  
    ds.ConversionType = "WSD"  # WSD = Workstation, standard value for derived images
    ds.SeriesDescription = ds.SeriesDescription[:64]  # VR LO max length
    ds.Manufacturer = ds.Manufacturer[:64]  # VR LO max length
    ds.SoftwareVersions = ds.SoftwareVersions[:64]  # VR LO max length
    ds.InstitutionName = ds.InstitutionName[:64]  # VR LO max length
    # Add content date/time for this image
    ds.ContentDate = now.strftime('%Y%m%d')
    ds.ContentTime = now.strftime('%H%M%S.%f')
    ds.InstanceCreationDate = now.strftime('%Y%m%d')
    ds.InstanceCreationTime = now.strftime('%H%M%S')

    # Pixel data: image size and type
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
    
    # Save
    # Modify output_path to append 'PNG' and SeriesNumber before the file extension
    base, ext = os.path.splitext(output_path)
    new_output_path = f"{base}_{ds.SeriesNumber}_PNG{ext}"
    ds.save_as(new_output_path, enforce_file_format=True)

    # logging of key metadata
    logging.info(f"StudyID PNG-derived DICOM: {ds.StudyID}")
    logging.info(f"StudyInstanceUID PNG-derived DICOM: {ds.StudyInstanceUID}")
    logging.info(f"PatientID PNG-derived DICOM: {ds.PatientID}")
    logging.info(f"StudyDate PNG-derived DICOM: {ds.StudyDate}")
    logging.info(f"StudyTime PNG-derived DICOM: {ds.StudyTime}")

