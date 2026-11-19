"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

Shared DICOM helpers for the grayscale (save_data_dicom_grayscale) and PALETTE COLOR (save_data_dicom_color)
map exporters: scaling, UIDs, patient/study defaults, common metadata, source reference, template frame/file
selection, timing cleanup and series numbering.

Repository: https://github.com/JSIERO/ClinicalASL
Author: Jeroen Siero, UMCU (University Medical Center Utrecht), The Netherlands
License: BSD 3-Clause License
"""
import os
import logging
import datetime
import numpy as np
import pydicom
from pydicom.uid import generate_uid, ExplicitVRLittleEndian, UID
from pydicom.sequence import Sequence
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.dataelem import DataElement
from pydicom.tag import Tag
from clinical_asl_pipeline.__version__ import __version__ as TOOL_VERSION

IMPLEMENTATION_UID_ROOT = "1.3.6.1.4.1.54321.1.1"  # ClinicalASL UID root (not a registered PEN)
UID_PREFIX = IMPLEMENTATION_UID_ROOT + '.'
UNITS = {'CBF': 'ml/100g/min', 'CVR': 'ml/100g/min', 'ATA': 'ml/100g/min', 'AAT': 's'}


def unit_for(type_tag):
    if not type_tag:
        raise ValueError("Please supply a type_tag such as 'CBF', 'CVR', 'AAT', or 'ATA'")
    return UNITS.get(type_tag.upper(), '')


def scale_image(image, signed):
    # Scale a map to the full int16 (signed, CVR) or uint16 range; returns (scaled, scalingfactor),
    # real value = stored / scalingfactor. An all-zero map gets factor 1 instead of inf.
    image = np.nan_to_num(image, nan=0.0)
    peak = float(np.max(np.abs(image)) if signed else np.max(image))
    if not np.isfinite(peak) or peak <= 0:
        logging.warning("Map has no positive values: stored as zeros (scaling factor 1)")
        peak = 2**15 - 1 if signed else 2**16 - 1
    if signed:
        scalingfactor = (2**15 - 1) / peak
        return (image * scalingfactor).clip(-2**15, 2**15 - 1).astype(np.int16), scalingfactor
    scalingfactor = (2**16 - 1) / peak
    return (image * scalingfactor).clip(0, 2**16 - 1).astype(np.uint16), scalingfactor


def series_number(ds, series_number_incr, type_tag, cvr_number):
    # Map series: source SeriesNumber + 10 + increment; CVR has a fixed number (888 grayscale, 988 colour)
    if type_tag.upper() == 'CVR':
        return cvr_number
    return int(getattr(ds, 'SeriesNumber', 0)) + 10 + series_number_incr


def set_instance_uids(ds, sop_class_uid):
    # New SOP instance (file meta + dataset) of the given class; patient/study tags kept from the template,
    # with defaults only where the template lacks them.
    uid = generate_uid(prefix=UID_PREFIX)
    ds.file_meta = ds.file_meta or FileMetaDataset()
    ds.file_meta.MediaStorageSOPInstanceUID = uid
    ds.file_meta.MediaStorageSOPClassUID = UID(sop_class_uid)
    ds.file_meta.FileMetaInformationVersion = b'\x00\x01'
    ds.file_meta.ImplementationClassUID = IMPLEMENTATION_UID_ROOT
    ds.file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    ds.SOPInstanceUID = uid
    ds.SOPClassUID = UID(sop_class_uid)
    ds.StudyInstanceUID = getattr(ds, 'StudyInstanceUID', generate_uid(prefix=UID_PREFIX))
    ds.StudyID = getattr(ds, 'StudyID', 'Unknown')
    ds.PatientName = getattr(ds, 'PatientName', 'Anonymous^Patient')
    ds.PatientID = getattr(ds, 'PatientID', '000000')
    ds.PatientBirthDate = getattr(ds, 'PatientBirthDate', '')
    ds.PatientSex = getattr(ds, 'PatientSex', '')
    ds.ReferringPhysicianName = getattr(ds, 'ReferringPhysicianName', 'Unknown^Referring Physician')
    ds.StudyDate = getattr(ds, 'StudyDate', '')
    ds.StudyTime = getattr(ds, 'StudyTime', '')
    ds.AccessionNumber = getattr(ds, 'AccessionNumber', '')
    ds.ImageType = ['DERIVED', 'SECONDARY', 'QUANTITATIVE']


def set_common_metadata(ds, name, unit_str, type_tag, conversion_type=None):
    now = datetime.datetime.now()
    ds.ContentDate = now.strftime('%Y%m%d')
    ds.ContentTime = now.strftime('%H%M%S.%f')
    ds.ContentLabel = type_tag.upper()
    ds.ContentDescription = f"ClinicalASL-Siero: {type_tag.upper()} [{unit_str}]"
    ds.SeriesDescription = f"{name} (ClinicalASL, NOT CE MARKED)"
    ds.DerivationDescription = f"NOT CE MARKED - ClinicalASL-Siero: {type_tag.upper()}"

    # DCM 126302: post-processed derivative (PS3.16 Table D-1), so PACS identifies these as derived images
    code_item = Dataset()
    code_item.CodeValue = "126302"
    code_item.CodingSchemeDesignator = "DCM"
    code_item.CodeMeaning = "Post-processing"
    ds.DerivationCodeSequence = Sequence([code_item])

    ds.SoftwareVersions = f'ClinicalASL v{TOOL_VERSION}, https://github.com/JSIERO/ClinicalASL'
    ds.InstitutionName = getattr(ds, 'InstitutionName', 'University Medical Center Utrecht')
    ds.Manufacturer = f"ClinicalASL v{TOOL_VERSION}, https://github.com/JSIERO/ClinicalASL"
    ds.ImageComments = "ClinicalASL, NOT CE MARKED"
    if conversion_type:
        ds.ConversionType = conversion_type  # SC Equipment module, Type 1
    for tag in ('SeriesDescription', 'Manufacturer', 'SoftwareVersions', 'InstitutionName'):
        setattr(ds, tag, getattr(ds, tag)[:64])  # VR LO max length
    ds.ProtocolName = name


def add_source_dicom_reference(ds, source_dicom_path):
    # ReferencedSeriesSequence (+ ReferencedSOPInstance/ClassUID) pointing at the source ASL DICOM
    try:
        ref = pydicom.dcmread(source_dicom_path, stop_before_pixels=True)
        series_uid = getattr(ref, 'SeriesInstanceUID', None)
        instance_uid = getattr(ref, 'SOPInstanceUID', None)
        class_uid = getattr(ref, 'SOPClassUID', None)
        if not (series_uid and instance_uid and class_uid):
            logging.warning("Source DICOM missing UID fields, skipping ReferencedSeriesSequence.")
            return
        instance = Dataset()
        instance.ReferencedSOPClassUID = class_uid
        instance.ReferencedSOPInstanceUID = instance_uid
        series = Dataset()
        series.SeriesInstanceUID = series_uid
        series.ReferencedInstanceSequence = Sequence([instance])
        ds.ReferencedSeriesSequence = Sequence([series])
        ds.ReferencedSOPInstanceUID = instance_uid
        ds.ReferencedSOPClassUID = class_uid
    except Exception as e:
        logging.warning(f"Failed to add ReferencedSeriesSequence: {e}")


def select_multiframe_frames(ds):
    # Keep the frames of PLD 1, dynamic 1, condition 1 (one per slice) of a Philips Enhanced MR template;
    # set NumberOfFrames, PixelSpacing and SliceThickness. Returns the number of slices.
    try:
        private_seq = ds.PerFrameFunctionalGroupsSequence[0].get((0x2005, 0x140f))
        nplds = int(ds.get((0x2001, 0x1017)).value)
        ndyns = int(private_seq[0].get((0x0020, 0x0105)).value)
        nslices = int(ds.get((0x2001, 0x1018)).value)
        nconditions = int(ds.NumberOfFrames) // (nplds * ndyns * nslices)
        logging.info(f"Template structure: {nslices} slices x {nplds} PLDs x {ndyns} dynamics x {nconditions} conditions")
    except Exception as e:
        raise ValueError(f"Failed to extract dimensions from private tags: {e}")

    frames = ds.PerFrameFunctionalGroupsSequence
    ds.PerFrameFunctionalGroupsSequence = Sequence([frames[s * ndyns * nplds] for s in range(nslices)])
    ds.NumberOfFrames = nslices

    measures = ds.PerFrameFunctionalGroupsSequence[0]
    try:
        spacing = measures.PixelMeasuresSequence[0].PixelSpacing
        ds.PixelSpacing = [float(spacing[0]), float(spacing[1])]
    except Exception as e:
        logging.warning(f"Failed to extract PixelSpacing: {e}")
        ds.PixelSpacing = [1.0, 1.0]
    try:
        ds.SliceThickness = float(measures.PixelMeasuresSequence[0].SliceThickness)
    except Exception:
        ds.SliceThickness = 1.0
    return nslices


def multiframe_pixels(image_scaled):
    # (x, y, z) map -> Enhanced MR frame order (z, rows, cols), Philips orientation
    return np.flip(np.transpose(image_scaled, (2, 1, 0)), axis=1).tobytes()


def singleframe_pixels(slice_2d):
    # (x, y) slice -> (rows, cols), Philips orientation
    return np.flipud(slice_2d.T)


def clean_frame_timing(frame):
    # Per-frame functional group items that describe the source acquisition, not the derived map
    if 'CardiacSynchronizationSequence' in frame:
        for item in frame.CardiacSynchronizationSequence:
            if 'TriggerTime' in item:
                del item.TriggerTime
    if 'MRVelocityEncodingSequence' in frame:
        for item in frame.MRVelocityEncodingSequence:
            if 'VelocityEncodingDirection' in item:
                delattr(item, 'VelocityEncodingDirection')
    if 'RealWorldValueMappingSequence' in frame:
        del frame.RealWorldValueMappingSequence  # inherited Philips RWVM refers to the source data


def select_singleframe_templates(template_dicom_path, num_slices):
    # Source files of the same series (same filename prefix) with TemporalPositionIdentifier 1 and Philips phase
    # (2001,1008) 1, one per unique ImagePositionPatient z, sorted by z; at most num_slices. Returns full paths.
    template_dir = os.path.dirname(template_dicom_path)
    prefix = os.path.basename(template_dicom_path).rsplit('_', 1)[0] + '_'
    by_z = {}
    for f in sorted(os.listdir(template_dir)):
        path = os.path.join(template_dir, f)
        if not f.startswith(prefix) or not os.path.isfile(path):
            continue
        try:
            ds = pydicom.dcmread(path, stop_before_pixels=True, force=True)
            if getattr(ds, 'TemporalPositionIdentifier', None) == 1 and ds.get((0x2001, 0x1008)).value == 1:
                by_z.setdefault(ds.ImagePositionPatient[2], path)
        except Exception as e:
            logging.error(f"Error reading template DICOM {f}: {e}")
    if not by_z:
        raise FileNotFoundError(f"No single-frame template DICOMs (prefix {prefix}) in {template_dir}")
    return [by_z[z] for z in sorted(by_z)][:num_slices]


def clean_singleframe_timing(ds):
    # Remove/neutralise source timing so PACS sorts the derived slices as one volume
    if 'TriggerTime' in ds:
        del ds.TriggerTime
    private_tag, velocity_tag = Tag(0x2005, 0x140f), Tag(0x0018, 0x9090)  # Philips per-frame seq, VelocityEncodingDirection
    if private_tag in ds:
        for item in ds[private_tag].value:
            if velocity_tag in item:
                del item[velocity_tag]
    if 'NumberOfTemporalPositions' in ds:
        ds.NumberOfTemporalPositions = 1
    if Tag(0x2001, 0x1081) in ds:  # Philips number of dynamics
        ds[Tag(0x2001, 0x1081)] = DataElement(Tag(0x2001, 0x1081), 'IS', 1)
    if 'TemporalPositionIdentifier' in ds:
        ds.TemporalPositionIdentifier = 1


def set_16bit_pixels(ds, signed):
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 1 if signed else 0
    ds.SamplesPerPixel = 1


def output_path(output_dicom_dir, name, ds, singleframe):
    base = f"{name.replace(' ', '_')}_{ds.SeriesNumber}"
    return os.path.join(output_dicom_dir, f"{base}_{ds.InstanceNumber}.dcm" if singleframe else f"{base}.dcm")


def log_series_saved(kind, name, ds, n_files):
    # One line per series; no patient or study identifiers in the log
    logging.info(f"Saved {kind} DICOM series {ds.SeriesNumber} '{name}': {n_files} file(s)")
