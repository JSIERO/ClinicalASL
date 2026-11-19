"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

DICOM saving utility module (MONOCHROME2 - quantitative grayscale output).

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Saves a quantitative ASL map (CBF, AAT, ATA, CVR) as grayscale DICOM using a source ASL DICOM as template:
    Enhanced MR multi-frame for a multi-frame template, MR Image Storage per slice for a single-frame (PACS) template.
    Real values = stored x RescaleSlope (RealWorldValueMapping set from the map's own scaling). Shared helpers in
    dicom_common.py.

License: BSD 3-Clause License
"""
import numpy as np
import pydicom
from pydicom.uid import generate_uid
from pydicom.sequence import Sequence
from pydicom.dataset import Dataset
from pydicom.dataelem import DataElement
from pydicom.tag import Tag
from clinical_asl_pipeline.utils import dicom_common as dc

ENHANCED_MR = "1.2.840.10008.5.1.4.1.1.4.1"  # Enhanced MR Image Storage (multi-frame)
MR_IMAGE = "1.2.840.10008.5.1.4.1.1.4"  # MR Image Storage (single-frame)
CVR_SERIES_NUMBER = 888


def set_real_world_value_mapping(ds, type_tag, unit_str):
    """Real World Value mapping for the derived map (replaces the inherited Philips RWVM of the source data).
    Uses RescaleSlope/Intercept already set on ds."""
    mu = Dataset()
    mu.CodeValue = unit_str
    mu.CodingSchemeDesignator = "UCUM"
    mu.CodeMeaning = unit_str

    rwv = Dataset()
    if ds.PixelRepresentation == 1:  # signed (CVR): [-32768, 32767]
        rwv.RealWorldValueFirstValueMapped = -(2 ** (ds.BitsStored - 1))
        rwv.RealWorldValueLastValueMapped = 2 ** (ds.BitsStored - 1) - 1
    else:  # unsigned (CBF, AAT, ATA): [0, 65535]
        rwv.RealWorldValueFirstValueMapped = 0
        rwv.RealWorldValueLastValueMapped = 2 ** ds.BitsStored - 1
    rwv.RealWorldValueSlope = float(ds.RescaleSlope)
    rwv.RealWorldValueIntercept = float(ds.RescaleIntercept)
    rwv.MeasurementUnitsCodeSequence = Sequence([mu])
    rwv.LUTExplanation = f"{type_tag.upper()} [{unit_str}]"
    rwv.LUTLabel = type_tag.upper()  # Type 1 in the RWV Mapping Item Macro (SH, max 16 chars)
    ds.RealWorldValueMappingSequence = Sequence([rwv])


def save_data_dicom(image, source_dicom_path, output_dicom_dir, name, value_range, type_tag, series_number_incr):
    unit_str = dc.unit_for(type_tag)
    signed = type_tag.upper() == 'CVR'
    image_scaled, scalingfactor = dc.scale_image(image, signed)
    slope = f"{1.0 / scalingfactor:.10g}"
    center, width = f"{np.mean(value_range):.10g}", f"{np.ptp(value_range):.10g}"

    ds = pydicom.dcmread(source_dicom_path, force=True)

    if hasattr(ds, 'PerFrameFunctionalGroupsSequence'):  # --- MULTIFRAME (Enhanced MR) ---
        dc.select_multiframe_frames(ds)
        ds.PixelData = dc.multiframe_pixels(image_scaled)
        ds.SeriesInstanceUID = generate_uid(prefix=dc.UID_PREFIX)
        dc.set_instance_uids(ds, ENHANCED_MR)
        ds.SeriesNumber = dc.series_number(ds, series_number_incr, type_tag, CVR_SERIES_NUMBER)
        dc.set_16bit_pixels(ds, signed)
        ds.PhotometricInterpretation = "MONOCHROME2"
        vr = 'SS' if signed else 'US'
        ds[Tag(0x0028, 0x0106)] = DataElement(Tag(0x0028, 0x0106), vr, int(np.min(image_scaled)))  # Smallest
        ds[Tag(0x0028, 0x0107)] = DataElement(Tag(0x0028, 0x0107), vr, int(np.max(image_scaled)))  # Largest
        dc.set_common_metadata(ds, name, unit_str, type_tag)
        dc.add_source_dicom_reference(ds, source_dicom_path)
        for frame in ds.PerFrameFunctionalGroupsSequence:
            if hasattr(frame, "PixelValueTransformationSequence"):
                frame.PixelValueTransformationSequence[0].RescaleSlope = slope
                frame.PixelValueTransformationSequence[0].RescaleIntercept = f"{0.0:.10g}"
            if hasattr(frame, "FrameVOILUTSequence"):
                frame.FrameVOILUTSequence[0].WindowCenter = center
                frame.FrameVOILUTSequence[0].WindowWidth = width
            dc.clean_frame_timing(frame)
        set_real_world_value_mapping(ds, type_tag, unit_str)
        ds.save_as(dc.output_path(output_dicom_dir, name, ds, singleframe=False), enforce_file_format=True)
        dc.log_series_saved('multi-frame grayscale', name, ds, 1)
        return

    # --- SINGLEFRAME (MR Image Storage, one file per slice) ---
    templates = dc.select_singleframe_templates(source_dicom_path, image.shape[2])
    series_instance_uid = generate_uid(prefix=dc.UID_PREFIX)
    for i, template_file in enumerate(templates):
        ds = pydicom.dcmread(template_file, force=True)
        ds.decompress()
        ds.SeriesInstanceUID = series_instance_uid
        dc.set_instance_uids(ds, MR_IMAGE)
        ds.SeriesNumber = dc.series_number(ds, series_number_incr, type_tag, CVR_SERIES_NUMBER)
        dc.clean_singleframe_timing(ds)
        ds.InstanceNumber = i + 1
        ds.RescaleSlope = slope
        ds.RescaleIntercept = 0.0
        dc.set_16bit_pixels(ds, signed)
        ds.PixelSpacing = getattr(ds, 'PixelSpacing', [1.0, 1.0])
        ds.SliceThickness = getattr(ds, 'SliceThickness', 1.0)
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.WindowCenter = center
        ds.WindowWidth = width
        slice_img = dc.singleframe_pixels(image_scaled[:, :, i])
        ds.Rows, ds.Columns = slice_img.shape
        ds.PixelData = slice_img.tobytes()
        dc.set_common_metadata(ds, name, unit_str, type_tag)
        dc.add_source_dicom_reference(ds, source_dicom_path)
        set_real_world_value_mapping(ds, type_tag, unit_str)
        ds.save_as(dc.output_path(output_dicom_dir, name, ds, singleframe=True), enforce_file_format=True)
    dc.log_series_saved('single-frame grayscale', name, ds, len(templates))
