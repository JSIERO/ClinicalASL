"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

DICOM saving utility module (PALETTE COLOR output).

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Saves an ASL map (CBF, AAT, ATA, CVR) as PALETTE COLOR DICOM with an embedded 65536-entry colormap LUT (same
    colormap and range as the PNG output) for direct colour display in PACS viewers (e.g. Sectra IDS7, MicroDICOM).
    Visual only: no Rescale, W/L or RealWorldValueMapping (quantitative values are in the grayscale series).
    Multi-frame template -> Multi-frame True Color SC; single-frame (PACS) template -> Secondary Capture per slice.
    Shared helpers in dicom_common.py.

License: BSD 3-Clause License
"""
import logging
import numpy as np
import pydicom
import matplotlib.pyplot as plt
import importlib.resources as pkg_resources
from matplotlib.colors import ListedColormap
from scipy.io import loadmat
from pydicom.uid import generate_uid
from pydicom.dataelem import DataElement
from pydicom.tag import Tag
from clinical_asl_pipeline.utils import dicom_common as dc

MULTIFRAME_SC = "1.2.840.10008.5.1.4.1.1.7.4"  # Multi-frame True Color Secondary Capture Image Storage
SECONDARY_CAPTURE = "1.2.840.10008.5.1.4.1.1.7"  # Secondary Capture Image Storage
CVR_SERIES_NUMBER = 988  # grayscale CVR (888) + 100, like the other colour series

# Default colormap per map type (matches save_figure_to_png usage)
DEFAULT_COLORMAPS = {'CBF': 'viridis', 'CVR': 'vik', 'AAT': 'devon', 'ATA': 'viridis'}


def load_colormap(colormap_name):
    # matplotlib colormap by name, or the custom .mat colormaps shipped with ClinicalASL (vik, devon)
    name = colormap_name.lower()
    if name not in ('vik', 'devon'):
        return plt.get_cmap(name)
    try:
        with pkg_resources.files('clinical_asl_pipeline.colormaps').joinpath(f'{name}.mat').open('rb') as f:
            mat = loadmat(f)
    except FileNotFoundError:
        raise ValueError(f"Missing colormap file: {name}.mat")
    cmap_data = mat.get(name)
    if cmap_data is None:
        raise ValueError(f"Could not find colormap data in {name}.mat")
    return ListedColormap(cmap_data)


def generate_palette_lut(cmap, scalingfactor, value_range, rescale_intercept=0.0):
    # 65536-entry 16-bit LUT: stored value -> real value (stored / scalingfactor + intercept) -> normalised over
    # value_range (as the PNG) -> colormap. Index 0 is black (background / outside the mask).
    real_values = np.arange(65536, dtype=np.float64) / scalingfactor + rescale_intercept
    vmin, vmax = value_range
    normalized = np.full(65536, 0.5) if vmax == vmin else np.clip((real_values - vmin) / (vmax - vmin), 0, 1)
    colors = cmap(normalized)[:, :3]
    colors[0] = [0, 0, 0]
    lut16 = (colors * 65535).astype(np.uint16)
    return lut16[:, 0], lut16[:, 1], lut16[:, 2]


def set_palette_color_tags(ds, red, green, blue):
    ds.PhotometricInterpretation = "PALETTE COLOR"
    # Descriptor [0, 0, 16]: 0 encodes 65536 entries (PS3.3 C.7.6.3.1.6); VR US set explicitly ('US or SS' ambiguity)
    for tag in (0x1101, 0x1102, 0x1103):
        ds[Tag(0x0028, tag)] = DataElement(Tag(0x0028, tag), 'US', [0, 0, 16])
    ds.RedPaletteColorLookupTableData = red.tobytes()
    ds.GreenPaletteColorLookupTableData = green.tobytes()
    ds.BluePaletteColorLookupTableData = blue.tobytes()


def remove_display_scaling(ds):
    # The LUT defines the display mapping: no Rescale, W/L or RWVM on PALETTE COLOR Secondary Capture
    for tag_name in ('RescaleSlope', 'RescaleIntercept', 'RescaleType', 'WindowCenter', 'WindowWidth',
                     'RealWorldValueMappingSequence'):
        if tag_name in ds:
            delattr(ds, tag_name)


def save_data_dicom(image, source_dicom_path, output_dicom_dir, name, value_range, type_tag, series_number_incr,
                    colormap_name=None, mask=None):
    # image: 3D map (x, y, slices); mask: 3D, voxels that are NaN or 0 become pixel 0 (black). The mask is required
    # for signed data (CVR): after the unsigned offset, background would otherwise map to the colormap midpoint.
    unit_str = dc.unit_for(type_tag)
    colormap_name = colormap_name or DEFAULT_COLORMAPS.get(type_tag.upper(), 'viridis')
    signed = type_tag.upper() == 'CVR'
    image_scaled, scalingfactor = dc.scale_image(image, signed)
    rescale_intercept = 0.0
    if signed:  # PALETTE COLOR needs unsigned pixels: int16 [-32768, 32767] -> uint16 [0, 65535]
        image_scaled = (image_scaled.astype(np.int32) + 32768).astype(np.uint16)
        rescale_intercept = -32768.0 / scalingfactor
    if mask is not None:
        image_scaled[~(np.isfinite(mask) & (mask != 0))] = 0

    lut = generate_palette_lut(load_colormap(colormap_name), scalingfactor, value_range, rescale_intercept)
    logging.info(f"Applying PALETTE COLOR with colormap '{colormap_name}', range {value_range} for {type_tag}")

    ds = pydicom.dcmread(source_dicom_path, force=True)

    if hasattr(ds, 'PerFrameFunctionalGroupsSequence'):  # --- MULTIFRAME ---
        dc.select_multiframe_frames(ds)
        ds.PixelData = dc.multiframe_pixels(image_scaled)
        ds.SeriesInstanceUID = generate_uid(prefix=dc.UID_PREFIX)
        dc.set_instance_uids(ds, MULTIFRAME_SC)
        ds.SeriesNumber = dc.series_number(ds, series_number_incr, type_tag, CVR_SERIES_NUMBER)
        dc.set_16bit_pixels(ds, signed=False)
        set_palette_color_tags(ds, *lut)
        dc.set_common_metadata(ds, name, unit_str, type_tag, conversion_type="WSD")
        dc.add_source_dicom_reference(ds, source_dicom_path)
        for frame in ds.PerFrameFunctionalGroupsSequence:
            for seq in ('PixelValueTransformationSequence', 'FrameVOILUTSequence'):
                if seq in frame:
                    delattr(frame, seq)
            dc.clean_frame_timing(frame)
        remove_display_scaling(ds)
        ds.save_as(dc.output_path(output_dicom_dir, name, ds, singleframe=False), enforce_file_format=True)
        dc.log_series_saved('multi-frame PALETTE COLOR', name, ds, 1)
        return

    # --- SINGLEFRAME (Secondary Capture, one file per slice) ---
    templates = dc.select_singleframe_templates(source_dicom_path, image.shape[2])
    series_instance_uid = generate_uid(prefix=dc.UID_PREFIX)
    for i, template_file in enumerate(templates):
        ds = pydicom.dcmread(template_file, force=True)
        ds.decompress()
        ds.SeriesInstanceUID = series_instance_uid
        dc.set_instance_uids(ds, SECONDARY_CAPTURE)
        ds.SeriesNumber = dc.series_number(ds, series_number_incr, type_tag, CVR_SERIES_NUMBER)
        dc.clean_singleframe_timing(ds)
        ds.InstanceNumber = i + 1
        dc.set_16bit_pixels(ds, signed=False)
        ds.PixelSpacing = getattr(ds, 'PixelSpacing', [1.0, 1.0])
        ds.SliceThickness = getattr(ds, 'SliceThickness', 1.0)
        set_palette_color_tags(ds, *lut)
        remove_display_scaling(ds)
        slice_img = dc.singleframe_pixels(image_scaled[:, :, i])
        ds.Rows, ds.Columns = slice_img.shape
        ds.PixelData = slice_img.tobytes()
        dc.set_common_metadata(ds, name, unit_str, type_tag, conversion_type="WSD")
        dc.add_source_dicom_reference(ds, source_dicom_path)
        ds.save_as(dc.output_path(output_dicom_dir, name, ds, singleframe=True), enforce_file_format=True)
    dc.log_series_saved('single-frame PALETTE COLOR', name, ds, len(templates))
