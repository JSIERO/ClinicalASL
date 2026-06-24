"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

DICOM saving utility module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Utility function to save data as DICOM files using header information from a template file.
    Uses PALETTE COLOR photometric interpretation for direct color rendering in PACS viewers.

License: BSD 3-Clause License
"""
import os
import logging
import numpy as np
import pydicom
import datetime
import matplotlib.pyplot as plt
import importlib.resources as pkg_resources
from matplotlib.colors import ListedColormap
from scipy.io import loadmat
from pydicom.uid import generate_uid
from pydicom.uid import ExplicitVRLittleEndian
from pydicom.sequence import Sequence
from pydicom.dataset import Dataset
from pydicom.dataelem import DataElement
from pydicom.tag import Tag
from clinical_asl_pipeline.__version__ import __version__ as TOOL_VERSION

# Default colormap per map type (matches save_figure_to_png usage)
DEFAULT_COLORMAPS = {
    'CBF': 'viridis',
    'CVR': 'vik',       # custom .mat, diverging
    'AAT': 'devon',     # custom .mat
    'ATA': 'viridis',
}


def load_colormap(colormap_name):
    """
    Load a colormap by name, supporting both matplotlib built-ins and
    custom .mat colormaps shipped with ClinicalASL (vik, devon).

    Matching save_figure_to_png colormap loading logic.

    Parameters
    ----------
    colormap_name : str
        'viridis', 'jet', 'vik', 'devon', or any matplotlib colormap name.

    Returns
    -------
    cmap : matplotlib.colors.Colormap
    """
    name_lower = colormap_name.lower()

    if name_lower in ['vik', 'devon']:
        try:
            with pkg_resources.files('clinical_asl_pipeline.colormaps').joinpath(f'{name_lower}.mat').open('rb') as f:
                mat = loadmat(f)
        except FileNotFoundError:
            raise ValueError(f"Missing colormap file: {name_lower}.mat")

        cmap_data = mat.get(name_lower)
        if cmap_data is None:
            raise ValueError(f"Could not find colormap data in {name_lower}.mat")
        cmap = ListedColormap(cmap_data)
    else:
        cmap = plt.get_cmap(name_lower)

    return cmap


def generate_palette_lut(cmap, scalingfactor, value_range, rescale_intercept=0.0):
    """
    Generate 65536-entry 16-bit Palette Color LUT that maps stored uint16
    pixel values through the colormap, respecting value_range for normalization.

    Index 0 is set to black (background/zero values), matching PNG output.

    Parameters
    ----------
    cmap : matplotlib.colors.Colormap
        Colormap to apply.
    scalingfactor : float
        Factor used to scale real values to stored uint16 values.
        real_value = stored_value * (1/scalingfactor) + rescale_intercept
    value_range : tuple
        (vmin, vmax) for colormap normalization (same as used for PNG output).
    rescale_intercept : float
        RescaleIntercept (nonzero for signed data offset to unsigned).

    Returns
    -------
    red, green, blue : np.ndarray (uint16, length 65536)
    """
    n = 65536
    indices = np.arange(n, dtype=np.float64)

    # Map stored pixel values back to real values
    real_values = indices / scalingfactor + rescale_intercept

    # Normalize to [0, 1] using value_range (matching PNG output)
    vmin, vmax = value_range
    if vmax == vmin:
        normalized = np.full(n, 0.5)
    else:
        normalized = np.clip((real_values - vmin) / (vmax - vmin), 0, 1)

    # Apply colormap
    colors = cmap(normalized)[:, :3]  # (65536, 3) float [0,1]

    # Set index 0 to black (background / zero values), matching PNG output
    colors[0] = [0, 0, 0]

    lut16 = (colors * 65535).astype(np.uint16)
    return lut16[:, 0], lut16[:, 1], lut16[:, 2]


def set_palette_color_tags(ds, red, green, blue):
    """
    Set DICOM Palette Color LUT tags on a dataset.

    Parameters
    ----------
    ds : pydicom.Dataset
    red, green, blue : np.ndarray (uint16, length 65536)
    """
    ds.PhotometricInterpretation = "PALETTE COLOR"
    # 0 encodes 65536 entries per DICOM convention (PS3.3 C.7.6.3.1.6)
    # Explicitly set VR to 'US' to resolve the ambiguous 'US or SS' VR for these tags
    ds[Tag(0x0028, 0x1101)] = DataElement(Tag(0x0028, 0x1101), 'US', [0, 0, 16])
    ds[Tag(0x0028, 0x1102)] = DataElement(Tag(0x0028, 0x1102), 'US', [0, 0, 16])
    ds[Tag(0x0028, 0x1103)] = DataElement(Tag(0x0028, 0x1103), 'US', [0, 0, 16])
    ds.RedPaletteColorLookupTableData = red.tobytes()
    ds.GreenPaletteColorLookupTableData = green.tobytes()
    ds.BluePaletteColorLookupTableData = blue.tobytes()


def set_common_metadata(ds, name, unit_str, type_tag, TOOL_VERSION):
    now = datetime.datetime.now()
    ds.ContentDate = now.strftime('%Y%m%d')
    ds.ContentTime = now.strftime('%H%M%S.%f')
    ds.ContentLabel = type_tag.upper()
    ds.ContentDescription = f"ClinicalASL-Siero: {type_tag.upper()} [{unit_str}]"
    ds.SeriesDescription = f"{name} - (RESEARCH ONLY - ClinicalASL)"
    ds.DerivationDescription =  f"RESEARCH ONLY - ClinicalASL-Siero: {type_tag.upper()}"
    
    # Document the derivation as post-processing in the DICOM metadata. This is important for traceability and to ensure that the DICOM files are correctly identified as derived products in PACS and other DICOM viewers.
    code_item = Dataset()
    # DCM:126302 identifies this image as a post-processed derivative.: Perfusion analysis by Arterial Spin Labeling (ASL) MR techniques, see DICOM PS3.16 2026b, Table D-1
    code_item.CodeValue = "126302"  
    code_item.CodingSchemeDesignator = "DCM"
    code_item.CodeMeaning = "Post-processing"

    ds.DerivationCodeSequence = Sequence([code_item])
    ds.SoftwareVersions = f'ClinicalASL v{TOOL_VERSION}, https://github.com/JSIERO/ClinicalASL'
    ds.InstitutionName = getattr(ds, 'InstitutionName', 'University Medical Center Utrecht')
    ds.Manufacturer = f"ClinicalASL v{TOOL_VERSION}, https://github.com/JSIERO/ClinicalASL"
    ds.ImageComments = "FOR RESEARCH PURPOSES ONLY"

    ds.SeriesDescription = ds.SeriesDescription[:64]  # VR LO max length
    ds.Manufacturer = ds.Manufacturer[:64]  # VR LO max length
    ds.SoftwareVersions = ds.SoftwareVersions[:64]  # VR LO max length
    ds.InstitutionName = ds.InstitutionName[:64]  # VR LO max length

    ds.ProtocolName = name

def add_source_dicom_reference(ds, source_dicom_path):
    # Add a ReferencedSeriesSequence to the DICOM dataset based on a source DICOM file.
    # This function reads the source DICOM file to extract SeriesInstanceUID, SOPInstanceUID,
    # and SOPClassUID, and adds them to the dataset's ReferencedSeriesSequence.
    # Parameters
    # ----------
    # ds : pydicom.dataset.Dataset
    #     The DICOM dataset to which the reference sequence will be added.
    # source_dicom_path : str
    #     Path to the source DICOM file from which to extract the reference information.        
    try:
        ref_ds = pydicom.dcmread(source_dicom_path, stop_before_pixels=True)

        ref_series_uid = getattr(ref_ds, 'SeriesInstanceUID', None)
        ref_instance_uid = getattr(ref_ds, 'SOPInstanceUID', None)
        ref_sop_class_uid = getattr(ref_ds, 'SOPClassUID', None)

        if ref_series_uid and ref_instance_uid and ref_sop_class_uid:
            referenced_instance = Dataset()
            referenced_instance.ReferencedSOPClassUID = ref_sop_class_uid
            referenced_instance.ReferencedSOPInstanceUID = ref_instance_uid

            referenced_series = Dataset()
            referenced_series.SeriesInstanceUID = ref_series_uid
            referenced_series.ReferencedInstanceSequence = Sequence([referenced_instance])

            ds.ReferencedSeriesSequence = Sequence([referenced_series])
            ds.ReferencedSOPInstanceUID = ref_instance_uid
            ds.ReferencedSOPClassUID = ref_sop_class_uid
            
        else:
            logging.warning("Source DICOM missing UID fields, skipping ReferencedSeriesSequence.")
    except Exception as e:
        logging.warning(f"Failed to add ReferencedSeriesSequence: {e}")

def save_data_dicom(image, source_dicom_path, output_dicom_dir, name, value_range, type_tag, series_number_incr, colormap_name=None, mask=None):
    #
    # Save a 3D ASL-derived image as either a multiframe or single-frame DICOM series,
    # based on the structure of the provided reference DICOM.

    # This function scales and inserts a quantitative ASL image (e.g., CBF, CVR, AAT, ATA)
    # into a DICOM file or series, updating relevant metadata for quantitative analysis
    # and PACS compatibility.

    # Uses PALETTE COLOR photometric interpretation with an embedded colormap LUT
    # for direct color rendering in PACS viewers (e.g., Sectra IDS7, MicroDICOM).
    # Quantitative values are preserved via RescaleSlope/RescaleIntercept.

    # The output format (Enhanced MR multiframe or classic single-frame DICOM) is determined
    # from the template DICOM's structure.

    # Parameters
    # ----------
    # image : np.ndarray
    #     3D ASL image array (Height, Width, Slices), typically from a NIfTI file.
    # source_dicom_path : str
    #     Path to a source DICOM file for referencing and as template in the output DICOM.
    # output_dicom_dir : str
    #     Directory where the output DICOM file(s) will be saved.
    # name : str
    #     Base for filename and Series/Protocol name in DICOM metadata.
    # value_range : tuple
    #     (min, max) specifying the value range for colormap normalization and VOI LUT.
    # type_tag : str
    #     Label describing the quantitative map type, e.g., 'CBF', 'CVR', 'AAT', or 'ATA'.
    # series_number_incr : int
    #     Incremental value to set the SeriesNumber in the DICOM metadata.
    # colormap_name : str or None
    #     Colormap name. If None, uses DEFAULT_COLORMAPS[type_tag].
    # mask : np.ndarray or None
    #     3D mask array (same shape as image). Non-brain voxels (NaN or 0) are set
    #     to pixel value 0, mapping to black in the Palette Color LUT.
    #     Required for correct background rendering of signed data (e.g., CVR).

    # Raises
    # ------
    # ValueError
    #     If type_tag is not provided or invalid.

    # Notes
    # -----
    # - The image is scaled to 16-bit integer range for DICOM compatibility.
    # - Metadata such as RescaleSlope, SeriesDescription, VOI LUT, and labeling are set.
    # - DICOM manipulation uses `pydicom`.
    # - Orientation/layout is matched to Philips DICOM conventions.
    # - Output format is chosen based on whether the template DICOM is multiframe or single-frame.
    # - For multiframe, the image is stored as a single DICOM file with multiple frames.
    # - For single-frame, each slice is saved as a separate DICOM file with InstanceNumber.
    # - For signed data (CVR), pixel values are offset to unsigned (PALETTE COLOR requirement).
    # - The function assumes the input image is in the correct orientation and shape for ASL data.
    #

    template_dicom_path = source_dicom_path

    if not type_tag:
        raise ValueError("Please supply a type_tag such as 'CBF', 'CVR', 'AAT', or 'ATA'")
    
    if colormap_name is None:
        colormap_name = DEFAULT_COLORMAPS.get(type_tag.upper(), 'viridis')

    IMPLEMENTATION_UID_ROOT = "1.3.6.1.4.1.54321.1.1" # Example root UID for ClinicalASL, fake PEN

    use_signed = type_tag.upper() == 'CVR'
    unit_str = {
        'CBF': 'ml/100g/min',
        'CVR': 'ml/100g/min',
        'ATA': 'ml/100g/min',
        'AAT': 's'
    }.get(type_tag.upper(), '')

    image = np.nan_to_num(image, nan=0.0)

    if use_signed:
        abs_max = np.max(np.abs(image))
        scalingfactor = (2**15 - 1) / abs_max
        image_scaled = (image * scalingfactor).clip(-2**15, 2**15 - 1).astype(np.int16)
        # PALETTE COLOR requires unsigned pixel data (PixelRepresentation=0)
        # Offset int16 [-32768, 32767] to uint16 [0, 65535]
        image_scaled = (image_scaled.astype(np.int32) + 32768).astype(np.uint16)
        rescale_intercept = -32768.0 / scalingfactor
    else:
        scalingfactor = (2**16 - 1) / np.nanmax(image)
        image_scaled = (image * scalingfactor).clip(0, 2**16 - 1).astype(np.uint16)
        rescale_intercept = 0.0

    # Apply mask: set non-brain voxels to 0 so they map to LUT index 0 (black)
    # This is essential for signed data (CVR) where background voxels would otherwise
    # map to the colormap midpoint after the unsigned offset.
    if mask is not None:
        brain_mask = np.isfinite(mask) & (mask != 0)
        image_scaled[~brain_mask] = 0

    # Load colormap and generate 65536-entry Palette Color LUT
    cmap = load_colormap(colormap_name)
    lut_red, lut_green, lut_blue = generate_palette_lut(cmap, scalingfactor, value_range, rescale_intercept)
    logging.info(f"Applying PALETTE COLOR with colormap '{colormap_name}', range {value_range} for {type_tag}")

    ds = pydicom.dcmread(template_dicom_path, force=True)
    is_multiframe = hasattr(ds, 'PerFrameFunctionalGroupsSequence')
    logging.info(f"Template DICOM is {'multiframe' if is_multiframe else 'single-frame'}.")

    if is_multiframe:
        if not hasattr(ds, "file_meta") or not hasattr(ds.file_meta, "TransferSyntaxUID"):
            ds.file_meta = ds.file_meta or pydicom.dataset.FileMetaDataset()
            ds.file_meta.TransferSyntaxUID = ExplicitVRLittleEndian

        # Step 1: Extract scan dimensions from private tags
        try:
            private_seq = ds.PerFrameFunctionalGroupsSequence[0].get((0x2005, 0x140f))
            nplds = int(ds.get((0x2001, 0x1017)).value)
            ndyns = int(private_seq[0].get((0x0020, 0x0105)).value)
            nslices = int(ds.get((0x2001, 0x1018)).value)
            total_frames = int(ds.NumberOfFrames)
            nconditions = total_frames // (nplds * ndyns * nslices)

            logging.info(f"Template structure: {nslices} slices x {nplds} PLDs x {ndyns} dynamics x {nconditions} conditions")
        except Exception as e:
            raise ValueError(f"Failed to extract dimensions from private tags: {e}")
        
        # Step 2: Helper to index into frames
        def get_frame_indices(pld, dynamic, condition, nslices, ndyns, nplds):
            return [
                condition * nslices * ndyns * nplds +
                s * ndyns * nplds +
                dynamic * nplds +
                pld
                for s in range(nslices)
            ]
        selected_indices = get_frame_indices(pld=0, dynamic=0, condition=0, nslices=nslices, ndyns=ndyns, nplds=nplds)

        # Step 4: Slice pixel data and assign
        image_fordicom = np.flip(np.transpose(image_scaled, (2, 1, 0)), axis=1)
        ds.PixelData = image_fordicom.tobytes()
        ds.NumberOfFrames = nslices

        # Step 5: Slice PerFrameFunctionalGroupsSequence
        original_pffs = ds.PerFrameFunctionalGroupsSequence

        # Filter frames based on Z-position
        ds.PerFrameFunctionalGroupsSequence = pydicom.sequence.Sequence(
            [original_pffs[i] for i in selected_indices]
        )
    
        # Step 6: Pixel spacing and thickness
        first_frame = ds.PerFrameFunctionalGroupsSequence[0]
        try:
            spacing = first_frame.PixelMeasuresSequence[0].PixelSpacing
            ds.PixelSpacing = [float(spacing[0]), float(spacing[1])]
        except Exception as e:
            logging.warning(f"Failed to extract PixelSpacing: {e}")
            ds.PixelSpacing = [1.0, 1.0]

        try:
            thickness = first_frame.PixelMeasuresSequence[0].SliceThickness
            ds.SliceThickness = float(thickness)
        except Exception:
            ds.SliceThickness = 1.0

        ds.SeriesInstanceUID = generate_uid(prefix=IMPLEMENTATION_UID_ROOT + '.')
        uid = generate_uid(prefix=IMPLEMENTATION_UID_ROOT + '.')
        ds.file_meta = ds.file_meta or pydicom.dataset.FileMetaDataset()
        ds.file_meta.MediaStorageSOPInstanceUID = uid
        ds.file_meta.MediaStorageSOPClassUID = pydicom.uid.UID("1.2.840.10008.5.1.4.1.1.4")
        ds.file_meta.FileMetaInformationVersion = b'\x00\x01'
        ds.file_meta.ImplementationClassUID = IMPLEMENTATION_UID_ROOT
        ds.file_meta.TransferSyntaxUID = ExplicitVRLittleEndian

        ds.SOPInstanceUID = uid
        ds.SOPClassUID = pydicom.uid.UID("1.2.840.10008.5.1.4.1.1.4.1")
        ds.StudyInstanceUID = getattr(ds, 'StudyInstanceUID', generate_uid(prefix=IMPLEMENTATION_UID_ROOT + '.'))
        ds.StudyID = getattr(ds, 'StudyID', 'Unknown')
        ds.PatientName = getattr(ds, 'PatientName', 'Anonymous^Patient')
        ds.PatientID = getattr(ds, 'PatientID', '000000')
        ds.PatientBirthDate = getattr(ds, 'PatientBirthDate', '')
        ds.PatientSex = getattr(ds, 'PatientSex', '')
        ds.ReferringPhysicianName = getattr(ds, 'ReferringPhysicianName', 'Unknown^Referring Physician')  
        ds.StudyDate = getattr(ds, 'StudyDate', '')
        ds.StudyTime = getattr(ds, 'StudyTime', '')
        ds.AccessionNumber = getattr(ds, 'AccessionNumber', '') 

        ds.SeriesNumber = int(getattr(ds, 'SeriesNumber', 0)) + 10 + series_number_incr # Increment SeriesNumber to avoid conflicts with existing series
        if type_tag.upper() == 'CVR':
            ds.SeriesNumber = 888  # Special case for CVR to avoid conflicts with CBF series

        ds.ImageType = ['DERIVED', 'SECONDARY', 'QUANTITATIVE']        
        ds.BitsAllocated = 16
        ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 0  # always unsigned for PALETTE COLOR

        # Choose a representative frame
        first_frame = ds.PerFrameFunctionalGroupsSequence[selected_indices[0]]

        ds.SamplesPerPixel = 1
        # PALETTE COLOR with embedded LUT (replaces MONOCHROME2)
        set_palette_color_tags(ds, lut_red, lut_green, lut_blue)

        smallest = int(np.min(image_scaled))
        largest = int(np.max(image_scaled))
        ds[Tag(0x0028, 0x0106)] = DataElement(Tag(0x0028, 0x0106), 'US', smallest)
        ds[Tag(0x0028, 0x0107)] = DataElement(Tag(0x0028, 0x0107), 'US', largest)

        set_common_metadata(ds, name, unit_str, type_tag, TOOL_VERSION)
        add_source_dicom_reference(ds, source_dicom_path)

        # Update per-frame sequences
        for frame in ds.PerFrameFunctionalGroupsSequence:
            if hasattr(frame, "PixelValueTransformationSequence"):
                frame.PixelValueTransformationSequence[0].RescaleSlope = f"{1.0 / scalingfactor:.10g}"
                frame.PixelValueTransformationSequence[0].RescaleIntercept = f"{rescale_intercept:.10g}"
            if hasattr(frame, "FrameVOILUTSequence"):
                frame.FrameVOILUTSequence[0].WindowCenter = "32768" # Set window center for PALETTE COLOR, set W/L to span the full uint16 stored pixel range so all 65536 LUT entries are accessible:
                frame.FrameVOILUTSequence[0].WindowWidth = "65536" #  Set window width for PALETTE COLOR, set W/L to span the full uint16 stored pixel range so all 65536 LUT entries are accessible:
                if 'CardiacSynchronizationSequence' in frame:
                    cs_seq = frame.CardiacSynchronizationSequence
                    for item in cs_seq:
                        if 'TriggerTime' in item:
                            del item.TriggerTime # Remove TriggerTime if present, not relevant for ASL-derived maps, PACS compatibility
                if 'MRVelocityEncodingSequence' in frame:
                    for item in frame.MRVelocityEncodingSequence:
                        if 'VelocityEncodingDirection' in item:
                            delattr(item, 'VelocityEncodingDirection')  # Remove VelocityEncodingDirection if present, not relevant for ASL-derived maps, PACS compatibility

        output_filename = f"{name.replace(' ', '_')}_{ds.SeriesNumber}.dcm"
        output_path = os.path.join(output_dicom_dir, f"{output_filename}")

        ds.save_as(output_path, enforce_file_format=True)
        logging.info(f"Saved multi-frame PALETTE COLOR DICOM: {output_path}")

    else:
        # --- SINGLEFRAME EXPORT ---

        num_slices_needed = image.shape[2]

        # Get directory and prefix
        template_dir = os.path.dirname(template_dicom_path)
        template_prefix = os.path.basename(template_dicom_path).rsplit('_', 1)[0] + '_'

        # Get all matching DICOM files
        template_files = [
            f for f in os.listdir(template_dir)
            if f.startswith(template_prefix) and os.path.isfile(os.path.join(template_dir, f))
        ]

        # Read files and extract ImagePostionPatient, only select files with TemporalPositionIdentifier == 1, and Phase number (private tag 2001,1008, PLD) == 1
        file_slice_pairs = []
        for f in template_files:
            path = os.path.join(template_dir, f)
            try:
                ds = pydicom.dcmread(path, force=True)
                temporal_pos = getattr(ds, 'TemporalPositionIdentifier', None)
                phase_number = ds.get((0x2001, 0x1008), None).value
                if temporal_pos == 1 and phase_number == 1:
                    image_position = getattr(ds, 'ImagePositionPatient', None)
                    image_z_coord = image_position[2]
                    if image_z_coord is not None:
                        file_slice_pairs.append((f, image_z_coord))
                    else:
                        logging.warning(f"Template DICOM {f} missing SliceLocation.")
            except Exception as e:
                logging.error(f"Error reading DICOM {f}: {e}")

        # Sort by unique slice locations and pick the top N slices
        unique_slices = {}
        for f, image_z_coord in file_slice_pairs:
            if image_z_coord not in unique_slices:
                unique_slices[image_z_coord] = f
        # Sort by slice location and select up to num_slices_needed
        template_files_sorted = [
            unique_slices[loc] for loc in sorted(unique_slices.keys())
        ][:num_slices_needed]

        series_instance_uid = generate_uid(prefix=IMPLEMENTATION_UID_ROOT + '.')

        for i, fname in enumerate(template_files_sorted):
            template_file = os.path.join(template_dir, fname)
            ds = pydicom.dcmread(template_file, force=True)
            ds.decompress()
            ds.SeriesInstanceUID = series_instance_uid
            
            uid = generate_uid(prefix=IMPLEMENTATION_UID_ROOT + '.')
            ds.file_meta = ds.file_meta or pydicom.dataset.FileMetaDataset()
            ds.file_meta.MediaStorageSOPInstanceUID = uid
            ds.file_meta.MediaStorageSOPClassUID = pydicom.uid.UID("1.2.840.10008.5.1.4.1.1.4")
            ds.file_meta.FileMetaInformationVersion = b'\x00\x01'
            ds.file_meta.ImplementationClassUID = IMPLEMENTATION_UID_ROOT
            ds.file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
            
            ds.SOPInstanceUID = uid
            ds.SOPClassUID = pydicom.uid.UID("1.2.840.10008.5.1.4.1.1.4")
            ds.StudyInstanceUID = getattr(ds, 'StudyInstanceUID', generate_uid(prefix=IMPLEMENTATION_UID_ROOT + '.'))
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

            ds.SeriesNumber = int(getattr(ds, 'SeriesNumber', 0)) + 10 + series_number_incr
            if type_tag.upper() == 'CVR':
                ds.SeriesNumber = 888  # Special case for CVR to avoid conflicts with CBF series

            if 'TriggerTime' in ds: 
                del ds.TriggerTime # Removing TriggerTime from DICOM as it is not applicable for ASL-derived map

            private_tag = Tag(0x2005, 0x140f)  # Private Per-Frame Sequence
            velocity_tag = Tag(0x0018, 0x9090)  # VelocityEncodingDirection
            if private_tag in ds:
                sequence = ds[private_tag].value
                for item in sequence:
                    if velocity_tag in item:
                        del item[velocity_tag]

            if 'NumberOfTemporalPositions' in ds: # Ensure NumberOfTemporalPositions is set to 1 or it will cause issues with PACS slice ordering
                ds.NumberOfTemporalPositions = 1
            if Tag(0x2001, 0x1081) in ds: # Ensure number of dynamics, Tag(0x2001, 0x1081) is set to 'IS' and has a value of 1
                ds[Tag(0x2001, 0x1081)] = DataElement(Tag(0x2001, 0x1081), 'IS', 1)
            if 'TemporalPositionIdentifier' in ds: # Ensure TemporalPositionIdentifier is set to 1 or it will cause issues with PACS slice ordering
                ds.TemporalPositionIdentifier = 1

            ds.InstanceNumber = i + 1
            ds.RescaleSlope = f"{1.0 / scalingfactor:.10g}"
            ds.RescaleIntercept = f"{rescale_intercept:.10g}"
            ds.BitsAllocated = 16
            ds.BitsStored = 16
            ds.HighBit = 15
            ds.PixelRepresentation = 0  # always unsigned for PALETTE COLOR
            ds.PixelSpacing = getattr(ds, 'PixelSpacing', [1.0, 1.0])
            ds.SliceThickness = getattr(ds, 'SliceThickness', 1.0)
            ds.SamplesPerPixel = 1
            # PALETTE COLOR with embedded LUT (replaces MONOCHROME2)
            set_palette_color_tags(ds, lut_red, lut_green, lut_blue)
            ds.WindowCenter = "32768" # Set window center for PALETTE COLOR, set W/L to span the full uint16 stored pixel range so all 65536 LUT entries are accessible:
            ds.WindowWidth = "65536" # Set window width for PALETTE COLOR, set W/L to span the full uint16 stored pixel range so all 65536 LUT entries are accessible:

            slice_img = np.flipud(image_scaled[:, :, i].T)
            ds.Rows, ds.Columns = slice_img.shape
            ds.PixelData = slice_img.tobytes()

            smallest = int(np.min(slice_img))
            largest = int(np.max(slice_img))

            set_common_metadata(ds, name, unit_str, type_tag, TOOL_VERSION)
            add_source_dicom_reference(ds, source_dicom_path)

            output_filename = f"{name.replace(' ', '_')}_{ds.SeriesNumber}_{ds.InstanceNumber}.dcm"
            output_path = os.path.join(output_dicom_dir, output_filename)
            
            ds.save_as(output_path, enforce_file_format=True)
            logging.info(f"Saved single-frame PALETTE COLOR DICOM: {output_path}")
            # logging of key metadata
            logging.info(f"StudyID singleframe-derived DICOM: {ds.StudyID}")
            logging.info(f"StudyInstanceUID singleframe-derived DICOM: {ds.StudyInstanceUID}")
            logging.info(f"PatientID singleframe-derived DICOM: {ds.PatientID}")
            logging.info(f"StudyDate singleframe-derived DICOM: {ds.StudyDate}")
            logging.info(f"StudyTime singleframe-derived DICOM: {ds.StudyTime}")
