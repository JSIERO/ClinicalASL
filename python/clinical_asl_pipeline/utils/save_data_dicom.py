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
import datetime
from pydicom.uid import generate_uid
from pydicom.uid import ExplicitVRLittleEndian
from pydicom.sequence import Sequence
from pydicom.dataset import Dataset
from pydicom.dataelem import DataElement
from pydicom.tag import Tag
from clinical_asl_pipeline.__version__ import __version__ as TOOL_VERSION

def set_common_metadata(ds, name, unit_str, type_tag, TOOL_VERSION):
    now = datetime.datetime.now()
    ds.ContentDate = now.strftime('%Y%m%d')
    ds.ContentTime = now.strftime('%H%M%S.%f')
    ds.ContentLabel = type_tag.upper()
    ds.ContentDescription = f"ClinicalASL-Siero derived: {type_tag.upper()} [{unit_str}]"
    ds.SeriesDescription = f"{name} [{unit_str}] - derived (ClinicalASL-Siero)"
    ds.DerivationDescription =  f"ClinicalASL-Siero derived: {type_tag.upper()} [{unit_str}]"
    code_item = Dataset()
    code_item.CodeValue = "113072"
    code_item.CodingSchemeDesignator = "DCM"
    code_item.CodeMeaning = "Post-processing"

    ds.DerivationCodeSequence = Sequence([code_item])
    ds.SoftwareVersions = f'ClinicalASL v{TOOL_VERSION}, https://github.com/JSIERO/ClinicalASL'
    ds.InstitutionName = getattr(ds, 'InstitutionName', 'University Medical Center Utrecht')
    ds.Manufacturer = f"ClinicalASL v{TOOL_VERSION}, https://github.com/JSIERO/ClinicalASL"

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

def save_data_dicom(image, source_dicom_path, output_dicom_dir, name, value_range, type_tag, series_number_incr):
    #
    # Save a 3D ASL-derived image as either a multiframe or single-frame DICOM series,
    # based on the structure of the provided reference DICOM.

    # This function scales and inserts a quantitative ASL image (e.g., CBF, CVR, AAT, ATA)
    # into a DICOM file or series, updating relevant metadata for quantitative analysis
    # and PACS compatibility.

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
    #     (min, max) specifying the value range for VOI LUT (windowing).
    # type_tag : str
    #     Label describing the quantitative map type, e.g., 'CBF', 'CVR', 'AAT', or 'ATA'.
    # series_number_incr : int
    #     Incremental value to set the SeriesNumber in the DICOM metadata.

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
    # - The function handles both signed (e.g., CVR) and unsigned (e.g., CBF, AAT, ATA) data types.
    # - The function assumes the input image is in the correct orientation and shape for ASL data.
    #

    template_dicom_path = source_dicom_path

    if not type_tag:
        raise ValueError("Please supply a type_tag such as 'CBF', 'CVR', 'AAT', or 'ATA'")
    
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
    else:
        scalingfactor = (2**16 - 1) / np.nanmax(image)
        image_scaled = (image * scalingfactor).clip(0, 2**16 - 1).astype(np.uint16)

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

            logging.info(f"Structure: {nslices} slices x {nplds} PLDs x {ndyns} dynamics x {nconditions} conditions")
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

        # Step 3: Compute normal vector for spacing
        try:
            orientation = ds.PerFrameFunctionalGroupsSequence[0].PlaneOrientationSequence[0].ImageOrientationPatient
            row_cosines = np.array(orientation[:3])
            col_cosines = np.array(orientation[3:])
            normal_vector = np.cross(row_cosines, col_cosines)
        except Exception as e:
            raise ValueError(f"Failed to extract ImageOrientationPatient: {e}")

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
            ds.SliceThickness = float(np.mean(spacings)) if spacings else 1.0


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
        ds.PixelRepresentation = 1 if use_signed else 0

        # Choose a representative frame — e.g., the first of the selected ones
        first_frame = ds.PerFrameFunctionalGroupsSequence[selected_indices[0]]

        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"        
        smallest = int(np.min(image_scaled))
        largest = int(np.max(image_scaled))
        vr = 'SS' if use_signed else 'US'
        ds[Tag(0x0028, 0x0106)] = DataElement(Tag(0x0028, 0x0106), vr, smallest)
        ds[Tag(0x0028, 0x0107)] = DataElement(Tag(0x0028, 0x0107), vr, largest)

        set_common_metadata(ds, name, unit_str, type_tag, TOOL_VERSION)
        add_source_dicom_reference(ds, source_dicom_path)

        # Clear existing VOI LUT sequences if present
        for frame in ds.PerFrameFunctionalGroupsSequence:
            if hasattr(frame, "PixelValueTransformationSequence"):
                frame.PixelValueTransformationSequence[0].RescaleSlope = f"{1.0 / scalingfactor:.10g}"
                frame.PixelValueTransformationSequence[0].RescaleIntercept = 0.0
            if hasattr(frame, "FrameVOILUTSequence"):
                frame.FrameVOILUTSequence[0].WindowCenter = f"{np.mean(value_range):.10g}"
                frame.FrameVOILUTSequence[0].WindowWidth = f"{np.ptp(value_range):.10g}"
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
        logging.info(f"Saved multi-frame DICOM: {output_path}")

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

        # Read files and extract SliceLocation 
        file_slice_pairs = []
        for f in template_files:
            path = os.path.join(template_dir, f)
            try:
                ds = pydicom.dcmread(path, force=True)
                slice_location = getattr(ds, 'SliceLocation', None)
                if slice_location is not None:
                    file_slice_pairs.append((f, slice_location))
                else:
                    logging.warning(f"Template DICOM {f} missing SliceLocation .")
            except Exception as e:
                logging.error(f"Error reading DICOM {f}: {e}")

        # Sort by slice location and pick the top N slices
        template_files_sorted = [
            f for f, _ in sorted(file_slice_pairs, key=lambda x: x[1])
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
            velocity_tag = Tag(0x0018, 0x9090)  # VelocityEncodingDirectio
            if private_tag in ds:
                sequence = ds[private_tag].value
                for item in sequence:
                    if velocity_tag in item:
                        del item[velocity_tag]

            ds.InstanceNumber = i + 1
            ds.RescaleSlope = f"{1.0 / scalingfactor:.10g}"
            ds.RescaleIntercept = 0.0
            ds.BitsAllocated = 16
            ds.BitsStored = 16
            ds.HighBit = 15
            ds.PixelRepresentation = 1 if use_signed else 0
            ds.PixelSpacing = getattr(ds, 'PixelSpacing', [1.0, 1.0])
            ds.SliceThickness = getattr(ds, 'SliceThickness', 1.0)
            ds.SamplesPerPixel = 1
            ds.PhotometricInterpretation = "MONOCHROME2"
            ds.WindowCenter = f"{np.mean(value_range):.10g}"
            ds.WindowWidth = f"{np.ptp(value_range):.10g}"

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
            logging.info(f"Saved single-frame DICOM: {output_path}")
