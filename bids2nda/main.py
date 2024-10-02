#!/usr/bin/env python
#

# import modules used here -- sys is a very standard one
from __future__ import print_function
import argparse
import csv
import logging
import zipfile
from collections import OrderedDict
from glob import glob
import os
import sys

import nibabel as nb
import json
import pandas as pd
import numpy as np


# Gather our code in a main() function
from shutil import copy


def get_metadata_for_nifti(bids_root, path):

    #TODO support .nii
    sidecarJSON = path.replace(".nii.gz", ".json")

    pathComponents = os.path.split(sidecarJSON)
    filenameComponents = pathComponents[-1].split("_")
    sessionLevelComponentList = []
    subjectLevelComponentList = []
    topLevelComponentList = []
    ses = None;
    sub = None;

    for filenameComponent in filenameComponents:
        if filenameComponent[:3] != "run":
            sessionLevelComponentList.append(filenameComponent)
            if filenameComponent[:3] == "ses":
                ses = filenameComponent
            else:
                subjectLevelComponentList.append(filenameComponent)
                if filenameComponent[:3] == "sub":
                    sub = filenameComponent
                else:
                    topLevelComponentList.append(filenameComponent)

    topLevelJSON = os.path.join(bids_root, "_".join(topLevelComponentList))
    potentialJSONs = [topLevelJSON]

    subjectLevelJSON = os.path.join(bids_root, sub, "_".join(subjectLevelComponentList))
    potentialJSONs.append(subjectLevelJSON)

    if ses:
        sessionLevelJSON = os.path.join(bids_root, sub, ses, "_".join(sessionLevelComponentList))
        potentialJSONs.append(sessionLevelJSON)

    potentialJSONs.append(sidecarJSON)

    merged_param_dict = {}
    for json_file_path in potentialJSONs:
        if os.path.exists(json_file_path):
            param_dict = json.load(open(json_file_path, "r"))
            merged_param_dict.update(param_dict)

    return merged_param_dict


def dict_append(d, guid,lookup_fields,lookup_df,key, value):
    #if this key is one we were supposed to grab from the lookup csv file
    if key in lookup_fields:
        imported_value = np.unique(lookup_df[lookup_df['subjectkey']==guid][key])
        if len(imported_value)==0:
            print(f"Did not find any {key} values for {guid} in lookup csv. You will need to fill it in manually.")
            value = ''
        elif len(imported_value)>1:
            print(f"More than one {key} value for {guid} in lookup csv. You will need to fill it in manually.")
            value = ''
        else:
            print(f"Participant {guid}: Overwriting BIDS {key} {value} with {imported_value[0].astype(type(value))} from lookup csv.")
            value = imported_value[0].astype(type(value)) #make sure the data we're inserting is the correct type
        if value!=value:
            print(f"Populated {key} value for {guid} is nan. Verify that this is correct.")
    if key in d:
        d[key].append(value)
    else:
        d[key] = [value, ]


def cosine_to_orientation(iop):
    """Deduce slicing from cosines

    From http://nipy.org/nibabel/dicom/dicom_orientation.html#dicom-voxel-to
    -patient-coordinate-system-mapping

    From Section C.7.6.1.1.1 we see that the "positive row axis" is left to
    right, and is the direction of the rows, given by the direction of last
    pixel in the first row from the first pixel in that row. Similarly the
    "positive column axis" is top to bottom and is the direction of the columns,
    given by the direction of the last pixel in the first column from the first
    pixel in that column.

    Let's rephrase: the first three values of "Image Orientation Patient" are
    the direction cosine for the "positive row axis". That is, they express the
    direction change in (x, y, z), in the DICOM patient coordinate system
    (DPCS), as you move along the row. That is, as you move from one column to
    the next. That is, as the column array index changes. Similarly, the second
    triplet of values of "Image Orientation Patient" (img_ornt_pat[3:] in
    Python), are the direction cosine for the "positive column axis", and
    express the direction you move, in the DPCS, as you move from row to row,
    and therefore as the row index changes.

    Parameters
    ----------
    iop: list of float
       Values of the ImageOrientationPatient field

    Returns
    -------
    {'Axial', 'Coronal', 'Sagittal'}
    """
    # Solution based on https://stackoverflow.com/a/45469577
    iop_round = np.round(iop)
    plane = np.cross(iop_round[0:3], iop_round[3:6])
    plane = np.abs(plane)
    if plane[0] == 1:
        return "Sagittal"
    elif plane[1] == 1:
        return "Coronal"
    elif plane[2] == 1:
        return "Axial"
    else:
        raise RuntimeError(
            "Could not deduce the image orientation of %r. 'plane' value is %r"
            % (iop, plane)
        )


def run(args):

    guid_mapping = dict([line.split(" - ") for line in open(args.guid_mapping).read().split("\n") if line != ''])
    if args.expid_mapping is not None:
        expid_mapping = dict([line.split(" - ") for line in open(args.expid_mapping).read().split("\n") if line != ''])
    else:
        expid_mapping = False
    
    if (args.lookup_csv is not None) and (args.lookup_fields is not None): #if we have both a lookup csv and fields to grab
        #read lookup csv
        lookup_df = pd.read_csv(args.lookup_csv,header=0)
        lookup_fields = args.lookup_fields
        missing_fields = []
        for field in lookup_fields:
            if field not in lookup_df.columns:
                missing_fields.append(field)
        if len(missing_fields) > 0:
            missing_fields = []
            lookup_df = pd.read_csv(args.lookup_csv,header=1)
            for field in lookup_fields:
                if field not in lookup_df.columns:
                    missing_fields.append(field)
            if len(missing_fields) > 0:
                raise RuntimeError(f"Could not find these fields in the lookup csv: {' '.join(missing_fields)}")
    elif (args.lookup_csv is not None) and (args.lookup_fields is None):
        raise RuntimeError('If a lookup csv is provided, you need to provide a field or list of fields to grab')
    elif (args.lookup_csv is None) and (args.lookup_fields is not None):
        raise RuntimeError('If you specify the field(s) to grab, you must provide a lookup csv (i.e. ndar_subject01.csv)')
    else:
        lookup_fields = []
        lookup_df = []

    suffix_to_scan_type = {"dwi": "MR diffusion",
                           "bold": "fMRI",
                           "sbref": "fMRI",
                           #""MR structural(MPRAGE)",
                           "T1w": "MR structural (T1)",
                           "PD": "MR structural (PD)",
                           #"MR structural(FSPGR)",
                           "T2w": "MR structural (T2)",
                           "inplaneT2": "MR structural (T2)",
                           "FLAIR": "FLAIR",
                           "FLASH": "MR structural (FLASH)",
                           #PET;
                            #ASL;
                            #microscopy;
                            #MR structural(PD, T2);
                            #MR structural(B0 map);
                            #MR structural(B1 map);
                            #single - shell DTI;
                            #multi - shell DTI;
                           "epi": "Field Map",
                           "phase1": "Field Map",
                           "phase2": "Field Map",
                           "phasediff": "Field Map",
                           "magnitude1": "Field Map",
                           "magnitude2": "Field Map",
                           "fieldmap": "Field Map"
                           #X - Ray
                           }

    units_dict = {"mm": "Millimeters",
                  "sec": "Seconds",
                  "msec": "Milliseconds"}

    participants_df = pd.read_csv(os.path.join(args.bids_directory, "participants.tsv"), header=0, sep="\t")

    image03_dict = OrderedDict()
    for file in glob(os.path.join(args.bids_directory, "sub-*", "*", "sub-*.nii.gz")) + \
            glob(os.path.join(args.bids_directory, "sub-*", "ses-*", "*", "sub-*_ses-*.nii.gz")):

        metadata = get_metadata_for_nifti(args.bids_directory, file)

        bids_subject_id = os.path.split(file)[-1].split("_")[0][4:]
        guid = guid_mapping[bids_subject_id]
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'subjectkey', guid)
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'src_subject_id', bids_subject_id)

        sub = file.split("sub-")[-1].split("_")[0]
        if "ses-" in file:
            ses = file.split("ses-")[-1].split("_")[0]
            scans_file = (os.path.join(args.bids_directory, "sub-" + sub, "ses-" + ses, "sub-" + sub + "_ses-" + ses + "_scans.tsv"))
        else:
            scans_file = (os.path.join(args.bids_directory, "sub-" + sub, "sub-" + sub + "_scans.tsv"))

        if os.path.exists(scans_file):
            scans_df = pd.read_csv(scans_file, header=0, sep="\t")
        else:
            print("%s file not found - information about scan date required by NDA could not be found." % scans_file)
            sys.exit(-1)
        for (_, row) in scans_df.iterrows():
            if file.endswith(row["filename"].replace("/", os.sep)):
                date = row.acq_time
                break

        sdate = date.split("-")
        ndar_date = sdate[1].zfill(2) + "/" + sdate[2].split("T")[0] + "/" + sdate[0]
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'interview_date', ndar_date)

        interview_age = int(round(list(participants_df[participants_df.participant_id == "sub-" + sub].age)[0]*12, 0))
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'interview_age', interview_age)

        sex = list(participants_df[participants_df.participant_id == "sub-" + sub].sex)[0]
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'gender', sex)

        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_file', file)

        suffix = file.split("_")[-1].split(".")[0]
        if suffix == "bold":
            # description = suffix + " " + metadata["TaskName"]
            if "_task-" in file:
                task_name = file.split("_task-")[1].split("_")[0]
            else:
                task_name = metadata["TaskName"]
            description = suffix + " " + task_name
            if expid_mapping:
                dict_append(image03_dict,guid,lookup_fields,lookup_df, 'experiment_id', expid_mapping[task_name])
            else:
                dict_append(image03_dict,guid,lookup_fields,lookup_df, 'experiment_id', "")
            dict_append(image03_dict,guid,lookup_fields,lookup_df, 'taskname', task_name)
        else:
            description = suffix
            dict_append(image03_dict,guid,lookup_fields,lookup_df, 'experiment_id', "")
            dict_append(image03_dict,guid,lookup_fields,lookup_df, 'taskname', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'abbrev_taskname', "")

        # Shortcut for the global.const section -- apparently might not be flattened fully
        metadata_const = metadata.get('global', {}).get('const', {})
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_description', description)
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'scan_type', suffix_to_scan_type[suffix])
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'scan_object', "Live")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_file_format', "NIFTI")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_modality', "MRI")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'scanner_manufacturer_pd', metadata.get("Manufacturer", ""))
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'scanner_type_pd', metadata.get("ManufacturersModelName", ""))
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'scanner_software_versions_pd', metadata.get("SoftwareVersions", ""))
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'magnetic_field_strength', metadata.get("MagneticFieldStrength", ""))
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'mri_echo_time_pd', metadata.get("EchoTime", ""))
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'flip_angle', metadata.get("FlipAngle", ""))
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'receive_coil', metadata.get("ReceiveCoilName", ""))
        # ImageOrientationPatientDICOM is populated by recent dcm2niix,
        # and ImageOrientationPatient might be provided by exhastive metadata
        # record done by heudiconv
        iop = metadata.get(
            'ImageOrientationPatientDICOM',
            metadata_const.get("ImageOrientationPatient", None)
        )
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_orientation', cosine_to_orientation(iop) if iop else '')

        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'transformation_performed', 'Yes')
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'transformation_type', 'BIDS2NDA')

        nii = nb.load(file)
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_num_dimensions', len(nii.shape))
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_extent1', nii.shape[0])
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_extent2', nii.shape[1])
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_extent3', nii.shape[2])
        if len(nii.shape) > 3:
            image_extent4 = nii.shape[3]
        else:
            image_extent4 = ""

        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_extent4', image_extent4)
        if suffix == "bold":
            extent4_type = "time"
        elif description == "epi" and len(nii.shape) == 4:
            extent4_type = "time"
        elif suffix == "dwi":
            extent4_type = "diffusion weighting"
        else:
            extent4_type = ""
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'extent4_type', extent4_type)

        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'acquisition_matrix', "%g x %g" %(nii.shape[0], nii.shape[1]))

        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_resolution1', nii.header.get_zooms()[0])
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_resolution2', nii.header.get_zooms()[1])
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_resolution3', nii.header.get_zooms()[2])
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_slice_thickness', metadata_const.get("SliceThickness", nii.header.get_zooms()[2]))
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'photomet_interpret', metadata.get("global",{}).get("const",{}).get("PhotometricInterpretation",""))
        if len(nii.shape) > 3:
            image_resolution4 = nii.header.get_zooms()[3]
        else:
            image_resolution4 = ""
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_resolution4', image_resolution4)

        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_unit1', units_dict[nii.header.get_xyzt_units()[0]])
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_unit2', units_dict[nii.header.get_xyzt_units()[0]])
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_unit3', units_dict[nii.header.get_xyzt_units()[0]])
        if len(nii.shape) > 3:
            image_unit4 = units_dict[nii.header.get_xyzt_units()[1]]
            if image_unit4 == "Milliseconds":
                TR = nii.header.get_zooms()[3]/1000.
            else:
                TR = nii.header.get_zooms()[3]
            dict_append(image03_dict,guid,lookup_fields,lookup_df, 'mri_repetition_time_pd', TR)
        else:
            image_unit4 = ""
            dict_append(image03_dict,guid,lookup_fields,lookup_df, 'mri_repetition_time_pd', metadata.get("RepetitionTime", ""))

        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'slice_timing', metadata.get("SliceTiming", ""))
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_unit4', image_unit4)

        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'mri_field_of_view_pd', "%g x %g %s" % (nii.header.get_zooms()[0],
                                                                          nii.header.get_zooms()[1],
                                                                          units_dict[nii.header.get_xyzt_units()[0]]))
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'patient_position', 'head first-supine')

        if file.split(os.sep)[-1].split("_")[1].startswith("ses"):
            visit = file.split(os.sep)[-1].split("_")[1][4:]
        else:
            visit = ""

        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'visit', visit)

        if len(metadata) > 0 or suffix in ['bold', 'dwi']:
            _, fname = os.path.split(file)
            zip_name = fname.split(".")[0] + ".metadata.zip"
            with zipfile.ZipFile(os.path.join(args.output_directory, zip_name), 'w', zipfile.ZIP_DEFLATED) as zipf:

                zipf.writestr(fname.replace(".nii.gz", ".json"), json.dumps(metadata, indent=4, sort_keys=True))
                if suffix == "bold":
                    #TODO write a more robust function for finding those files
                    events_file = file.split("_bold")[0] + "_events.tsv"
                    arch_name = os.path.split(events_file)[1]                        

                    if not os.path.exists(events_file):
                        events_file = os.path.join(args.bids_directory, "task-" + task_name + "_events.tsv")

                    if os.path.exists(events_file):
                        events_df = pd.read_csv(events_file, header=0, sep="\t")
                        if "stim_file" in events_df.columns:
                            for stim_filename in events_df.stim_file.unique().tolist():
                                stim_file = os.path.join(args.bids_directory, "stimuli",stim_filename)
                                if os.path.exists(stim_file):
                                    arc_name = os.path.split(stim_file)[-1]
                                    zipf.write(stim_file, arc_name)
                        zipf.write(events_file, arch_name)

            dict_append(image03_dict,guid,lookup_fields,lookup_df, 'data_file2', os.path.join(args.output_directory, zip_name))
            dict_append(image03_dict,guid,lookup_fields,lookup_df, 'data_file2_type', "ZIP file with additional metadata from Brain Imaging "
                                                                "Data Structure (http://bids.neuroimaging.io)")
        else:
            dict_append(image03_dict,guid,lookup_fields,lookup_df, 'data_file2', "")
            dict_append(image03_dict,guid,lookup_fields,lookup_df, 'data_file2_type', "")

        if suffix == "dwi":
            # TODO write a more robust function for finding those files
            bvec_file = file.split("_dwi")[0] + "_dwi.bvec"
            if not os.path.exists(bvec_file):
                bvec_file = os.path.join(args.bids_directory, "dwi.bvec")

            if os.path.exists(bvec_file):
                dict_append(image03_dict,guid,lookup_fields,lookup_df, 'bvecfile', bvec_file)
            else:
                dict_append(image03_dict,guid,lookup_fields,lookup_df, 'bvecfile', "")

            bval_file = file.split("_dwi")[0] + "_dwi.bval"
            if not os.path.exists(bval_file):
                bval_file = os.path.join(args.bids_directory, "dwi.bval")

            if os.path.exists(bval_file):
                dict_append(image03_dict,guid,lookup_fields,lookup_df, 'bvalfile', bval_file)
            else:
                dict_append(image03_dict,guid,lookup_fields,lookup_df, 'bvalfile', "")
            if os.path.exists(bval_file) or os.path.exists(bvec_file):
                dict_append(image03_dict,guid,lookup_fields,lookup_df, 'bvek_bval_files', 'Yes')
            else:
                dict_append(image03_dict,guid,lookup_fields,lookup_df, 'bvek_bval_files', 'No')
        else:
            dict_append(image03_dict,guid,lookup_fields,lookup_df, 'bvecfile', "")
            dict_append(image03_dict,guid,lookup_fields,lookup_df, 'bvalfile', "")
            dict_append(image03_dict,guid,lookup_fields,lookup_df, 'bvek_bval_files', "")

        # comply with image03 changes from 11/15/2023
        # https://nda.nih.gov/data_structure_history.html?short_name=image03
        
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'deviceserialnumber', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'procdate', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'visnum', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'manifest', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'emission_wavelingth', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'objective_magnification', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'objective_na', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'immersion', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'exposure_time', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'camera_sn', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'block_number', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'level', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'cut_thickness', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'stain', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'stain_details', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'pipeline_stage', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'deconvolved', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'decon_software', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'decon_method', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'psf_type', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'psf_file', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'decon_snr', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'decon_iterations', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'micro_temmplate_name', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'in_stack', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'decon_template_name', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'stack', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'slices', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'slice_number', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'slice_thickness', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'type_of_microscopy', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'comments_misc', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_thumbnail_file', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'transmit_coil', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_history', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'qc_outcome', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'qc_description', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'qc_fail_quest_reason', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'decay_correction', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'frame_end_times', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'frame_end_unit', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'frame_start_times', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'frame_start_unit', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'pet_isotope', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'pet_tracer', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'time_diff_inject_to_image', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'time_diff_units', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'pulse_seq', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'slice_acquisition', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'software_preproc', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'study', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'week', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'experiment_description', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'year_mta', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'timepoint_label', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'aqi', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'fd_mean', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'dvars_std', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'tsnr', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'fetal_age', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'fetal_age_type', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'accession_number', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'ageyears', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'iti_onset', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'stim1', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'stim2', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'stim1_side', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'stim1_magnitude', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'stim2_magnitude', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'choice_side', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'computer_choice', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'stim1_outcome', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'session_fmri', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'choice_fmri', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'outcome_fmri', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'rt_fmri', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'task__version', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'block_sv', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'trial_num', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'options_onset', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'cue_onset', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'interval_onset', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'monitor_onset', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'session_det', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'gbc', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'vtca', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'vtcan', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'eventname', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'vendor', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_extent5', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'extent5_type', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_unit5', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'image_resolution5', "")
        dict_append(image03_dict,guid,lookup_fields,lookup_df, 'excitation_wavelength', "")

    image03_df = pd.DataFrame(image03_dict)

    with open(os.path.join(args.output_directory, "image03.csv"), "w") as out_fp:
        out_fp.write('"image","3"\n')
        image03_df.to_csv(out_fp, index=False, quoting=csv.QUOTE_ALL, encoding='utf-8')

def main():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(
        description="BIDS to NDA converter, adapted from https://github.com/bids-standard/bids2nda/tree/master.",
        fromfile_prefix_chars='@')
    parser.add_argument(
        "bids_directory",
        help="Location of the root of your BIDS compatible directory",
        metavar="BIDS_DIRECTORY")
    parser.add_argument(
        "guid_mapping",
        help="Path to a text file with participant_id to GUID mapping. You will need to use the "
             "GUID Tool (https://ndar.nih.gov/contribute.html) to generate GUIDs for your participants.",
        metavar="GUID_MAPPING")
    parser.add_argument(
        "output_directory",
        help="Directory where NDA files will be stored",
        metavar="OUTPUT_DIRECTORY")
    parser.add_argument(
        "-e","--expid_mapping",
        metavar="EXPID_MAPPING",
        help="Path to a text file with experiment name to NDA experiment ID mapping.",
        required=False)
    parser.add_argument(
        "--lookup_csv",
        metavar="LOOKUP_CSV",
        help="Path to a csv with data we need for the image03.csv, i.e. ndar_subject01.csv",
        required=False)
    parser.add_argument(
        "--lookup_fields",
        nargs='+',
        help="List of column names to grab from the lookup csv. i.e. --lookup_fields interview_age sex",
        required=False)
    args = parser.parse_args()

    run(args)
    print("Metadata extraction complete.")


if __name__ == '__main__':
    main()
