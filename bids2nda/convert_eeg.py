#!/usr/bin/env python
#

from __future__ import print_function
import argparse
import csv
import zipfile
from collections import OrderedDict
from glob import glob
import os
import sys
import datetime
from dateutil.relativedelta import relativedelta
import json
import pandas as pd
import numpy as np
from shutil import copy


def get_metadata_for_eeg(bids_root, path):

    #TODO support .nii
    sidecarJSON = path.replace(".eeg", ".json")

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


def run(args):

    guid_mapping = dict([line.split(" - ") for line in open(args.guid_mapping).read().split("\n") if line != ''])
    if args.expid_mapping is not None:
        expid_mapping = dict([line.split(" - ") for line in open(args.expid_mapping).read().split("\n") if line != ''])
    else:
        expid_mapping = False

    if args.birthdates is not None:
        birthdate_mapping = dict([line.split(" - ") for line in open(args.birthdates).read().split("\n") if line != ''])
    else:
        birthdate_mapping = False
    
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

    participants_df = pd.read_csv(os.path.join(args.bids_directory, "participants.tsv"), header=0, sep="\t")

    eeg_sub_files01_dict = OrderedDict()
    eeg_details01_dict = OrderedDict()
    for file in glob(os.path.join(args.bids_directory, "sub-*", "*", "sub-*.eeg")) + \
            glob(os.path.join(args.bids_directory, "sub-*", "ses-*", "*", "sub-*_ses-*.eeg")):

        metadata = get_metadata_for_eeg(args.bids_directory, file)

        bids_subject_id = os.path.split(file)[-1].split("_")[0][4:]
        guid = guid_mapping[bids_subject_id]
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'subjectkey', guid)
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'subjectkey', guid)

        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'src_subject_id', bids_subject_id)
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'src_subject_id', bids_subject_id)


        sub = file.split("sub-")[-1].split("_")[0]
        if "ses-" in file:
            ses = file.split("ses-")[-1].split("_")[0]
            scans_file = (os.path.join(args.bids_directory, "sub-" + sub, "ses-" + ses, "sub-" + sub + "_ses-" + ses + "_scans.tsv"))
        else:
            scans_file = (os.path.join(args.bids_directory, "sub-" + sub, "sub-" + sub + "_scans.tsv"))

        if os.path.exists(scans_file):
            scans_df = pd.read_csv(scans_file, header=0, sep="\t")
        else:
            print("%s file not found - information about data collection date required by NDA could not be found." % scans_file)
            sys.exit(-1)
        for (_, row) in scans_df.iterrows():
            if file.endswith(row["filename"].replace("/", os.sep).replace(".vhdr",".eeg")):
                date = row.acq_time
                break

        sdate = date.split("-")
        int_month = sdate[1].zfill(2)
        int_day = sdate[2].split("T")[0]
        int_year = sdate[0]
        ndar_date = int_month + "/" + int_day + "/" + int_year
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'interview_date', ndar_date)
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'interview_date', ndar_date)

        # If birthdate is provided, we can calculate exact age in months at the time of the session
        if birthdate_mapping:
            try:
                birthdate = birthdate_mapping[bids_subject_id].split("/")
                formatted_birthdate = datetime.date(int(birthdate[2]), int(birthdate[0]), int(birthdate[1]))
                formatted_int_date = datetime.date(int(int_year), int(int_month), int(int_day))

                diff = relativedelta(formatted_int_date, formatted_birthdate)
                interview_age = diff.years * 12 + diff.months
            except:
                print(f"Birthdate for subject {bids_subject_id} cannot be parsed. Check MM/DD/YYYY format. "
                      "Approximating age in months from participants.tsv.")      
                interview_age = int(round(list(participants_df[participants_df.participant_id == "sub-" + sub].age)[0]*12, 0))
        else:
            interview_age = int(round(list(participants_df[participants_df.participant_id == "sub-" + sub].age)[0]*12, 0))


        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'interview_age', interview_age)
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'interview_age', interview_age)

        sex = list(participants_df[participants_df.participant_id == "sub-" + sub].sex)[0]
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'sex', sex)
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'sex', sex)
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'comments_misc', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'capused', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'ofc', '')

        if "_task-" in file:
            task_name = file.split("_task-")[1].split("_")[0]
        else:
            task_name = metadata["TaskName"]

        if expid_mapping:
            dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'experiment_id', expid_mapping[task_name])
        else:
            dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'experiment_id', "")

        if "run-" in file:
            runnum = file.split("run-")[-1].split("_")[0]
            experiment_notes = 'Run '+ runnum + '. EEG runs within each session should be concatenated before analysis'
        else:
            experiment_notes = ''
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'experiment_notes', experiment_notes)

        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'experiment_terminated', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'experiment_validity', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'data_behavioralperformance_acc', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'data_behavioralperformance_rt', '')

        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'data_file1', file)
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'data_file1_type', ".eeg datafile")

        vhdr_name = file.replace(".eeg", ".vhdr")
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'data_file2', vhdr_name)
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'data_file2_type', ".vhdr datafile")  

        vmrk_name = file.replace(".eeg", ".vmrk")
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'data_file3', vmrk_name)
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'data_file3_type', ".vmrk datafile")      

        _, fname = os.path.split(file)
        zip_name = fname.split(".")[0] + ".metadata.zip"
        with zipfile.ZipFile(os.path.join(args.output_directory, zip_name), 'w', zipfile.ZIP_DEFLATED) as zipf:

            zipf.writestr(fname.replace(".eeg", ".json"), json.dumps(metadata, indent=4, sort_keys=True))
            #TODO write a more robust function for finding those files
            events_file = file.split("_eeg")[0] + "_events.tsv"
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
            
            events_json = file.split("_eeg")[0] + "_events.json"
            if os.path.exists(events_json):
                events_json_arch_name=os.path.split(events_json)[-1]
                zipf.write(events_json,events_json_arch_name)
            else:
                print('no events json found - not including in zipfile')
        
            physio_files=glob(file.split("_eeg")[0]+"*_physio*")
            for phys in physio_files:
                if os.path.exists(phys):
                    phys_arch_name=os.path.split(phys)[-1]
                    zipf.write(phys,phys_arch_name)
                else:
                    print('no physio files found - not including in zipfile')
            
            channels_file = file.split("_eeg")[0] + "_channels.tsv"
            if os.path.exists(channels_file):
                channels_arch_name=os.path.split(channels_file)[-1]
                zipf.write(channels_file, channels_arch_name)
            else:
                print('no channels file found - not including in zipfile')


            dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'data_file4', os.path.join(args.output_directory, zip_name))
            dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'data_file4_type', "ZIP file with additional metadata from Brain Imaging "
                                                                "Data Structure (http://bids.neuroimaging.io)")


        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'data_includedtrials', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'data_validity', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'head_circum', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'study', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'week', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'site', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'qeeg', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'pqeeg', '')        
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'visit', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'completed', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'lab020', '')


        
        # fill out the remaining fields in eeg_sub_files01_dict, mostly with placeholders so they can be filled in manually later
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'task_name', task_name)
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'image_description', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'manifest', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'visitnum', '')        
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'drug_name', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'drug_dosage', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'custom_src_id', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'mmn_fz_duration_voltage', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'mmn_fz_frequency_voltage', '')        
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'mmn_fz_both_voltage', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'mmn_cz_duration_voltage', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'mmn_cz_frequency_voltage', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'mmn_cz_both_voltage', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'timepoint_label', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'eventname', '')
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'session_number', ses)
        dict_append(eeg_sub_files01_dict,guid,lookup_fields,lookup_df, 'eeg_data_collected', '')

        #fill out the remaining fields in eeg_details01_dict, mostly with placeholders so they can be filled in manually later
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'site', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'visit', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg00a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg00b', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg001', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg003a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg003b', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg003c', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg003d', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg003e', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg003f', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg003g', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg004', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg004a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg005', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg005a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg006', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg006a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg006b', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg007', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg008a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg008b', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg008c', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg009', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg010', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg011', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg011a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg011b', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg011c', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg011d', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg012', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg012a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg013', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'head_circum', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg015', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg016', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg016a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg016b', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg017', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg018', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'electrode_impedance', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg020', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg020a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg021', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg022', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg023', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg023a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg024', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg024a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg025', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg025a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg026', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg026a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg026b', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg026c', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg027', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg027a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg027c', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg028', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg028a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg028b', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg028c', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg029', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg029a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg029c', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg030', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg030a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg030d', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg030b', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg030c', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg031', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg032', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg033', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'phne_sc_40', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegadmin_2', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegadmin_2a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegadmin_3', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegadmin_3a', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegadmin_8', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegdfr_done', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegdfr_start', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegdfr_end', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegdfr_qual', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegdfr_qual2', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegdfr_qual2_other', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegeo_start', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegeo_end', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegeo_qual2', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegeo_qual2_other', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegdpx_done', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegdpx_start', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegdpx_end', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegdpx_qual', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegdpx_qual2', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegdpx_qual2_other', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegec_start', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegec_end', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegec_qual2', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegec_qual2_other', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegcda_done', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegcda_start', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegcda_end', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegcda_qual', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegcda_qual2', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegcda_qual2_other', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegscap_done', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegscap_start', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegscap_end', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegscap_qual', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegscap_qual2', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eegscap_qual2_other', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'eeg_4', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'capfit_03', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'comments_misc', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'sseegexpd03', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'athfmed_durationvb', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'sseegnapd03', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'sseegnapt03', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'sseegnetd03', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'ssasleeptime03', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'ssemg', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'version_form', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'qeeg', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'pqeeg', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'assessment_matrics', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'numsuspectmean_chan', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'numsuspectstd_chan', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'numsuspectmean_trls', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'numsuspectstd_trls', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'numinterp', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'maxhr', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'minhr', '') 
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'badeeg', '')
        dict_append(eeg_details01_dict,guid,lookup_fields,lookup_df, 'badhr', '')

    eeg_sub_files01_df = pd.DataFrame(eeg_sub_files01_dict)
    with open(os.path.join(args.output_directory, "eeg_sub_files01.csv"), "w") as out_fp:
        out_fp.write('"eeg_sub_files","1"\n')
        eeg_sub_files01_df.to_csv(out_fp, index=False, quoting=csv.QUOTE_ALL, encoding='utf-8')

    eeg_details01_df = pd.DataFrame(eeg_details01_dict)
    with open(os.path.join(args.output_directory, "eeg_details01.csv"), "w") as out_fp:
        out_fp.write('"eeg_details","1"\n')
        eeg_details01_df.to_csv(out_fp, index=False, quoting=csv.QUOTE_ALL, encoding='utf-8')

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
        "-b","--birthdates",
        metavar="BIRTHDATE",
        help="Path to a text file with BIDS subject ID to birthdate (MM/DD/YYYY) mappings. i.e. 101 - 05/24/1990 \n"
            "If provided, this will be used to calculate participant age at the time of each session, which is required by NDA. "
            "Otherwise, approximate age in months will be calculated from the participants.tsv file. "
            "Because birthdate is personally identifiable information, it will not be included in the output csv file or zipfiles. "
            "Delete this mapping file after running the converter.",
        required=False)
    parser.add_argument(
        "--lookup_csv",
        metavar="LOOKUP_CSV",
        help="Path to a csv with data we need for the eeg_sub_files01.csv, i.e. ndar_subject01.csv",
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
