# BIDS2NDA
Extract NIMH Data Archive compatible metadata from Brain Imaging Data Structure (BIDS) compatible datasets. Adapted and extended from https://github.com/bids-standard/bids2nda to also:
* zip and include stimulus files referenced in events.tsv files
* populate experiment ID column using a BIDS taskname to NDA experiment ID mapping text file

## Installation


    pip install https://github.com/brown-bnc/bids2nda/archive/master.zip


## Usage

    usage: bids2nda [-h] [-e EXPID_MAPPING] [--lookup_csv LOOKUP_CSV] [--lookup_fields LOOKUP_FIELDS [LOOKUP_FIELDS ...]] BIDS_DIRECTORY GUID_MAPPING OUTPUT_DIRECTORY
    BIDS to NDA converter.

    positional arguments:
      BIDS_DIRECTORY        Location of the root of your BIDS compatible directory
      GUID_MAPPING          Path to a text file with participant_id to GUID mapping. You will need to use the GUID Tool (https://ndar.nih.gov/contribute.html) to generate GUIDs for
                            your participants.
      OUTPUT_DIRECTORY      Directory where NDA files will be stored

    options:
      -h, --help            show this help message and exit
      -e EXPID_MAPPING, --expid_mapping EXPID_MAPPING
                            Path to a text file with experiment name to NDA experiment ID mapping.
      --lookup_csv LOOKUP_CSV
                            Path to a csv with data we need for the image03.csv, i.e. ndar_subject01.csv
      --lookup_fields LOOKUP_FIELDS [LOOKUP_FIELDS ...]
                            List of column names to grab from the lookup csv. i.e. --lookup_fields interview_age sex


## GUID_MAPPING file format
This is the file format produced by the GUID Tool, one line per subject in the format:

`<participant_id> - <GUID>`

## EXPID_MAPPING file format
This is a text file with one line per task present in the dataset in the format:

`<BIDS task name> - <NDA experiment ID number>`

The BIDS task name should match the value associated with the "task-" key in the .nii filename.

The NDA experiment ID number(s) are received from NDA after setting the study up through the NDA website [here](https://ndar.nih.gov/user/dashboard/collections.html).

## LOOKUP_CSV
This can be any csv file, including something like the `ndar_subject01.csv`, with data
you want to use to populate the `image03.csv` file.

This is useful in cases where the BIDS-derived data is missing or incorrect, such as
the calculation of interview_age in months, which is only approximated based on age in years
from the BIDS data.

## LOOKUP_FIELDS
If you pass a lookup csv, you need to specify which columns you want to use to populate
the `image03.csv`. The column names must match, so if you're trying to populate the "interview_age"
column in the `image03.csv`, your lookup csv must have a column with the same name.

You can pass more than one value, separated by spaces, like:
`--lookup_fields interview_age gender`.

## Example outputs
See [/examples](/examples)

