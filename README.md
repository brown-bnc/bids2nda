# BIDS2NDA
Extract NIMH Data Archive compatible metadata from Brain Imaging Data Structure (BIDS) compatible datasets. Adapted and extended from https://github.com/bids-standard/bids2nda to also:
* zip and include stimulus files referenced in events.tsv files
* populate experiment ID column using a BIDS taskname to NDA experiment ID mapping text file

## Installation


    pip install https://github.com/brown-bnc/bids2nda/archive/master.zip


## Usage

    usage: bids2nda [-h] [-e EXPID_MAPPING] BIDS_DIRECTORY GUID_MAPPING OUTPUT_DIRECTORY

    BIDS to NDA converter.

    positional arguments:
      BIDS_DIRECTORY    Location of the root of your BIDS compatible directory.
      GUID_MAPPING      Path to a text file with participant_id to GUID mapping.
                        You will need to use the GUID Tool
                        (https://ndar.nih.gov/contribute.html) to generate GUIDs
                        for your participants.
      OUTPUT_DIRECTORY  Directory where NDA files will be stored.

    optional arguments:
      -h, --help        Show this help message and exit.
      -e EXPID_MAPPING  Path to a text file with experiment name to NDA experiment ID mapping.


## GUID_MAPPING file format
This is the file format produced by the GUID Tool, one line per subject in the format:

`<participant_id> - <GUID>`

## EXPID_MAPPING file format
This is a text file with one line per task present in the dataset in the format:

`<BIDS task name> - <NDA experiment ID number>`

The BIDS task name should match the value associated with the "task-" key in the .nii filename.

The NDA experiment ID number(s) are received from NDA after setting the study up through the NDA website [here](https://ndar.nih.gov/user/dashboard/collections.html).

## Example outputs
See [/examples](/examples)

