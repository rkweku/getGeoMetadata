# getGeoMetadata
======== Description ========
getGeoMetadata functions by downloading metadata and parsing the
characteristics for each sample from a GEO data series or GEO dataset. Metadata
is saved and then can be uploaded to a project on synapse.org. There are three
files within this project. 

getGeoMetadataSingle.py downloads the metadata for one dataset and saves the
metadata as a csv file

getGeoMetadata.py downloads the metadata for multiple GEO datasets and saves
the metadata as a csv file. It then uploads all of these into a specified
table on synapse.org

getGeoMetadataBySample.py downloads the metadata for one dataset and saves
the metadata in a format different than the other metadata files. It takes
all unique characteristics and creates a column for it so that each sample
only has one entry per table.

uploadToSynapse.py uploads an entire folder to synapse in one large table.
This is useful if getGeometadata.py failed at some point but a large amount
of data has already been downloaded. This can be used later to upload the
full folder once all data has been downloaded

======== Requirements ========
getGeoMetadata requires the following python libraries to perform properly:
Bio - http://biopython.org/wiki/Main_Page
lxml.etree - http://lxml.de/
pandas - http://pandas.pydata.org/
synapseclient - http://python-docs.synapse.org/

These may easily be installed using (Python) PIP. Intructions to install PIP -
https://pip.pypa.io/en/stable/installing.html

======== Execution ========
Executing any of the getGeoMetadata programs require a few specific command
line arguments for a successful execution. They are listed below along with
examples further below.

========== Notes ==========
If you are running getGeoMetada.py for multiple GEO access numbers, the
format of the inputFile will just be a list of accession IDs in plain text.
There should be one accession ID per line immediately followed by a new line.
Do not end the file with a new line. An example is provided in the repo.

======= Arguments ========
*Arguments for getGeoMetdataSingle.py and getGeoMetadataBySample.py
id		GEO accession number for the dataset of interest

*Arguments for getGeoMetadata.py
inputFile	File holding a simple list of the GEO accession numbers for
...		all datasets and series of interest
directory	Name of the folder to output the data to

*Arguments for uploadToSynapse.py
diretory	Name of the folder holding the data to be uploaded to Synapse

==============Examples ==================
1. Execution on a single dataset or data series to provide metadata with just
   one characteristic and its value per row.

python getGeoMetadataSingle.py -id GDS4071

2. Execution on multiple datasets or data series to provide metdata with just
   one characteristic and its value per row.

python getGeoMetadata.py -inputFile accessionIDs.txt -directory outputDirectory

3. Execution on a single dataset or data series to have one sample per row

python getGeoMetadataBySample.py -id GDS4071

4. Execution to upload a directory to Synapse

python uploadToSynapse.py -directory outputDirectory

======= Contact ===========
Reza Hammond
RezaKweku@gmail.com
