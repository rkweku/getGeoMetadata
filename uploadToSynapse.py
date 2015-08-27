import argparse
import csv
import pandas as pd
import os
import synapseclient
import sys

from Bio import Entrez
from synapseclient import Schema, Column, Table, Row, RowSet, as_table_columns

syn=synapseclient.Synapse()
syn.login()

#############################Execution Variables###############################
# Name of the folder to be used to save the metadata files
directory = ''

# The synapse ID for the project where the table will be saved
synID = ''

# The name that you would like to name the table.
synName = ''
###############################################################################

def upload(directory, synID, synName, dataFrameList):
    """
    Upload the data to a Synapse table

    Input:
        directory: The name of the directory holding the data
        synID: Synapse ID of the project where the table will be stored
        synName: Name to be given to the new table
        dataFrameList: List of dataframes with all of the data

    """

    df = pd.DataFrame()
    print("Creating dataframe")
    for entry in dataFrameList:
        df = df.append(entry, ignore_index=True)

    # Each of these columns are longer than 1000 characters each.
    # Cut them down to 1000 chars max
    df = df.applymap(lambda x: str(x)[:1000])
    
    print("Writing to file")
    df.to_csv('%s/allData.csv' % directory, encodings='utf-8', index=False)
    
    print("Uploading to Synapse")
    schema = Schema(name=synName, columns=as_table_columns(df),
        parent=synID)
    syn.store(Table(schema, df))

def readFiles(directory):
    """
    Read all of the files from a directory

    Input:
        directory: The name of the directory holding the data

    Output:
        dataFramesList: List of all of the read files as data frames

    """

    dataFrameList = []

    # Delete merged file if it exists
    try:
        os.remove('%s/allData.csv' % directory)
    except OSError:
        pass

    for filename in os.listdir(directory):
        filename = os.path.join(directory, filename)
        if os.path.isfile(filename):
            with open(filename, 'r') as f:
                dataFrameList.append(pd.read_csv(filename))

    return(dataFrameList)

def main():
    dataFrameList = readFiles(directory)  

    # Upload to synapse
    print('Uploading the data to Synapse')
    upload(directory, synID, synName, dataFrameList)

if __name__ == "__main__":
    main()
