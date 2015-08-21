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

def cmdLineParser():
    """
    Parse the command line arguments

    Return:
        parser for all arguments

    """

    # Initialize parser
    parser = argparse.ArgumentParser()

    ### Add all arguments here
    parser.add_argument('-directory', required=True, help='Name of the '\
        'folder containing the metadata')

    return(parser)

def upload(directory, dataFrameList):
    """
    Upload the data to a Synapse table

    Input:
        directory: The name of the directory holding the data
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
    
    # Print lengths of each column
    for col in df.columns.values:
        print(col, df[col].map(lambda x: len(str(x))).max())

    print("Uploading to Synapse")
    schema = Schema(name='Diseases That Cause Sepsis Metadata', columns=as_table_columns(df),
        parent='syn4012977')
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

    for filename in os.listdir(directory):
        filename = os.path.join(directory, filename)
        if os.path.isfile(filename):
            with open(filename, 'r') as f:
                dataFrameList.append(pd.read_csv(filename))

    return(dataFrameList)

def main():
    # Parse the command line arguments to find the GEO dataset ID
    parser = cmdLineParser()
    args = parser.parse_args()

    dataFrameList = readFiles(args.directory)  

    # Upload to synapse
    print('Uploading the data to Synapse')
    upload(args.directory, dataFrameList)

if __name__ == "__main__":
    main()
