import argparse
import csv
import pandas as pd
import synapseclient
import sys
import time
import urllib2
import lxml.etree as ET

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
    parser.add_argument('-inputFile', required=True, help='Name of the '\
        'file containing all accession IDs')
    parser.add_argument('-directory', required=True, default='dataFiles',
        help='Name of folder to output data to')

    return(parser)

def getAccessionIDsFromFile(filename):
    """
    Open and read the accession IDs from the input file

    Input:
        filename: The name of the input file with the accession IDs
    Return:
        List containing all of the accession IDs

    """

    # Create empty array to hold full csv file
    accessionIDs = []

    # loop through file using csv.reader and store in wholeFile
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            # Row is a single element list as only accession IDs are
            # in the file so use 0 index of row for just the name
            accessionIDs.append(row[0])

    return accessionIDs

def getGDSMetadata(gdsID):
    """
    Get the GDS metadata for the datset of interest

    Input:
        gdsID: GDS accession number of the dataset of interest
    Return:
        Dictionary containing the metadata for the entire dataset

    """

    # Registered email address with Entrez
    Entrez.email = 'rkweku@udel.edu'

    # Get the metadata for the GEO dataset
    handle = Entrez.esummary(db='gds', id=gdsID)
    metadata = Entrez.read(handle)[0]

    return(metadata)

def getSamplesFromGDS(geoMetadata):
    """
    Get all sample IDs from the GEO metadata

    Input:
        getMetadata: Dictionary with the entire metadata
    Return:
        List of all GSM accession numbers

    """

    # Initialize samples list to store all GSM accession numbers
    samples = []

    for sample in geoMetadata['Samples']:
        samples.append(sample['Accession'])

    return(samples)

def getSamplesFromOthers(accessionID):
    """
    Get the sample names from non GDS datasets

    Input:
        accessionID: The ID of the dataset
    Return:
        A list of all samples for a given dataset

    """

    samples = []

    # Complete the URL by adding the sample and then the URL cap
    # Get the initial URL 
    url = '%s' % 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
    url = '%s%s%s' % (url, accessionID, '&targ=self&form=xml&view=brief')

    try:
        usock = urllib2.urlopen(url)
    except:
        print("Throttled. Sleeping for 10 seconds")
        time.sleep(10)
        usock = urllib2.urlopen(url)

    tree = ET.parse(usock)
    usock.close()
    root = tree.getroot()

    # Get the sample elements from the xml file
    sampleElements = (root.xpath('//t:Sample',
        namespaces = {'t':'http://www.ncbi.nlm.nih.gov/geo/info/'\
        'MINiML'}))

    for element in sampleElements:
        samples.append(element.attrib['iid'])

    return(samples)

def getSampleMetadata(samples):
    """
    Get the metadata for all samples

    Input:
        samples: A list of all sample IDs as pulled from the GEO dataset
    Return:
        Dictionary containing the metadata for each sample

    """

    metadata = []

    # Go through all samples to access the metadata for each sample.
    # Append that to the final dictionary
    for sample in samples:
        sampleDict = {}
        person = {}

        # Complete the URL by adding the sample and then the URL cap
        # Get the initial URL 
        url = '%s' % 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='
        url = '%s%s%s' % (url, sample, '&targ=self&form=xml&view=brief')
 
        try:
            usock = urllib2.urlopen(url)
        except:
            print("Throttled. Sleeping for 10 seconds")
            time.sleep(10)
            usock = urllib2.urlopen(url)

        tree = ET.parse(usock)
        usock.close()
        root = tree.getroot()

        # Get the channel count from the xml file
        channelCount = int(root.xpath('//t:Channel-Count/text()', 
            namespaces = {'t':'http://www.ncbi.nlm.nih.gov/geo/info/'\
            'MINiML'})[0])

        for channel in range(channelCount):
            # Channel can't be 0. Increment by 1
            channel += 1
        # Traverse through all elements of the tree with sample information
        # for parasing.
        for child in root.find('{http://www.ncbi.nlm.nih.gov/geo/'\
                'info/MINiML}Sample').iter():
            # Save attribute, tag and text information
            attrib = child.attrib
            tag = child.tag
            # Remove MINiML info from the tag
            tag = tag.replace('{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}',
                '')
            text = child.text

            # Get the parent of the child for label checking
            child.getparent()

            # If a tag doesn't exist for the person yet, create it
            if(tag not in person.keys()):
                person[tag] = []

            # Attrib not guaranteed for each child so scheck if it exists
            if(attrib):
                # Text not guaranteed for each child so check if it exists
                if(text):
                    text = text.replace('\n', '').rstrip()
                    # If after removing excess spaces, text still exists,
                    # save tag, attrib and text as a dictionary
                    if(text):

                        # Try to set text to a float. Leave as strin g
                        # if it fails
                        try:
                            text = float(text)
                        except:
                            pass

                        # Initialize dictIndex to -1
                        dictIndex = -1
                        # Flag to record if dicitonary exists
                        dictionaryFlag = False

                        # Go through all indices in person[tag] to find
                        # if a dictionary exists
                        for i in range(len(person[tag])):
                            # If dictinary is found, set found flag to true
                            # and set dictIndex to i
                            if(type(person[tag][i]) is dict):
                                dictionaryFlag = True
                                dictIndex = i
                        
                        # If person[tag] does not exist yet, create
                        # empty dictionary in the list
                        if(not dictionaryFlag):
                            dictIndex += len(person[tag]) + 1
                            person[tag].append({})
                            dictionaryFlag = True

                        # Append to dictionary.
                        person[tag][dictIndex][attrib.get(attrib.keys()\
                            [len(attrib.keys())-1])] = text

                    # If text doesn't exist after replacing excess spaces,
                    # save only tag and attrib
                    else:
                         # Get the value of the attribute dict
                        attribValue = attrib.get(attrib.keys()[len(
                            attrib.keys())-1])

                        # Try to set attribValue to float
                        try:
                            attribValue = float(attribValue)
                        except:
                            pass
                        
                        person[tag].append(attribValue)

                # If text doesn't exist, just save the tag and attrib
                else:
                     # Get the value of the attribute dict
                    attribValue = attrib.get(attrib.keys()[len(
                        attrib.keys())-1])

                    # Try to set attribValue to float
                    try:
                        attribValue = float(attribValue)
                    except:
                        pass

                    person[tag].append(attribValue)

            # if no attrib and text, save the tag and attrib
            elif(text):
                person[tag].append(text.replace('\n','').rstrip())
        
        metadata.append(person)
    
    return(metadata)

def writeToFile(metadata, accessionID, directory):
    """
    Write the metadata to a csv file

    Input
        metadata: Dictionary containing all metadata information
        addessionID: ID of the dataset the data belongs to
        directory: Name of the folder to output data to        

    """

    # Empty dictionary to keep track of header names
    header = ['Sample ID', 'Accession ID']

    # Empty list to retain all characteristic keys
    characteristicKeys = []

    # Empty list for final output lists as dataframe
    allRows = []

    # Go through all metadata and identify all key names.
    for person in metadata:
        # Go through all keys for this specific person
        for key in person.keys():
            if(key == 'Contact-Ref' or key == 'Sample'):
                pass

            elif(key != 'Characteristics'):
                if(key not in header):
                    header.append(key)

        if('Characteristics' in person.keys()):
            # Gather all possible characteristics while here
            # Characteristic can be either a string or dictionary.
            for i in range(len(person['Characteristics'])):
                # If the type of the current characteristic is a dict,
                # add each key to the header if it does not already exist
                if(type(person['Characteristics'][i]) is dict):
                    # Go through all indices in characteristics
                    for characteristic in person['Characteristics'][i].\
                            keys():

                        # If characteristic is not yet in the
                        # list of characteristic keys, add it.
                        if(characteristic not in characteristicKeys):
                            characteristicKeys.append(characteristic)
                
                # If the type of the characteristic is not a dict,
                # it is a string, so add it to the characterList
                else:
                    if(not person['Characteristics'][i]):
                        pass
                    elif(person['Characteristics'][i] not in characteristicKeys):
                        characteristicKeys.append(person['Characteristics'][i])

    header.append('Characteristics')
    header.append('Value')

    # Go through all metadata and write to file
    for person in metadata:
        # Empty list that will hold all person lists
        dataList = []
        # Empty list to hold each person with just on charactertic
        # Each will be one row in the final output file
        personRow = []

        personRow.append(person['Sample'][len(person['Sample'])-1])
        personRow.append(accessionID)
        
        # Go through all tags in the header to ensure order
        for tag in header:
            if(tag == 'Characteristics' or tag == 'Value' or
                    tag == 'Sample' or tag == 'Contact-Ref' or
                    tag == 'Sample ID' or tag == 'Accession ID'):
                pass
            # If a tag exists for this person, it must be added to the
            # dataList
            elif(tag in person.keys()):
                featureList = person[tag]

                for feature in featureList:
                    # Need features to contain 1 element. If there is
                    # more than 1, just write list for now
                    if(len(featureList) > 1):
                        personRow.append(str(featureList))
                        break

                    elif(type(feature) == unicode):
                        personRow.append(feature.encode('utf-8'))

                    elif(type(feature) == dict):
                        personRow.append(str(feature))

                    else:
                        personRow.append(feature)

            # If tag isn't one of the blacklisted ones or characteristic
            # or value, append a comma to the personRow so that it has
            # same dimensions as the rest
            else:
                personRow.append('')

        # Characteristic is a special case so must
        # create separate loop for it
        if('Characteristics' in person.keys()):
            # Find the dictionary in the characteristics list
            # Initialize variable to -1 so loops don't happen
            # when dictionary doesn't exist
            dictIndex = -1
            for i in range(len(person['Characteristics'])):
                if(type(person['Characteristics'][i]) is dict):
                    dictIndex = i

            # save the character list for this specific person
            personCharList = person['Characteristics']

            # Go through all characteristics possible and find
            # which this person has
            for characteristic in characteristicKeys:
                # First search if the characteristic exists
                # in the base list. If it does, it is just a 
                # string not attached to a dictionary
                if(characteristic in personCharList):
                    tempList = personRow[:]
                    # Append the characteristic to the new item
                    tempList.append(characteristic)

                    # Get index of the characteristic
                    charIndex = personCharList.index(characteristic)

                    # Append the characteristic to the tempList
                    tempList.append('')
                    
                    # Append tempList to the dataList
                    dataList.append(tempList)
                    
                # If it's not in the list, it either doesn't exist
                # or is a key in the dictionary.
                else:
                    # Check if characteristic exists in the dict
                    # If it does, add it to the dataList, otherwise
                    # just skip it
                    if(dictIndex > -1 and characteristic in 
                            person['Characteristics'][dictIndex].\
                            keys()):

                        tempList = personRow[:]
                        # Append the characteristic to the new item
                        tempList.append(characteristic)
                        
                        tempList.append(person\
                            ['Characteristics'][dictIndex]\
                            [characteristic])

                        # Append tempList to the dataList
                        dataList.append(tempList)

        # If there is no characteristic for the current person
        # dataList must be updated and blanks must be assigned to 
        # Characteristic and Value
        else:
            personRow.append(['', ''])
            dataList.append(personRow)

        # Write the sample data to the allRows variable
        for row in dataList:
            allRows.append(row)

    # Create a panda dataframe for the allRows
    df = pd.DataFrame(allRows)

    # Set the columns for the dataframe to header
    df.columns = header

    # Write the metadata to the file
    df.to_csv('%s/%s.csv' % (directory, accessionID), encoding='utf-8', 
        index=False)

    # If description is in the header, limit it to 1000 characters
#    if('Description' in header):
#        df[['Description']].applymap(lambda x: x[:1000])

    return df

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
#    for col in df.columns.values:
#        print(col, df[col].map(lambda x: len(str(x))).max())

    print("Uploading to Synapse")
    schema = Schema(name='Sepsis Causing Diseases Metadata', columns=as_table_columns(df),
        parent='syn4012977')
    syn.store(Table(schema, df))

def main():
    # Parse the command line arguments to find the GEO dataset ID
    parser = cmdLineParser()
    args = parser.parse_args()

    # Get the accession IDs from the input file
    accessionIDs = getAccessionIDsFromFile(args.inputFile)

    # Create dataframe for the data
    dataFrameList = []

    # Loop through all accession IDs
    for accessionID in accessionIDs:
        metadata = []
        samples = []

        print('Currently processing %s' % accessionID)

        # If the dataset is GDS, use the GDS parsing functions
        if(accessionID.startswith('GDS')):
            # Get the GDS metadata for the given GEO ID (without 'GDS')
            gdsMetadata = getGDSMetadata(accessionID.split('GDS')[1])

            # Get all samples from the GEO 
            samples = getSamplesFromGDS(gdsMetadata)

        else:
            samples = getSamplesFromOthers(accessionID)
    
        # Only write samples if samples existed for the dataset
        if(samples):
            metadata = getSampleMetadata(samples)

        # Only write to the file if metadata existed for the dataset
        if(metadata):
            # Write to csv file
            dataFrameList.append(writeToFile(metadata, accessionID,
                args.directory))

    # Uplaod to synapse
    print('Uploading the data to Synapse')
    upload(args.directory, dataFrameList)


if __name__ == "__main__":
    main()
