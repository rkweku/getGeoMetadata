import argparse
import pandas as pd
import synapseclient
import sys
import urllib2
import lxml.etree as ET

from Bio import Entrez
from synapseclient import Schema, Column, Table, Row, RowSet, as_table_columns

def cmdLineParser():
    """
    Parse the command line arguments

    Return:
        parser for all arguments

    """

    # Initialize parser
    parser = argparse.ArgumentParser()

    ### Add all arguments here
    parser.add_argument('-id', required=True, help='GEO accession number for '\
        'the dataset of interest')

    return(parser)

def getGDSMetadata(gdsID):
    """
    Get the GEO metadata for the datset of interest

    Input:
        gdsID: GDS accession number of the dataset of intest
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
    Get all sample IDs from the GDS metadata

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
        sleep(10)
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
                    # If after removing excess spaces, text still exists,
                    # save tag, attrib and text as a dictionary
                    if(text.replace('\n', '').rstrip()):
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
                        # empty dictionary at the end of the list
                        if(not dictionaryFlag):
                            dictIndex += len(person[tag]) + 1
                            person[tag].append({})
                            dictionaryFlag = True

                        # Append to dictionary.
                        person[tag][dictIndex][attrib.get(attrib.keys()\
                            [len(attrib.keys())-1])] = text.replace(
                            '\n','').rstrip()

                    # If text doesn't exist after replacing excess spaces,
                    # save only tag and attrib
                    else:
                        person[tag].append(attrib.get(attrib.keys()[len(
                            attrib.keys())-1]))
                # If text doesn't exist, just save the tag and attrib
                else:
                    person[tag].append(attrib.get(attrib.keys()[len(
                        attrib.keys())-1]))

            # if no attrib and text, save the tag and attrib
            elif(text):
                person[tag].append(text.replace('\n','').rstrip())
        
        metadata.append(person)
    
    return(metadata)

def writeToFile(metadata, gdsID):
    """
    Write the metadata to a csv file

    Input
        metadata: Dictionary containing all metadata information
        gdsID: ID of the GDS dataset the data belongs to

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
        personRow.append(gdsID)
        
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
                        personRow.append(featureList)
                        break

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
            dataList.append(personRow)

        for row in dataList:
            allRows.append(row)

    df = pd.DataFrame(allRows)
    df.columns = header
    df.to_csv('%s.csv' % gdsID, encoding='utf-8',  index=False)

def main():
    # Parse the command line arguments to find the GEO dataset ID
    parser = cmdLineParser()
    args = parser.parse_args()

    accessionID = args.id
    print accessionID
    print('Getting the samples from the dataset')
    if(accessionID.startswith('GDS')):
        # Get the GDS metadata for the given GEO ID (without 'GDS')
        gdsMetadata = getGDSMetadata(accessionID.split('GDS')[1])

        # Get all samples from the GEO 
        samples = getSamplesFromGDS(gdsMetadata)

    else:
        samples = getSamplesFromOthers(accessionID)

    print('Getting the metadata for all samples')
    if(samples):
        metadata = getSampleMetadata(samples)

    print('Writing the data to a file')
    # Write to csv file
    writeToFile(metadata, accessionID)

if __name__ == "__main__":
    main()
