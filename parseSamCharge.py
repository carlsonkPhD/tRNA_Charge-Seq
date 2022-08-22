#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 11:35:49 2021

@author: carlsokr
"""

######################
#arguement parser for passing filenames to script from command line

import argparse
import csv

parser = argparse.ArgumentParser(description='Parse input sam files then write counts to outfile')
parser.add_argument('file1', metavar='<Input1>', type=str, nargs=1, help='name of first sam file')
#parser.add_argument('file2', metavar='<Input2>', type=str, nargs=1, help='name of second sam file')
parser.add_argument('outfile', metavar='<Outfile>', type=str, nargs=1, help='name of output file for count table')
args = parser.parse_args()
print(args.file1)
#print(args.file2)
print(args.outfile)
####################

#construct empty array to hold values
geneNames = []
myCounts = []
count = 0
with open("hg38tRNAsCCA.fa", "r") as fp:
	for line in fp:
		count += 1
		if count % 2 == 1:
			#adds gene name to gene list
			namefull = line.strip()
			name = namefull[1:]
			geneNames.append(name)

geneNames.append("*")
           
for item in geneNames:
    lstTmp = [item, 0, 0, 0]
    myCounts.append(lstTmp)





def parseSamCharge(file1, outfile):


    with open(file1,"r") as fp:
        Lines = fp.readlines()[262:]
        for index, line in enumerate(Lines):
            if index % 2 == 0:
                linesplit = line.split(maxsplit=10)
                #zero indexed value of 2 is third element readtmp
                readname = linesplit[2]
                readseq = linesplit[9]
                #if read is unmapped increase count then move on
                if readname == "*":
                    myCounts[geneNames.index(readname)][1] += 1
                else:
                    #slice last 3 characters from read and test
                    if readseq[-3:] == "CCA":
                        myCounts[geneNames.index(readname)][1] += 1
                        myCounts[geneNames.index(readname)][2] += 1
                    elif readseq[-2:] == "CC":
                        myCounts[geneNames.index(readname)][1] += 1
                        myCounts[geneNames.index(readname)][3] += 1
                    else:
                        #if last 3 chars not CCA or CCT look at pair of read
                        next_index = index + 1
                        linesplitPrev = Lines[next_index].split(maxsplit=10)
                        readnamePrev = linesplitPrev[2]
                        readseqPrev = linesplitPrev[9]
                        if readseqPrev[-3:] == "CCA":
                            myCounts[geneNames.index(readnamePrev)][1] += 1
                            myCounts[geneNames.index(readnamePrev)][2] += 1
                        elif readseqPrev[-2:] == "CC":
                            myCounts[geneNames.index(readnamePrev)][1] += 1
                            myCounts[geneNames.index(readnamePrev)][3] += 1
                        else:
                            myCounts[geneNames.index(readname)][1] += 1

            
            
    #with open(file2,"r") as fp:
    #    Lines = fp.readlines()[228:]
    #    for index, line in enumerate(Lines):
    #        if index % 2 == 1:
    #            linesplit = line.split(maxsplit=10)
                #zero indexed value of 2 is third element readtmp
    #            readname = linesplit[2]
    #            readseq = linesplit[9]
                #if read is unmapped increase count then move on
    #            if readname == "*":
    #                myCounts[geneNames.index(readname)][1] += 1
    #            else:
                    #slice last 3 characters from read and test
    #                if readseq[-3:] == "CCA":
    #                    myCounts[geneNames.index(readname)][1] += 1
    #                    myCounts[geneNames.index(readname)][2] += 1
    #                elif readseq[-3:] == "CCT":
    #                   myCounts[geneNames.index(readname)][1] += 1
    #                   myCounts[geneNames.index(readname)][3] += 1
    #                else:
                        #if last 3 chars not CCA or CCT look at pair of read
    #                    prev_index = index - 1
    #                    linesplitPrev = Lines[prev_index].split(maxsplit=10)
    #                    readnamePrev = linesplitPrev[2]
    #                    readseqPrev = linesplitPrev[9]
    #                    if readseqPrev[-3:] == "CCA":
    #                        myCounts[geneNames.index(readnamePrev)][1] += 1
    #                        myCounts[geneNames.index(readnamePrev)][2] += 1
    #                    elif readseqPrev[-3:] == "CCT":
    #                        myCounts[geneNames.index(readnamePrev)][1] += 1
    #                        myCounts[geneNames.index(readnamePrev)][3] += 1
    #                    else:
    #                        myCounts[geneNames.index(readname)][1] += 1

   
    for item in myCounts:
        print(item)
        
    with open(outfile, "w") as f:
        wr = csv.writer(f, delimiter = ":")
        wr.writerows(myCounts)
    
    return

parseSamCharge(args.file1[0], args.outfile[0])