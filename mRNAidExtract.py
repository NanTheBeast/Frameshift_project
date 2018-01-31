#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @IDE: PyCharm Edu
# @date: Jan-09-2018 11:10 AM
# @author: nhu

"""
Usage 	python mRNAidExtract.py -h for help.
		python mRNAidExtract.py -sp -tp -refseq -o .
"""


import sys
import os
import argparse
from Bio import SeqIO

class ExtractRNAid:
	# Create a new object to store and parse fasta format refseq RNA ids.
    def __init__(self, species, rna_type):
        # Create a new object to extract total rna id
        self.species = species
        self.rnatype = rna_type

    def parseRefseq(self, refseq_filename):
    	# Parse refeq and get total ids.
        record_id = []
        for record in SeqIO.parse(refseq_filename, "fasta"):
            record_id.append(record.id)
        return record_id

    def get1type(self,recordid):
    	# get one type(NM_ NR_ XM_ XR_) of id
        record_type_id = []
        for i in recordid:
            if self.rnatype in i:
                record_type_id.append(i)
        return record_type_id

    def idCount(self,ids):
    	# Count ids in the list
        count_id = 0
        for i in ids:
            count_id += 1
        print("Number of ids:", count_id)

    def writeout(self, ids, output_directory):
    	# Write out id list
        tmpfile = open(output_directory + self.species + '_' + self.rnatype + '_id.txt', 'w')
        for item in ids:
            tmpfile.write("%s\n" % item)
        tmpfile.close()

def printIntro():
	# Introduction of the program
    print("This program extract RNA id from refseq, you need to specify ")
    print("'species', 'rna_type', 'refseq_filename', 'output_directory',")
    print("The refseq_file shoule be 'fasta' format.\n")

# def getInput():
#     # Return the parameters
#     species = raw_input("species name? ")
#     rna_type = raw_input("which RNA type: NM, NR, XM, XR?")
#     refseq_filename = raw_input("refseq filename? ")
#     output_directory = raw_input("output directory? ")
#     return species, rna_type, refseq_filename, output_directory

######################################################################################
###########################    Program starts from here    ###########################
######################################################################################
def main():
    printIntro()
#     species, rna_type, refseq_filename, output_directory = getInput()
    
    # Extract RNA id
    idlist = ExtractRNAid(species, rna_type)
    recordid = idlist.parseRefseq(refseq_filename)
    typeid = idlist.get1type(recordid)
    idlist.idCount(recordid)
    idlist.idCount(typeid)

    # Write out id_list
    idlist.writeout(typeid, output_directory)


if __name__ == "__main__":
	# Parse the arguments, users need input -sp,-tp,-refseq,-o values when they run this program. The arguments are optional, if users don't input them, the program will use default values.
    parser = argparse.ArgumentParser()
    parser.add_argument("-sp","--species", type=str, help = "species name?")
    parser.add_argument("-tp","--rna_type", type=str, help = "which RNA type: NM, NR, XM, XR?")
    parser.add_argument("-refseq","--refseq_filename", type=str, help = "refseq filename?")
    parser.add_argument("-o","--output_directory", type=str ,help = "output directory?")
    args = parser.parse_args()
    
    species = args.species
    rna_type = args.rna_type 
    refseq_filename = args.refseq_filename
    output_directory = args.output_directory
    
    # start the program
    main()
   
