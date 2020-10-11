#!/usr/bin/env python3
########################################################################
# File: samJuncs.py
#  executable: samJuncs.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 02/12/2018 Created
#
########################################################################

########################################################################
# Hot Imports & Global Variable
########################################################################
import sys
# import os
# import numpy as np
from multiprocessing import Pool
import pysam
# from tqdm import *


########################################################################
# CommandLine
########################################################################
class CommandLine(object):
    """
    Handle the command line, usage and help requests.
    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    and a standard usage and help,
    attributes:
    myCommandLine.args is a dictionary which includes each of the available command line arguments as
    myCommandLine.args['option']

    methods:

    """
    
    def __init__(self, inOpts=None):
        """
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        """
        import argparse
        self.parser = argparse.ArgumentParser(description='samJuncs.py - lorem ipsium.',
                                              add_help=True,  # default is True
                                              prefix_chars='-',
                                              usage='%(prog)s -i sorted_indexed.bam ')
        # Add args
        self.parser.add_argument('-i', '--ibam', type=str, action='store', required=True,
                                 help='Input BAM file.')
        self.parser.add_argument('-p', '--threads', action='store', required=False, default=2,
                                 help='Num threads.')
        self.parser.add_argument('--quiet', action='store_true', required=False, default=True,
                                 help='Quiet stderr output.')

        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))


class SAM(object):
    """
    Handles sequence alignment format file input and output.
    """

    def __init__(self, inFile=None, isHISAT=False, fetch=''):
        # Attributes
        self.inFile = inFile
        self.fetch = fetch

        # Start pysam object as bam
        try:
            self.reader = pysam.AlignmentFile(self.inFile, 'rb')
        except:
            # File does not exist.
            print("ERROR: Cannot find file %s. Exiting!" % self.inFile, file=sys.stderr)
            sys.exit(1)
    
        self.strandInfo = {0: '+', 16: '-'}
        # print(self.reader.find_introns((read for read in self.reader.fetch() if read.is_reverse)))

        # sys.exit(1)

        if isHISAT:
            self.inferJunctionStrand = self.inferHISATJunctionStrand
        else:
            self.inferJunctionStrand = self.inferMM2JunctionStrand

    def inferMM2JunctionStrand(self, read):
        # minimap gives junction strand denoted as 'ts'
        # the sign corresponds to the alignment orientation, where + agrees and - disagrees
        orientation = read.flag
        try:
            junction_dir = read.get_tag('ts')
        except:
            junction_dir = None

        # Try to resolve strand by looking for polyA
        if not junction_dir:
            left, right = read.cigar[0], read.cigar[-1]
            s1, s2 = read.seq[:50], read.seq[-50:]
            # pa = str()
            if ("T"*10 in s1 and left[0] == 4 and left[1] >= 10) and ("A"*10 in s2 and right[0] == 4 and right[1] >= 10):
                # probably internal priming
                junction_dir = "ambig"

            elif "T"*10 in s1 and left[0] == 4 and left[1] >= 10:
                # maps to positive strand but has a rev comp polyA
                junction_dir = "-" if orientation == 16 else "+"
                # print("anti")
                # pa = "ppa"
            elif "A"*10 in s2 and right[0] == 4 and right[1] >= 10:
                # maps to positive strand but has a sense polyA
                junction_dir = "+" if orientation == 16 else "-"
                # print("sense")
                # pa = "ppa"
            else:
                # no polyA or polyT. Fragment?
                junction_dir = "ambig"
                # pa = "nan"

        else:
            if orientation == 0 and junction_dir == "+":
                junction_dir = "+"
            elif orientation == 0 and junction_dir == "-":
                junction_dir = "-"
            elif orientation == 16 and junction_dir == "+":
                junction_dir = "-"
            elif orientation == 16 and junction_dir == "-":
                junction_dir = "+"
        return junction_dir

    def inferHISATJunctionStrand(self, read):
        # Next will be junctions
        junctions = list()
        orientation = read.flag

        if orientation == 0:
            orientation = "+"
        elif orientation == 16:
            orientation = "-"
        
        tags = read.get_tags()

        junction_dir = [x[-1] for x in tags if x[0] == 'XS']

        return junction_dir

    def readJunctions(self):
        """
        Returns start, end and junctions from a single read.
        """

        for read in self.reader.fetch():
            try:
                # Skip unmapped or multimapped reads.
                strand = self.strandInfo[read.flag]
            except:
                continue

            qName = read.query_name
            chromosome = read.reference_name
            
            refPos = read.pos
            refEnd = read.pos

            startPos = read.pos
            cigar = read.cigar
            
            # Here is the read starts
            rstart = int(read.pos)

            # Next will be junctions
            junctions = list()
            orientation = read.flag

            junction_dir = self.inferJunctionStrand(read)

            for num, flagTuple in enumerate(cigar,1):
                flag, length = flagTuple 
                if flag not in [0,2,3,7,8]:
                    continue
                    
                if flag == 3:
                    junctions.append((refEnd, refEnd+length))

                refPos = refEnd+length
                refEnd = refPos

            # Last is the end
            rend = refEnd

            yield qName, chromosome, rstart, junctions, rend, orientation, junction_dir, read.mapq


def runCMD(x):
    fetch, alignType, bam = x
    bObj = SAM(bam, alignType, fetch)
    print(bObj)
    return bObj.countJuncs()


def main():
    '''
    TDB
    '''
    myCommandLine = CommandLine()
    
    alignmentFile = myCommandLine.args['ibam']
    threads = myCommandLine.args['threads']
    quiet = myCommandLine.args['quiet']
    header = pysam.view("-H", alignmentFile).split("\n")

    alignmentCommand = header[-2]
    if "minimap" in alignmentCommand.lower():
        alignType = "mm2"
    elif "hisat" in alignmentCommand.lower():
        alignType = "his"
    else:
        print("Alignment not done using minimap2 or hisat2. Exiting.", alignmentCommand, file=sys.stderr, sep="\n")
        sys.exit(1)
    
    referenceIDs = [(i.split()[1].split(":")[-1], alignType, alignmentFile) for i in header[1:-2]]
    p = Pool(threads)

    # results = p.imap_unordered(runCMD, tqdm(referenceIDs, desc="Parsing BAM for junctions", total=len(referenceIDs)))
    results = p.map(runCMD, referenceIDs)
    print(results)

    # print(results[0])
    # for c,j in d.items():
    #     for i in j:
    #         print(c,i[0]-1, i[1], ".", d[c][i], i[-1],  sep="\t")
    

########################################################################
# Main
# Here is the main program
# 
########################################################################

if __name__ == "__main__":
    main()      
