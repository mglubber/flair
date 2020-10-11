#!/usr/bin/env python3
########################################################################
# File: bam2Bed12.py
#  executable: bam2Bed12.py
# Purpose: Convert minimap2 aligned bam file to bed12
#
#          
# Author: Cameron M. Soulette
# History:      cms 03/22/2018 Created
#
########################################################################

########################################################################
# Hot Imports & Global Variable
########################################################################
import argparse
from flair.samJuncs import SAM
# import pysam


########################################################################
# Functions
########################################################################
def junctionsToBed12(start, end, coords):
    """
    junctionsToBed12 takes in alignment start position, end position, and genomic junction coordinates
    and converts them to start, end, and length blocks for bed12.
    """
    
    sizes, starts = [], []
    
    # coords with 0 length are reads without introns
    if len(coords) > 0:
        for num, junction in enumerate(coords, 0):
            # a junction is 2 Splice Sites
            ss1, ss2 = junction
    
            # initial start is 0
            if num == 0:
                st = 0
                size = abs(start-ss1)
            else:
                st = coords[num-1][1] - start
                size = ss1 - (st + start)
            starts.append(st)
            sizes.append(size)

        # Here is the computation for the BED end coordinate
        st = coords[-1][1] - start
        size = end - (st + start)
        starts.append(st)
        sizes.append(size)

        return len(starts), sizes, starts
    else:
        return 1, [end-start], [0] 


########################################################################
# Main
# Here is the main program
# 
########################################################################
def bamToBed12(alignment_file, positive_txn="27,158,119", negative_txn="217,95,2", unknown_txn="99,99,99"):
    """
    alignment_file: file path for the input BAM
    Default color codes for positive and negative stranded read transcripts
    positive_txn = "27,158,119"
    negative_txn = "217,95,2"
    unknown_txn = "99,99,99"
    """

    sObj = SAM(alignment_file)

    for num, readData in enumerate(sObj.readJunctions(), 0):
        read, chromosome, startPos, junctions, endPos, flags, tags, score = readData
        blocks, sizes, starts = junctionsToBed12(startPos, endPos, junctions)
        name = read + ";" + str(flags)

        if tags == "+":
            txn = positive_txn
        elif tags == "-":
            txn = negative_txn
        else:
            tags = "+" if flags == "0" else "-"
            txn = unknown_txn

        yield [chromosome, startPos, endPos, name, score, tags, startPos, endPos, txn, blocks,
               ",".join(str(x) for x in sizes) + ",", ",".join(str(x) for x in starts) + ","]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='A tool to convert minimap2 BAM to Bed12.',
                                          add_help=True,  # default is True
                                          prefix_chars='-',
                                          usage='%(prog)s -i sorted.aligned.bam ')
    # Add args
    parser.add_argument('-i', "--input_bam", action='store', required=True, help='Input bam file.')
    parser.add_argument("output", type=str, default="-", nargs=1, help="Output location. Default is stdout (-).")

    args = parser.parse_args()

    if args.output != "-":
        with open(args.output, "w") as out_file:
            [out_file.write("\t".join(record)) for record in bamToBed12(args.input_bam)]
    else:
        [print("\t".join(record)) for record in bamToBed12(args.input_bam)]
