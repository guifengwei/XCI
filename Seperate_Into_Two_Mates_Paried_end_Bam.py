#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 01 Apr 2016 10:51:35

import sys, os, argparse, string
import pysam


def main():
    ''' main scripts '''
    bamfilename = sys.argv[1]
    samfile = pysam.AlignmentFile(bamfilename, 'rb')
    outsamfile = bamfilename.strip('.bam')
    outsamfile_A  = pysam.AlignmentFile(outsamfile+'_1.bam', 'wb', template=samfile)
    outsamfile_B  = pysam.AlignmentFile(outsamfile+'_2.bam', 'wb', template=samfile)
    for row in samfile:
        # mate_pair 1
        if row.flag & 0x0040:
            n = outsamfile_A.write(row)
        elif row.flag & 0x0080:
            n = outsamfile_B.write(row)
    print >>sys.stderr, '## Split done !'
    outsamfile_A.close()
    outsamfile_B.close()

    
if __name__ == "__main__":
    if len(sys.argv) == 1 or len(sys.argv) > 2:
        print ' Usage: python Seperate_Into_Two_Mates_Paried_end_Bam.py CLIP-seq.bam'
    elif sys.argv[1][-4:] != '.bam':
        print ' Usage: python Seperate_Into_Two_Mates_Paried_end_Bam.py CLIP-seq.bam'
    else:
        main()

