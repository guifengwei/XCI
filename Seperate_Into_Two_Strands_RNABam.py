#!/usr/bin/python
# Programmer: Wei Guifeng, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 12 Jan 2016 17:19:18

import sys, os, argparse, string

def parse_argument():
    ''' argument parser'''
    p=argparse.ArgumentParser(description='example: %(prog)s --bam *.bam --prefix Name', epilog="dependency: samtools")
    #
    p.add_argument('--bam', dest='bamfile',metavar='BAM', type=str, help="RNA-base sequencing, stranded. The Bam could be RNA-seq, CLIP-seq, 4sU-seq, etc.")
    #
    p.add_argument('--prefix',dest='prefix', default = 'output', type=str, help = "the prefix name for output files, the output will be Name_1.bam, Name_2.bam. default:output")
    if len(sys.argv) == 1 :
        sys.exit(p.print_help())
    args=p.parse_args()
    return args

def main():
    ''' main scripts '''
    args = parse_argument()
    bam = args.bamfile
    name = args.prefix

    ## ww no meaning here
    ww = os.system('samtools view -b -f 147 %s -o %s' %(bam, name+'_1_1.bam')  )
    ww = os.system('samtools view -b -f 99 %s -o %s ' %(bam, name+'_1_2.bam')  )
    #
    ww = os.system('samtools merge %s %s %s' %(name+'_1.bam', name+'_1_1.bam', name+'_1_2.bam') )
    ww = os.system('rm %s %s' %(name+'_1_1.bam', name+'_1_2.bam') )
    ## one strand, Normally minus_strand
    ww = os.system('samtools view -b -f 83 %s -o %s' %(bam, name+'_2_1.bam')  )
    ww = os.system('samtools view -b -f 163 %s -o %s'%(bam, name+'_2_2.bam') )
    #
    ww = os.system('samtools merge %s %s %s' %(name+'_2.bam', name+'_2_1.bam', name+'_2_2.bam') )
    ww = os.system('rm %s %s' %(name+'_2_1.bam', name+'_2_2.bam') )

if __name__ == "__main__":
    main()

