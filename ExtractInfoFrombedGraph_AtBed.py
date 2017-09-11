#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 25 Jul 2017 17:42:38

import sys, os, argparse, string
from collections import defaultdict,namedtuple


def parse_argument():
    ''' argument parser'''
    p=argparse.ArgumentParser(description='example: %(prog)s --genomesize *.tab --bedGraph *.bedGraph --bed *.bed -o result.txt')
    p.add_argument('--genomesize',dest='genomesize', metavar='Genome.txt',   required=True, type=str, help="The genome size file")
    p.add_argument('--bedGraph',  dest='bedGraph',   metavar='*.bedGraph',required=True, type=str, help="The data is stored in bedGraph format, which we want to extract")
    p.add_argument('--bed', dest='bed',  metavar="bed", type=str, required=True, help="Provide the bed interval to retrieve the infomation from bedGraph")
    p.add_argument('-o', '--output', dest="output", metavar="output", type=str,  help="output file" )
    ###########
    if len(sys.argv) == 1 :
        sys.exit(p.print_help())
    args=p.parse_args()
    return args

def genomesize(GenomeSizeFile = ""):
    ''' 
        read genome size
    '''
    genomesize = {}
    try:
        fi = open(GenomeSizeFile, "r")
        for line in fi:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                genomesize[ line[0] ] = int(line[1])
    except:
        print >>sys.stderr, "Error: GenomeSizeFile"
    return genomesize

def read_bedGraph(bedGraphFile='', gsize = {}):
    '''
        read bedgraph and store it in python_dict
    '''
    bedGraph = defaultdict(list)
    if bedGraphFile == "" or gsize == {}:
        print >>sys.stderr, "Error: Please check the GenomeSizeFile and bedGraph data file"
    #### 
    print >>sys.stderr, "Initializing ..."
    for chrs, length in gsize.iteritems():
        bedGraph[chrs] = [0 for row in range(length)]
    for line in open(bedGraphFile, "r"):
        if not line.startswith("#"):
            line = line.strip().split("\t")
            chrs, start, stop, value = line[0], int(line[1]), int(line[2]), float(line[3])
            print >>sys.stderr, chrs, "\r",
            bedGraph[chrs][start:stop] = [value] * (stop - start)
    print 
    print >>sys.stderr, "Initializing ... Done "
    return bedGraph

def read_bed(bedfile=''):
    '''
        return the bedlist
    '''
    Bed = namedtuple('BED', ['chr', 'start', 'stop', 'name', 'score', 'strand'])
    if bedfile:
        for line in open(bedfile, "r"):
            if not line.startswith("#") and not line.startswith("track"):
                line = line.strip().split("\t")
                Bed = Bed._make(line[0:6])
                yield Bed

def main():
    ''' main scripts '''
    args=parse_argument()
    out=open(args.output,'w')
    ##
    bedGraphFile  = args.bedGraph
    ####
    gsize = genomesize(GenomeSizeFile = args.genomesize)
    
    data_bedGraph = read_bedGraph(bedGraphFile = args.bedGraph, gsize = gsize)

    for bed in read_bed(bedfile = args.bed):
        ss = 1.0*sum(data_bedGraph[bed.chr][int(bed.start):int(bed.stop)])/(  int(bed.stop) - int(bed.start)  )
        print >>out, "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(bed.chr, bed.start, bed.stop, bed.name, ss, bed.strand)
    out.close()


if __name__ == "__main__":
    main()



