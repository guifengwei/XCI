#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 24 Jul 2017 19:28:52

import sys, os, argparse, string
from collections import defaultdict


def parse_argument():
    ''' argument parser'''
    p=argparse.ArgumentParser(description='example: %(prog)s --genomesize *.tab --bedGraph *.bedGraph --bin 10000 -o result.txt')
    p.add_argument('--genomesize',dest='genomesize', metavar='Genome.txt',   required=True, type=str, help="The genome size file")
    p.add_argument('--bedGraph',  dest='bedGraph',   metavar='*.bedGraph',required=True, type=str, help="The data is stored in bedGraph format, which we want to extract")
    p.add_argument('--bin', dest='binsize',  metavar="binsize", default=10000, type=int, help="The resultion: default=%(default)d")
    p.add_argument('-o', '--output', dest="output", metavar="output", type=str, help="output file" )
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

def main():
    ''' main scripts '''
    args=parse_argument()
    out=open(args.output,'w')
    ##
    bedGraphFile  = args.bedGraph
    binsize  = args.binsize
    ####
    gsize = genomesize(GenomeSizeFile = args.genomesize)
    
    data_bedGraph = read_bedGraph(bedGraphFile = args.bedGraph, gsize = gsize)

    for chrs, chromdata in data_bedGraph.iteritems():
        n_bin = gsize[chrs]/binsize
        if gsize[chrs]%binsize != 0: n_bin += 1
        x = 1
        while x < n_bin :
            print >>out, "{0}\t{1}\t{2}\t{3}".format(chrs, (x-1)*binsize, x*binsize,  1.0*sum(data_bedGraph[chrs][ (x-1)*binsize:x*binsize ])/binsize )
            x+=1
        print >>out, "{0}\t{1}\t{2}\t{3}".format(chrs, (x-1)*binsize, gsize[chrs], 1.0*sum(data_bedGraph[chrs][ (x-1)*binsize:x*binsize ])/(gsize[chrs] - x*binsize) )
    out.close()

if __name__ == "__main__":
    main()
