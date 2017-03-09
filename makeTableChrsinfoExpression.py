#!/usr/bin/python
# Programmer: Wei Guifeng, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 10 Aug 2016 23:57:40

import sys, os, argparse, string

def readGeneBed(bedfile):
    ''' '''
    isoformID2site = {}
    # for exampe geneid2site['Rnf181'] = ['chr6', 18500000]
    # left-most coordinates
    print >>sys.stderr, '## reading the genebed file'
    for line in open(bedfile, 'r'):
        if not line.startswith('#'):
            line = line.strip().split('\t')
            isoformid, chrs, start, stop, strand = line[3], line[0], int(line[1]), int(line[2]), line[5]
            isoformID2site[isoformid] = [chrs, start, stop, strand]
    print >>sys.stderr, '## reading the genebed file done ! '
    return isoformID2site


def readGeneTranscripts(GeneTranscriptsFile):
    ''' '''
    Gene2Transcripts = {}
    for line in open(GeneTranscriptsFile, 'r'):
        if not line.startswith('#'):
            line = line.strip().split('\t')
            Gene2Transcripts[line[0]] = line[1:]
    return Gene2Transcripts

def readExpFold(ExpressionFoldFile):
    ''' '''
    Gene2ExpressionFold = {}
    for line in open(ExpressionFoldFile, 'r'):
        if not line.startswith('#'):
            line = line.strip().split('\t')
            Gene2ExpressionFold[line[0]] = [ str(round(float(e),4)) for e in line[1:] ]
    return Gene2ExpressionFold


def main():
    ''' main scripts '''
    isoformID2site = readGeneBed('/usr/people/bioc1387/Project/mm10/Annotation/archive-2014-05-23-16-05-10/Genes/Genes.Genebed')
    Gene2Transcripts = readGeneTranscripts('/usr/people/bioc1387/Project/mm10/Annotation/archive-2014-05-23-16-05-10/Genes/Gene.Transcripts')
    Gene2ExpressionFold = readExpFold(sys.argv[1])
    print >>sys.stderr, '## reading the genes from certain chromosome '
    #Chr6Genes = []
    #for line in open('/usr/people/bioc1387/Project/mm10/Annotation/archive-2014-05-23-16-05-10/Genes/chr8.genes', 'r'):
    #    if not line.startswith('#'):
    #        Chr6Genes.append(line.strip())
    #print >>sys.stderr, '## reading the genes from chr6 done ! '
    for g in Gene2ExpressionFold.iterkeys():
        isoforms = Gene2Transcripts[g]
        left_coordinates = []
        right_coordinates = []
        for e in isoforms:
            left_coordinates.append( isoformID2site[e][1] )
            right_coordinates.append( isoformID2site[e][2])
        sites = min(left_coordinates)
        stop = max(right_coordinates)
        if isoformID2site[e][3] == "+":
            print '{0}\t{1}\t{2}\t{3}\t{4}'.format(isoformID2site[e][0], sites, sites+1, g, "\t".join(Gene2ExpressionFold[g]) )
        else:
            print '{0}\t{1}\t{2}\t{3}\t{4}'.format(isoformID2site[e][0], stop-1, stop, g, "\t".join(Gene2ExpressionFold[g]) )


if __name__ == "__main__":
    main()

