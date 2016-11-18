#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 12 Aug 2016 15:26:22


""" allele-specific RNA-seq analysis
    
    To take full advantage of the unassigned reads, which donot cover any SNP between two alleles, 
    we just split the reads into G1 and G2, based on the G1/G2 SNP reads ratio.
"""
import sys, os, argparse, string
from collections import defaultdict

def parse_argument():
    ''' argument parser'''
    p=argparse.ArgumentParser(description='example: %(prog)s -g1 *genome1.counts -g2 *.genome2.counts -g_all *.counts')
    p.add_argument('-g1',dest='g1',metavar='G1',type=str, required=True, help="genome1 counts file")
    p.add_argument('-g2',dest='g2',metavar='G2',type=str, required=True, help="genome2 counts file")
    p.add_argument('-g_all',dest='g_all',metavar='G_ALL',type=str, required=True, help="overall counts file, take all the G1,G2,UA,CF reads together, so G1 and G2 reads are used for calculating the two alleles' ratio, and G_all can be the normalized counts")
    if len(sys.argv) == 1 :
        sys.exit(p.print_help())
    args=p.parse_args()
    return args

def main():
    ''' main scripts '''
    args = parse_argument()
    GeneExpression = defaultdict(list)
    ### genome1
    for line in open(args.g1, "r"):
        if not line.startswith("#") and not line.startswith("_"):
            line = line.strip().split('\t')
            ######### construct the dict(GeneExpression)
            GeneExpression[line[0]].append(float(line[1]))
    ### genome2 and g_all
    filenname=[args.g2, args.g_all]
    for f in filenname:
        for line in open(f, "r"):
            if not line.startswith("#") and not line.startswith("_"):
                line = line.strip().split("\t")
                if GeneExpression.has_key(line[0]):
                    GeneExpression[line[0]].append(float(line[1]))
                else:
                    print >>sys.exit("Error happened for Gene" + line[0])
    ###########
    ## for example: Xist    123    567    5566
    print >>sys.stderr, "\n## reading over \n"
    print >>sys.stderr, "## Starting correction \n"
    fo1 = open(args.g1+"_cnts2", "w")
    fo2 = open(args.g2+"_cnts2", "w")
    for k in sorted(GeneExpression.keys()):
        #print k,v
        if GeneExpression.has_key(k):
            v = GeneExpression[k]
            if len(v)>2:
                (g1, g2, g_all) = v
                try:
                    cg1 = round(g_all * g1/(g1+g2),4)
                    cg2 = round(g_all * g2/(g1+g2),4)
                except ZeroDivisionError:
                    cg1 = cg2 = 0
                print >>fo1, "{0}\t{1}".format(k, cg1)
                print >>fo2, "{0}\t{1}".format(k, cg2)
    fo1.close()
    fo2.close()
    print >>sys.stderr, "## Done! \n" 

if __name__ == "__main__":
    main()

