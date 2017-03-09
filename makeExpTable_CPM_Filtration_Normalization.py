#!/usr/bin/python
# Programmer: Wei Guifeng, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 05 Dec 2016 12:09:37

import sys, os, argparse, string
from collections import defaultdict

'''
    pipeline
    1st, Expression Table making, according to the htseq-count, just making one table for analysis
    2nd, CPM calculation. 
    3rd, filtration of those lowly exprssed genes
    4th, Normalization (sequencing-depth adjust)
    5th, Combination into one Expression Table

    First version of the script
'''

def main():
    ''' main scripts '''
    ######################### Part I : Allele-Specific Gene Counts
    print >>sys.stderr, '## Part I: Expression table making '
    GeneList = []
    GeneExpression = defaultdict(list)
    ## file order is very important
    filenames = ["CL1_NoDox.counts", "CL1_DoxA.counts", "CL1_DoxB.counts"]
    fo = open('GenesCounts.txt', 'w')
    for f in filenames:
        for line in open(f, 'r'):
            if not line.startswith('#') and not line.startswith('_'):
                line = line.strip().split('\t')
                GeneExpression[line[0]].append(line[1])
    ## output
    for k in sorted(GeneExpression.iterkeys()):
        # both Rn45s and Rn4.5s are the annotated rRNA here
        if not k.startswith('Rn45s') and not k.startswith('Rn4.5s'):
            print >>fo, '{0}\t{1}'.format( k, "\t".join(GeneExpression[k]) )
    fo.close()
    #################################################### Part II: calculating the CPM
    ### using the edgeR function: cpm()
    print >>sys.stderr, '## Part II: CPM Calculating '
    ww = os.system('Rscript ForCPM.R') ## ww is used to receive the system output, no useful meanings here.
    ##
    ###################################################### Part III: filtration of the lowly expressed genes
    print >>sys.stderr, '## Part III: Filtration '
    ## test
    cut_off = 0
    Reserved_Genes = []
    GeneExpression2 = defaultdict(list)
    fo2 = open('CPM_data_filtered.txt', 'w')
    for line in open('CPM_data.txt', 'r'):
        if not line.startswith('#'):
            line = line.strip().split('\t')
            gene_name = line[0].split('_')[0]
            if not GeneExpression2.has_key(gene_name):
                GeneExpression2[gene_name] = line[1:]
            else:
                GeneExpression2[gene_name] += line[1:]
    ### filter the genes according to the expression value
    for g, exp in GeneExpression2.iteritems():
        exp = map(float, exp)
        if sum(exp) > cut_off:
            Reserved_Genes.append(g)
    #### screen
    for line in open('CPM_data.txt', 'r'):
        if not line.startswith('#'):
            line = line.strip().split('\t')
            gene_name = line[0].split('_')[0]
            if gene_name in Reserved_Genes:
                print >>fo2, "\t".join(line)
    fo2.close()
    ################################################  Part IV: Normalization (depth-adjusted reads per million)
    ###########
    print >>sys.stderr, '## Part IV: Normalization '
    ww = os.system('Rscript ForNormalization.R')
    ###
    ############################ Part V: Expression Table combination
    # Reserved_Genes
    print >>sys.stderr, '## Part V: Expression Table combination '
    GeneExpression3 = defaultdict(list)
    fo3 = open('GeneExpression.Normalized.Counts.txt', 'w')
    for line in open('CPM_data_filtered_libsizeNormed.txt', 'r'):
        if not line.startswith('#'):
            line = line.strip().split('\t')
            gene_name = line[0].split('_')[0]
            if gene_name in Reserved_Genes:
                if not GeneExpression3.has_key(gene_name):
                    GeneExpression3[gene_name] = line[1:]
                else:
                    GeneExpression3[gene_name] += line[1:]
    for k,v in GeneExpression3.iteritems():
        print >>fo3, '{0}\t{1}'.format(k, '\t'.join(v))
    fo3.close()
    
if __name__ == "__main__":
    main()

