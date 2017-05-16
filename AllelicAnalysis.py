#/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 04 Apr 2017 17:06:45

import sys, os, argparse, string
from collections import defaultdict


def parse_cfgfile(cfg_file):
    ''' 
        parse the configure file
    '''
    NoDox_Counts_file=''
    DoxA_Counts_file=''
    DoxB_Counts_file=''
    ########
    NoDox_g1_Counts_file=''
    NoDox_g2_Counts_file=''
    #########
    DoxA_g1_Counts_file=''
    DoxA_g2_Counts_file=''
    #########
    DoxB_g1_Counts_file=''
    DoxB_g2_Counts_file=''
    
    cut_off = 0

    for line in open(cfg_file,'r'):
        if line[0] == '\n' or line[0] == '#':
            comments = line;
        else:
            line = line.rstrip();
            command = line.split('=');
            if command[0] == "NoDox":
                NoDox_Counts_file = command[1].strip()
            if command[0] == "DoxA":
                DoxA_Counts_file = command[1].strip()
            if command[0] == "DoxB":
                DoxB_Counts_file = command[1].strip()
            ### parameters
            if command[0].strip() == '_cut_off':
                cut_off = int(command[1].strip())
            ##################
            if command[0] == 'NoDox_genome1':
                NoDox_g1_Counts_file = command[1].strip()
            if command[0] == 'NoDox_genome2':
                NoDox_g2_Counts_file = command[1].strip()
            ##
            if command[0] == 'DoxA_genome1':
                DoxA_g1_Counts_file = command[1].strip()
            if command[0] == 'DoxA_genome2':
                DoxA_g2_Counts_file = command[1].strip()
            ##
            if command[0] == 'DoxB_genome1':
                DoxB_g1_Counts_file = command[1].strip()
            if command[0] == 'DoxB_genome2':
                DoxB_g2_Counts_file = command[1].strip()

    if NoDox_Counts_file=='':
        print >> sys.stderr, "No Control file! Please provide the NoDox Counts file "
        exit(1)
    if DoxA_Counts_file == '' or DoxA_Counts_file == '':
        print >> sys.stderr, "Please note: only one replicate !"
    if DoxA_g1_Counts_file == '' or DoxA_g2_Counts_file == '' or DoxB_g1_Counts_file == '' or DoxB_g2_Counts_file == '' or NoDox_g1_Counts_file == '' or NoDox_g2_Counts_file == '':
        print >>sys.stderr, " There are not any allelic Counts! "
        exit(1)
    ###
    print >>sys.stderr, '\n## Parsing the configure file! Done! '
    return NoDox_Counts_file, DoxA_Counts_file, DoxB_Counts_file, cut_off, NoDox_g1_Counts_file, NoDox_g2_Counts_file, DoxA_g1_Counts_file, DoxA_g2_Counts_file, DoxB_g1_Counts_file, DoxB_g2_Counts_file


def makeExpTable_CPM_Filtration_Normalization(NoDox_Counts_file, DoxA_Counts_file, DoxB_Counts_file):

    '''
        pipeline
        1st, Expression Table making, according to the htseq-count, just making one table for analysis
        2nd, CPM calculation. 
        3rd, filtration of those lowly exprssed genes
        4th, Normalization (sequencing-depth adjust)
        5th, Combination into one Expression Table

        First version of the script
    '''
    ######################### Part I : Allele-Specific Gene Counts
    print >>sys.stderr, '.................................................'
    print >>sys.stderr, '## Part I: Expression table making, RawCounts ... '
    GeneList = []
    GeneExpression = defaultdict(list)
    ## file order is very important
    filenames = [NoDox_Counts_file, DoxA_Counts_file, DoxB_Counts_file]
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
    ######### Now GeneExpression dict stores the RawCounts info
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
    ############################## creating the RawCounts and Normalized Counts file
    fo4 =  open('RawCounts_NormedCounts.txt', 'w')
    print >>fo4, '#GeneName\tRawCounts_NoDox\tRawCounts_DoxA\tRawCounts_DoxB\tNormedCounts_NoDox\tNormedCounts_DoxA\tNormedCounts_DoxB'
    for k,v in GeneExpression3.iteritems():
        print >>fo4, "{0}\t{1}\t{2}".format(k, '\t'.join(GeneExpression[k]), '\t'.join(v) )
    fo4.close()


def Allelic_Expression_Calling(genome1_counts='', genome2_counts='', genome_all_counts=''):
    '''
       Dividing the normed expression value on the basis of ration of genome1 to genome2
    '''
    GeneExpression = defaultdict(list)
    ### genome1
    for line in open(genome1_counts, "r"):
        if not line.startswith("#") and not line.startswith("_"):
            line = line.strip().split('\t')
            ######### construct the dict(GeneExpression)
            GeneExpression[line[0]].append(float(line[1]))
    ### genome2 and g_all
    filenname=[ genome2_counts, genome_all_counts]
    for f in filenname:
        for line in open(f, "r"):
            if not line.startswith("#") and not line.startswith("_"):
                line = line.strip().split("\t")
                if GeneExpression.has_key(line[0]):
                    GeneExpression[line[0]].append(float(line[1]))
                else:
                    print >>sys.exit(" Error happened for Gene" + line[0])
    ###########
    ## for example: Xist    123    567    5566
    print >>sys.stderr, "## reading over"
    print >>sys.stderr, "## Starting correction ... "
    fo1 = open(genome1_counts+"_cnts2", "w")
    fo2 = open(genome2_counts+"_cnts2", "w")
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
    print >>sys.stderr, "## Starting correction Done!" 


def make_allelic_Exp_Table(NoDox_G1, NoDox_G2, DoxA_G1, DoxA_G2, DoxB_G1, DoxB_G2, theshold = 0):
    ######################### Part I : Allele-Specific Gene Counts
    print >>sys.stderr, '.................................................................'
    print >>sys.stderr, '## Part I: Allelic Expression table making, normed allelic counts'
    GeneList = []
    GeneExpression = defaultdict(list)
    ## file order is very important
    filenames = [NoDox_G1, NoDox_G2, DoxA_G1, DoxA_G2, DoxB_G1, DoxB_G2]
    fo = open('GenesCounts.txt', 'w')
    for f in filenames:
        for line in open(f, 'r'):
            if not line.startswith('#') and not line.startswith('_'):
                line = line.strip().split('\t')
                if 'genome1' in f:
                    GeneExpression[line[0]+'_g1'].append(line[1])
                else:
                    GeneExpression[line[0]+'_g2'].append(line[1])
    ## output
    for k in sorted(GeneExpression.iterkeys()):
        # both Rn45s and Rn4.5s are the annotated rRNA here
        if not k.startswith('Rn45s') and not k.startswith('Rn4.5s'):
            print >>fo, '{0}\t{1}'.format( k, "\t".join(GeneExpression[k]) )
    fo.close()
    #################################################### Part II: calculating the CPM
    ### using the edgeR function: cpm()
    #print >>sys.stderr, '## Part II: CPM Calculating '
    #ww = os.system('Rscript ForCPM.R') ## ww is used to receive the system output, no useful meanings here.
    ##
    ###################################################### Part III: filtration of the lowly expressed genes
    print >>sys.stderr, '## Part III: Filtration '
    ## test
    cut_off = theshold
    Reserved_Genes = []
    GeneExpression2 = defaultdict(list)
    fo2 = open('CPM_data_filtered.txt', 'w')
    for line in open('GenesCounts.txt', 'r'):
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
    for line in open('GenesCounts.txt', 'r'):
        if not line.startswith('#'):
            line = line.strip().split('\t')
            gene_name = line[0].split('_')[0]
            if gene_name in Reserved_Genes:
                print >>fo2, "\t".join(line)
    fo2.close()
    ################################################  Part IV: Normalization (depth-adjusted reads per million)
    ###########
    #print >>sys.stderr, '## Part IV: Normalization '
    #ww = os.system('Rscript ForNormalization.R')
    ###
    ############################ Part V: Expression Table combination
    # Reversed_Genes
    print >>sys.stderr, '## Part V: Expression Table combination '
    GeneExpression3 = defaultdict(list)
    fo3 = open('GeneExpressionTable.Normalized.txt', 'w')
    for line in open('CPM_data_filtered.txt', 'r'):
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
    ##############################################
    fo4=open('RawCounts_NormedCounts_NormedAllelicExpression.txt', 'w')
    RawNormCounts = {}
    for line in open("RawCounts_NormedCounts.txt", 'r'):
        if not line.startswith("#"):
            line = line.strip().split('\t')
            RawNormCounts[line[0]] = line[1:]
    print >>fo4, '#GeneName\tRawCounts_NoDox\tRawCounts_DoxA\tRawCounts_DoxB\tNormedCounts_NoDox\tNormedCounts_DoxA\tNormedCounts_DoxB\tNoDox_G1\tNoDox_G2\tDoxA_G1\tDoxA_G2\tDoxB_G1\tDox_G2'
    for k,v in GeneExpression3.iteritems():
        print >>fo4, "{0}\t{1}\t{2}".format(k, '\t'.join(RawNormCounts[k]), '\t'.join(v) )
    fo4.close()


if __name__ == "__main__":
    Args = parse_cfgfile('configure.file')
    #print >>sys.stderr, 'Testing', Args
    makeExpTable_CPM_Filtration_Normalization(NoDox_Counts_file=Args[0], DoxA_Counts_file=Args[1], DoxB_Counts_file=Args[2])
    ww = os.system('cut -f1,5 RawCounts_NormedCounts.txt > %s' %(Args[0]+'_Normed') )
    ww = os.system('cut -f1,6 RawCounts_NormedCounts.txt > %s' %(Args[1]+'_Normed') )
    ww = os.system('cut -f1,7 RawCounts_NormedCounts.txt > %s' %(Args[2]+'_Normed') )
    ##
    ## calling of allelic expression
    print >>sys.stderr, '............................................'
    print >>sys.stderr, '## Starting to call allelic expression ... '
        
    Allelic_Expression_Calling(genome1_counts=Args[4], genome2_counts=Args[5], genome_all_counts=Args[0]+'_Normed')
    Allelic_Expression_Calling(genome1_counts=Args[6], genome2_counts=Args[7], genome_all_counts=Args[1]+'_Normed')
    Allelic_Expression_Calling(genome1_counts=Args[8], genome2_counts=Args[9], genome_all_counts=Args[2]+'_Normed')

    make_allelic_Exp_Table(NoDox_G1 = Args[4]+'_cnts2', NoDox_G2 = Args[5]+'_cnts2', DoxA_G1 = Args[6]+'_cnts2', DoxA_G2 = Args[7]+'_cnts2', DoxB_G1 = Args[8] + '_cnts2', DoxB_G2 = Args[9]+'_cnts2', theshold = Args[3])
    
    #######


