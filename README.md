# XCI
Scripts used for Xist-mediated chromosome silencing analysis

# Pipeline
1. Mapping the Fastq reads into N-masked genome, where SNPs were masked by"N", with STAR_run.sh
2. Split the alignment into allelic bamfile with Split.sh
3. Count the reads number for each gene with HTseq
4. Normalisation and making expression table with python AllelicAnalysis.py configure.file
5. Calculating the Repression Score with CalculateSilencing.R
6. Permutation and Calculate the significance with Analysis_Permutaion.R
