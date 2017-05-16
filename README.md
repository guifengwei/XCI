# XCI
Scripts used for Xist-mediated chromosome silencing analysis

Pipeline
1. Fastq reads mapping with STAR_run.sh
2. Split the alignment into allelic bamfile
3. Count the reads number for each gene
4. Normalisation and making expression table with python AllelicAnalysis.py configure.file
5. Calculating the Repression Score with CalculateSilencing.R
6. Permutation and Calculate the significance.

