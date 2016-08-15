
## run STAR alignment
##

## usage: sh RunSTAR.sh (reads1.fq) (reads2.fq or "") (/path/of/output/)

reads1=$1  # gzipped fastq file for read1
reads2=$2  # gzipped fastq file for read1, use "" if single-end

mismatch=2
max_multihits=1

genome_index_dir='/usr/people/bioc1387/Project/mm10/Sequences/STAR_Index/'
genes='/usr/people/bioc1387/Project/mm10/Annotation/archive-2014-05-23-16-05-10/Genes/Genes.gtf'
prefix=$3

/home/software/STAR-STAR_2.4.2a/source/STAR --runThreadN 8 --genomeDir $genome_index_dir --readFilesIn $reads1 $reads2 --readFilesCommand zcat --outSAMattributes NH HI NM MD AS XS nM --outFilterMultimapNmax $max_multihits --outFilterMismatchNmax $mismatch --alignEndsType EndToEnd --sjdbGTFfile $genes --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $prefix
