
## 
##

## required: samtools bedtools in PATH

## Usage: sh Generate_BigWig_from_RNA_seq_Bam.sh your.bam 4.88
## here 4.88 is your library size 4.88 million mapped reads
## for mm10 genome version

Bamfile=$1
LibrarySize=$2
# for minus strand, we could change the library size to -4.88 to make the bigwig file

chromsize='/usr/people/bioc1387/Project/mm10/Sequences/mm10.chromsizes.txt'

Bamname=`echo $Bamfile | sed s/.bam//`

## bedGraph
echo '# genomeCoverageBed ... '
genomeCoverageBed -bg -ibam $Bamfile -g $chromsize -split > ${Bamname}.bedGraph

echo '# Normalization ... '
## normalization into 1 million mapped reads
awk -vLIBSIZE=$LibrarySize '{FS=OFS="\t"}{print $1,$2,$3,$4/LIBSIZE}' ${Bamname}.bedGraph > ${Bamname}.Norm.bedGraph

echo '# making bigwig ... '
## bedGraph to BigWig
bedGraphToBigWig ${Bamname}.Norm.bedGraph $chromsize ${Bamname}.Norm.bw

echo '# done ! '

