# ATAC_analyses  

## 下载文件：  

merged234.Snap.filt.Rds (983.9 MB)  

The scATA-Seq data are stored in this Seurat object. This object contains the filtered cells after QC (3611 cells) from the 3 batches. There are 2 assays (peaks and bins) with raw data in the data slot and the binarized one in the counts slot;   

此文件包含三批经过过滤的3611个细胞，两个层次cell&bin；


## upstream
用到的软件：  

* samtools 
* bedtools  
* picard
* Bwa/bowtie2


### align.dup.sh 从比对开始

```bash
# general directory where all the processing is located
dir=""
# get names of the files with extension (1 file name per cell)
# i had cram files located in the untar folder 
variable=$(ls $dir/untar/ )
# should look like this f.e.
# variable="scATAC-seq8125199.cram"

for x in $variable
do
	# get the basename of the file without extension
        basenamefile=$(basename "$dir/untar/$x" .cram)
	# use bwa to align 2 fastq files (r1 and r2) to the hg38 genome 
        bwa mem -t 8 hg38.fa $dir/fastq/$basenamefile.r1.fastq.gz $dir/fastq/$basenamefile.r2.fastq.gz | samtools view -bS - | samtools sort > $dir/bam2/$basenamefile.bam
	# get the stats for the initial reads per cell
        samtools flagstat $dir/bam2/$basenamefile.bam > $dir/log2/$basenamefile.initial.log
	# mark duplicates 
        java -Xmx30g -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=$dir/bam2/$basenamefile.bam OUTPUT=$dir/tmp/$basenamefile.markdup.bam ASSUME_SORTED=true METRICS_FILE=$dir/log2/$basenamefile.metrics.log VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false 
	# remove duplicates and bad quality reads 
        samtools view -@ 8 -F 1804 -f 2 -b $dir/tmp/$basenamefile.markdup.bam  > $dir/bam2/$basenamefile.clean.bam
	# stats after removing duplicates
        samtools flagstat $dir/bam2/$basenamefile.clean.bam > $dir/log2/$basenamefile.rmdup.log
        # indexing bam file
        samtools index $dir/bam2/$basenamefile.clean.bam
	# removing unnecessary intermediate files 
        rm $dir/bam2/$basenamefile.bam
        rm $dir/tmp/$basenamefile.markdup.bam

done 

```
### makeFragmentsBam.sh 调整格式
```bash
# directory where the processing is going on
dir=""
# get the names of the bam files
variable=$(ls $dir/bam2  | grep "\.bam$" )


for x in $variable
do
  # basename of the bam files
  basenamefile=$(basename "$dir/bam2/$x" .bam)
  # sort specifically for the bamtobed
  samtools sort -n $dir/bam2/$x > $dir/bed/$basenamefile.sort.bam
  # bam to bed
  bedtools bamtobed -bedpe -i $dir/bed/$basenamefile.sort.bam | awk -v var="$x" '{ print $1, $2, $6,var}' |  sort -k1,1 -k2,2n | uniq -c | awk '{OFS="\t"};{print $2,$3,$4,$5,$1}'  > $dir/bed/$basenamefile.bed
  # remove intermediate files
  rm $dir/bed/$basenamefile.sort.bam

done


# concatenate files from each cell
 cat *.bed > ../merged.bed
# get the paired-end reads and transform to fragments
less merged.bed | awk '{OFS="\t"};{if( $3-$2 > 0 && $3-$2 < 1000) print $1,$2+4,$3-5,$4,$5}'  > merged.filt.bed
sort -k1,1 -k2,2n merged.filt.bed > merged.filt.sort.bed
# gzip bed file
bgzip merged.filt.sort.bed
# index it
tabix -p bed merged.filt.sort.bed.gz
```


### CreateSnap.sh   

snaptools软件对象 [https://github.com/r3fang/SnapTools](https://github.com/r3fang/SnapTools)  

```bash
# sort filtered fragments 
sort -k4,4 merged.filt.sort.bed > merged.filt.sort.RG.bed
# make snap object
snaptools snap-pre  \
        --input-file=merged.filt.sort.RG.bed.gz  \
        --output-snap=merged.snap  \
        --genome-name=hg38  \
        --genome-size=hg38.chrom.sizes  \
        --min-mapq=30  \
        --min-flen=50  \
        --max-flen=1000  \
        --keep-chrm=TRUE  \
        --keep-single=FALSE  \
        --keep-secondary=False  \
        --overwrite=True  \
        --max-num=20000  \
        --min-cov=500  \
        --verbose=True

# This command creates two files .snap and .snap.qc which contains the library quality control metrics as shown below.
#  cat demo.snap.qc

# Total number of unique barcodes:             3217
# TN - Total number of fragments:              576676
# UM - Total number of uniquely mapped:        540307
# SE - Total number of single ends:            0
# SA - Total number of secondary alignments:   1
# PE - Total number of paired ends:            540306
# PP - Total number of proper paired:          539772
# PL - Total number of proper frag len:        539772
# US - Total number of usable fragments:       539772
# UQ - Total number of unique fragments:       537336
# CM - Total number of chrM fragments:         0



# add bin matrix to the snap object
# create the cell-by-bin matrix
snaptools snap-add-bmat \
    --snap-file=merged.snap \
    --bin-size-list 5000 10000 \
    --verbose=True
```


## downstream
用到的软件：  

* samtools 
* bedtools  
* picard
* Bwa/bowtie2

### IndividualSamples.Rmd  

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Matrix", force = TRUE)
BiocManager::install("tidyverse", force = TRUE)
BiocManager::install("chromVAR", force = TRUE)
BiocManager::install("motifmatchr", force = TRUE)
BiocManager::install("SummarizedExperiment", force = TRUE)
BiocManager::install("BiocParallel", force = TRUE)
BiocManager::install("TFBSTools", force = TRUE)
BiocManager::install("pheatmap", force = TRUE)
BiocManager::install("data.table", force = TRUE)
BiocManager::install("irlba", force = TRUE)
BiocManager::install("Seurat", force = TRUE)
BiocManager::install("harmony", force = TRUE)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", force = TRUE)
BiocManager::install("EnsDb.Hsapiens.v86", force = TRUE)

library(Matrix)
library(tidyverse)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BiocParallel)
library(TFBSTools)
library(pheatmap)
library(data.table)
library(irlba)
library(Seurat)
library(Signac)
library(harmony)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)

```




















































