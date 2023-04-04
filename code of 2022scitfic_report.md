# 复现Sci Reports中图片  
# [文章](https://pubmed.ncbi.nlm.nih.gov/36509798/)
# 数据：[GSE167531](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA704903&o=acc_s%3Aa&s=SRR13783625)  

* 取样：Cells at passage 5 were used for experiments. For the RA treatment group, 5 wells of cells were seeded on MEF feeder cells and treated with 0.5 μM RA. On the day of the scATAC-seq experiment, the cells were trypsinized and collected. Cells were purifed using FACS to remove 
MEF feeder cells, cell debris and cell aggregate, and pooled together in equal numbers (total 120,000 cells for 
library construction). For the control group, cells were cultured on the MEF feeder cells and were treated with 
vehicle (DMSO) for 48 h.  

* 测序：Libraries were sequenced on HiSeq 2000 with 150 bp paired-end reads.  

# 下载软件  

1. bap
[bapmannul](https://github.com/caleblareau/bap/wiki/Installation-&--dependencies) 前处理并不推荐用bap，确定数据后使用cell ranger 更好。

```bash
python3 -m venv venv3
source venv3/bin/activate
pip install bap-atac
bap2 --help
```

# 预处理
[参考 Bio-Rad ATAC-Seq Analysis Toolkit](https://www.bio-rad.com/webroot/web/pdf/lsr/literature/Bulletin_7191.pdf)    

FASTQ debarcoding, read trimming, alignment, bead fltration, bead deconvolution, cell fltration and peak calling

```bash
# docker安装好后，每次重新打开ubuntu都要再激活
sudo service docker start

# sequence

echo "<=== downloading sequence ===>"
prefetch SRR13783625
cd ~/data/sra
ls -lh

echo "<=== sra2fz ===>"
# paired
fastq-dump --split-3 --gzip SRR13783625.sra
gzip -dc SRR13783625.fastq.gz | head -n 8
rm *.sra

echo "<=== move to target dir ===>"
mkdir -p /mnt/d/scATAC/sci_reports2022/sequence
cd /mnt/d/scATAC/sci_reports2022/sequence
cp  ~/data/sra/SRR13783625_?.fastq.gz  /mnt/d/scATAC/sci_reports2022/sequence/
```

```bash
# genome 

echo "<=== bowtie2 genome index ===>"
mkdir -p /mnt/d/scATAC/sci_reports2022/genome
cd /mnt/d/scATAC/sci_reports2022/genome
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/mm10.zip 
unzip mm10.zip
rm mm10.zip
```

```bash
# quality control

echo "<=== 1.quality control ===>"
mkdir -p /mnt/d/scATAC/sci_reports2022/fastqc
fastqc -t 8 -o /mnt/d/scATAC/sci_reports2022/fastqc /mnt/d/scATAC/sci_reports2022/sequence/*.gz 
cd /mnt/d/scATAC/sci_reports2022/fastqc 
multiqc .

docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
bioraddbg/atac-seq-fastqc \
-i /data/ -o /data/fastqc_results


echo "<=== 2.debarcoding ===>"
cd ~/data/sra
source venv3/bin/activate
mkdir -p /mnt/d/scATAC/sci_reports2022/debarcode
cd /mnt/d/scATAC/sci_reports2022/sequence/
bap-barcode v2.1 -a ./SRR13783625_1.fastq.gz -b ./SRR13783625_2.fastq.gz --nmismatches 1 -o ../debarcode/ -c 1 

docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
bioraddbg/atac-seq-debarcode-dbg \
-i /data/ -o /data/debarcoded_reads
# 注意此处，一定命名为 *R1_*.fastq.gz & *R2_*.fastq.gz，否则无法识别
# The outputs of debarcoding are debarcoded FASTQ files. Also included is a summary report of thedebarcoding process, 
# indicating how many reads were correctly debarcoded. All reads that faileddebarcoding will have been discarded.

echo "<=== 3.read trimming ===>"
mkdir -p /mnt/d/scATAC/sci_reports2022/trimmed
cd /mnt/d/scATAC/sci_reports2022/debarcode
trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired -o /mnt/d/scATAC/sci_reports2022/trimmed/  R1.fastq.gz R2.fastq.gz  

docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
bioraddbg/atac-seq-trim-reads \
-i /data/debarcoded_reads -o /data/trimmed_reads
# 太慢了

echo "<=== 4.可选quality control again ===>"
mkdir -p /mnt/d/scATAC/sci_reports2022/fastqc/again
cd /mnt/d/scATAC/sci_reports2022/trimmed
fastqc -t 8 -o ../fastqc/again   ./*.gz 
cd /mnt/d/scATAC/sci_reports2022/fastqc/again
multiqc .
```
* [基因组网站](https://hgdownload.soe.ucsc.edu/downloads.html)
```bash
# align

## mkdir
/mnt/d/scATAC/sci_reports2022/genome
mkdir -p /mnt/d/scATAC/sci_reports2022/genome/mm10/bwa
mkdir -p /mnt/d/scATAC/sci_reports2022/genome/mm10/bwa_chrX
mkdir -p /mnt/d/scATAC/sci_reports2022/genome/mm10/annotation/blacklist
mkdir -p /mnt/d/scATAC/sci_reports2022/genome/mm10/annotation/TSS

## 下载基因组
cd /mnt/d/scATAC/sci_reports2022/genome
# 只下载1-19
for i in {1..19} 
do
  echo $i
  wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/chr$i.fa.gz
done
wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/chrX.fa.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/chrY.fa.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/chrM.fa.gz
cat $(ls *.fa.gz) > mm10.fa.gz   
#gzip -dc chrX.fa.gz >> mm10.fa
#gzip -dc chrY.fa.gz >> mm10.fa
#gzip -dc chrM.fa.gz >> mm10.fa
gunzip mm10.fa.gz

# 按照tutrial中代码
mkdir -p /mnt/d/scATAC/sci_reports2022/genome/mm10/fasta
cd /mnt/d/scATAC/sci_reports2022/genome/mm10/fasta
wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/*.fa.gz # 下载不了
cat $(ls *.fa.gz | grep -v _) > mm10.fa.gz
gunzip mm10.fa.gz

## TSS and blacklist annotation
cd /mnt/d/scATAC/sci_reports2022/genome/mm10/annotation/TSS
wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/refGene.txt.gz
gunzip refGene.txt.gz
awk '{if($4=="-"){v=$6}else{v=$5} print $3,v,v,$4}' OFS='\t' refGene.txt | awk 'length($1) <10 {print $0}' | awk '$1 != "chrY" {print $0}' | sortBed -i stdin | uniq > mm10.refgene.tss.bed
wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/chromInfo.txt.gz
gunzip chromInfo.txt.gz
bedtools slop -i mm10.refgene.tss.bed -g chromInfo.txt -b 2000 > mm10.regene.tss.2k.bed
cd $genome/mm10/annotation/blacklist  # 在ATAC中复制



## mount the genome
docker run --rm -it -v /mnt/d/scATAC/sci_reports2022/genome/:/genome/ \
--entrypoint "/bin/bash" \
bioraddbg/atac-seq-bwa
## 建索引
cd /mnt/d/scATAC/sci_reports2022/genome/mm10/bwa_all
bwa index -p bwa_index /mnt/d/scATAC/sci_reports2022/genome/mm10.fa


echo "<=== 5.reads alignment ===>"
mkdir -p /mnt/d/scATAC/sci_reports2022/alignment
bowtie2_index=/mnt/d/scATAC/sci_reports2022/genome/bwa_index
align_dir=/mnt/d/scATAC/sci_reports2022/alignment

cd /mnt/d/scATAC/sci_reports2022/trimmed
bowtie2  -p 7 -x  $bowtie2_index --very-sensitive -X 2000 -1  1.fq.gz -2 2.fq.gz \
  2>$align_dir/sample.summary \
  -S $align_dir/sample.sam 


docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
-v /mnt/d/scATAC/sci_reports2022/genome/mm10/:/genome/ \
bioraddbg/atac-seq-bwa \
-i /data/trimmed_reads \
-o /data/alignments_all/ \
-r /genome/bwa_all/
```
```bash
# alignment QC  

echo "<=== 6.alignment QC ===>"
docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
-v /mnt/d/scATAC/sci_reports2022/genome/:/genome/ \
bioraddbg/atac-seq-alignment-qc \
-i /data/alignments_all/  \
-r /genome/mm10.fa \
-o /data/alignment_qc_all
```
* bead filtration  
The output of the Alignments Tool is used as the input to Bead Filtration. lf a different alignment method isused, a directory containing a position-sorted, indexed .bam file with bead barcodes annotated as XB:Z can be substituted.  


```bash
# bead filtration 
docker run --rm -it -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
-v /mnt/d/scATAC/sci_reports2022/genome/:/genome/ \
--entrypoint "/bin/bash" \
bioraddbg/atac-seq-filter-beads

docker run -it \
	-v /mnt/d/scATAC/sci_reports2022/sequence/:/data/  \
	--entrypoint "/bin/bash" \
	bioraddbg/atac-seq-filter-beads

# 加了chrX,Y,M
echo "<=== 7.bead filtration ===>"
docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
bioraddbg/atac-seq-filter-beads \
-i /data/alignments_all/  \
-o /data/bead_filtration/ \
-r mm10

# 因为一直出错，因此该步已经生成了一部份文件，剩下文件单行代码生成
# 作者源代码储存在/mnt/d/scATAC/sci_reports2022中
# 在本地跑有问题，R包无法使用
Rscript /mnt/d/scATAC/sci_reports2022/callBeadKnee.R -o /mnt/d/scATAC/sci_reports2022/sequence/bead_filtration/whitelist/ -bt NULL -ft NULL /mnt/d/scATAC/sci_reports2022/sequence/bead_filtration/whitelist/alignments.possorted.tagged.filtered.barcodequants.csv
# 下文在docker主机里面跑

# $output_dir = /data/bead_filtration/
# $Name = alignments.possorted.tagged
# $bamFile = /data/alignments_all/alignments.possorted.tagged.bam
# $bamIndex= /data/alignments_all/alignments.possorted.tagged.bam.bai
# $basename /data/alignments_all/alignments.possorted.tagged.bam
# $bedGenome=/mnt/d/scATAC/sci_reports2022/chrom_mm10.sizes

cd /data/bead_filtration/
Rscript ../../callBeadKnee.R -o /data/bead_filtration/whitelist/ -bt NULL -ft NULL /data/bead_filtration/whitelist/alignments.possorted.tagged.filtered.barcodequants.csv

# 有问题，因此直接用R
library(kneecallr)
library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)
barcodeQuants <- "/data/bead_filtration/whitelist/alignments.possorted.tagged.filtered.barcodequants.csv"
fragmentThreshold <- NULL
beadsToPass <- NULL
outDir <- "/data/bead_filtration/whitelist"
df <- fread(barcodeQuants, col.names=c("sample", "count")) %>% arrange(-count)
knee <- kneecallr::callKnee(df)
saveRDS(knee, file=paste0(outDir,"/bead_kneeCall.rds"))

samplePlot <- kneecallr::plotSampleCounts(df,knee) + theme_classic() + ylab("Reads per barcode") + xlab("Barcodes in rank-descending order") +
        annotation_logticks() +
        scale_x(breaks=c(1,c(1,10) %o% 10^(1:nchar(trunc(nrow(df)))-1)), labels=c(1,c(1,10) %o% 10^(1:nchar(trunc(nrow(df)))-1))) +
        scale_y(breaks=c(1,c(1,10) %o% 10^(1:nchar(trunc(max(df$count))))), labels=c(1,c(1,10) %o% 10^(1:nchar(trunc(max(df$count))))))

densityPlot <- kneecallr::plotCountDensity(knee) + theme_classic() + ylab("Density") + xlab("Count per barcode (log10)")

kneePlot <- kneecallr::plotKnee(df,knee) + theme_classic() + ylab("Cumulative read count") + xlab("Barcode rank")

# write.table(knee$filtered_data %>% filter(count >= knee$threshold) %>% arrange(-count), paste0(outDir,"/aboveKneeBarcodes.csv"), sep=",", row.names=F, col.names=T, quote=F)

cp <- cowplot::plot_grid(densityPlot, kneePlot, samplePlot, nrow=1)
title <- ggdraw() + draw_label(paste("Bead knee threshold =", signif(knee$threshold, 3), "fragments /", knee$num_accepted_samples, "beads"), fontface='bold')
pg <- plot_grid(title, cp, ncol=1, rel_heights=c(0.1, 1))
ggsave(pg, filename=paste0(outDir,"/bead_kneeCurve.png"), width=24, height=8)
# R结束


cat /data/bead_filtration/whitelist/aboveKneeBarcodes.csv | cut -d, -f1 | tail -n +2 - > /data/bead_filtration/whitelist/bap_bead_whitelist.csv

mkdir -p /data/bead_filtration/bap_out/.internal/samples/
mkdir -p /data/bead_filtration/bap_out/
mkdir -p /data/bead_filtration/bap_out/final/
mkdir -p /data/bead_filtration/bap_out/temp/
mkdir -p /data/bead_filtration/bap_out/temp/filt_split/
mkdir -p /data/bead_filtration/bap_out/logs/
mkdir -p /data/bead_filtration/bap_out/temp/frag_overlap/
mkdir -p /data/bead_filtration/bap_out/temp/drop_barcode/
mkdir -p /data/bead_filtration/bap_out/knee/

python /bap/bap/bin/python/10_quantBarcode_Filt.py -i /data/alignments_all/alignments.possorted.tagged.bam --name alignments.possorted.tagged --output /data/bead_filtration/bap_out/temp/filt_split \
 --barcode-tag XB --min-fragments 1 --bedtools-genome /bap/bap/anno/bedtools/chrom_mm10.sizes --ncores 6 --mapq 30 --barcode-whitelist /data/bead_filtration/whitelist/bap_bead_whitelist.csv

#compute NC per read; can be RAM intensive so using semaphore
inbams=$(ls /data/bead_filtration/bap_out/temp/filt_split | grep raw.bam)
for bam in $inbams
do
        sem --bg -j 10 /usr/bin/Rscript /bap/bap/bin/R/11a_computeNCperRead.R /data/bead_filtration/bap_out/temp/filt_split/$bam XB
done
sem --wait

#not ram intensive, can run in embarrassingly parallel fashion
for bam in $inbams
do
        outbam=$(basename $bam .raw.bam).bam
        outtsv=$(basename $bam .bam)_ncRead.tsv
        python /bap/bap/bin/python/11b_annoFiltNC.py --input /data/bead_filtration/bap_out/temp/filt_split/$bam --output /data/bead_filtration/bap_out/temp/filt_split/$outbam \
        --nc-filt 6 --dict-file /data/bead_filtration/bap_out/temp/filt_split/$outtsv --bead-barcode XB &
done
wait

#compute fragment overlaps; can be RAM intensive so using semaphore
bams=$(find /data/bead_filtration/bap_out/temp/filt_split/ | grep .bam | grep -v raw)
for file in $bams
do
        sem --bg -j 10 /usr/bin/Rscript /bap/bap/bin/R/12_fragOverlapMetricsChr.R $file XB /data/bead_filtration/bap_out/final/alignments.possorted.tagged.barcodequants.csv /bap/bap/anno/blacklist/mm10.full.blacklist.bed
done
sem --wait

#get implicated barcode pairs
/usr/bin/Rscript /getImplicatedBarcodes.R /data/bead_filtration/bap_out/temp/frag_overlap/ /data/bead_filtration/bap_out/final/alignments.possorted.tagged.barcodequants.csv /data/bead_filtration/bap_out/final/alignments.possorted.tagged.implicatedBarcodes.csv /data/bead_filtration/whitelist/aboveKneeBarcodes.csv

mkdir -p /data/bead_filtration/bap_out/final/jaccard

#call knee
zcat /data/bead_filtration/bap_out/final/alignments.possorted.tagged.implicatedBarcodes.csv | cut -d, -f1,6 | tail -n +2 - > /data/bead_filtration/bap_out/final/alignments.possorted.tagged.temp.implicatedBarcodes.csv

#use kneecallr to find threshold
/usr/bin/Rscript /callJaccardKnee.R /data/bead_filtration/bap_out/final/alignments.possorted.tagged.temp.implicatedBarcodes.csv /data/bead_filtration/bap_out/final/jaccard
rm /data/bead_filtration/bap_out/final/alignments.possorted.tagged.temp.implicatedBarcodes.csv

#move around outputs
cp -r /data/bead_filtration/bap_out/final/jaccard/ /data/bead_filtration/jaccard/
rm -r /data/bead_filtration/bap_out

if [[ -f "/imageInfo.txt" ]]; then
        cp /imageInfo.txt /data/bead_filtration/.filter-beads-version.txt
fi

```

* At a high level, this Tool takes the bead barcodes that were deemed to contain DNA fragments from cells and
merges bead barcodes that “see the same cell to a droplet barcode.  


```bash
# bead deconvolution

echo "<=== 8.bead deconvolution ===>"
docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
bioraddbg/atac-seq-deconvolute \
-i /data/alignments_all/ \
-f /data/bead_filtration/ \
-r mm10 \
-o /data/deconvoluted_data/ 
```

```bash
# cell filtratron

echo "<=== 9.cell filtratron ===>"
docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
bioraddbg/atac-seq-cell-filter \
-i /data/deconvoluted_data/ \
-r mm10 \
-o /data/cells_filtered/ 
```
## CALL PEAK
```bash
# peak calling

echo "<=== 10.peak calling ===>"
docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
bioraddbg/atac-seq-macs2 \
-i /data/deconvoluted_data/ \
-r mm10 \
-o /data/peaks 

# alignments/中数据尝试，错误的步骤
docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
bioraddbg/atac-seq-macs2 \
-i /data/alignments/ \
-r mm10 \
-o /data/err_peaks/ 
```
```bash
# ATAC-seq QC

echo "<=== 11.ATAC-seq QC ===>"
docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
bioraddbg/atac-seq-qc \
-r mm10 \
-d /data/deconvoluted_data/ \
-p /data/peaks \
-o /data/atac_qc
```

## 建立矩阵
```bash
# count matrix

echo "<=== 12.cells-by-peaks count matrix ===>"
docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
bioraddbg/atac-seq-chromvar \
-r mm10 \
-d /data/cells_filtered/ \
-p /data/peaks \
-o /data/count_matrix
```
## report
```bash
# report

echo "<=== 13.report===>"
docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
bioraddbg/atac-seq-report \
-i /data/ \
-o /data/report
```



For generation of the fragments fle, which contain the start and end 
genomic coordinates of all aligned sequenced fragments, sorted bam fles were further process with “bap-frag” 
module of BAP (https://github.com/caleblareau/bap, v0.6.0).  


sci-ATAC-Seq数据预处理后得到的是二进制BAM文件，10X Genomics Chromium Single Cell ATAC数据预处理后的得到的是Fragments文本文档。10X Genomics Chromium Single Cell ATAC配套有“御用”预处理工具箱CellRanger，可以完成所有Read的预处理步骤。
