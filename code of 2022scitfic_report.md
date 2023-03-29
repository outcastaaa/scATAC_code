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
trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired -o ../trimmed/  R1.fastq.gz R2.fastq.gz  

docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
bioraddbg/atac-seq-trim-reads \
-i /data/debarcoded_reads -o /data/trimmed_reads


echo "<=== 4.quality control again ===>"
mkdir -p /mnt/d/scATAC/sci_reports2022/fastqc/again
cd /mnt/d/scATAC/sci_reports2022/trimmed
fastqc -t 8 -o ../fastqc/again   ./*.gz 
cd /mnt/d/scATAC/sci_reports2022/fastqc/again
multiqc .
```

```bash
# align
## 下载基因组
cd /mnt/d/scATAC/sci_reports2022/genome
for i in {1..19} 
do
  echo $i
  wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/chromosomes/chr$i.fa.gz
done
## 建索引
docker run --rm -it -v ~/genomes/hg19/:/genome/
--entrypoint "/bin/bash"
bioraddbg/atac-seq-bwa

mkdir -p /mnt/d/scATAC/sci_reports2022/genome/bwa_index
bwa_index=/mnt/d/scATAC/sci_reports2022/genome/bwa_index
bwa index -p $bwa_index/mm10 





echo "<=== 5.reads alignment ===>"
mkdir -p /mnt/d/scATAC/sci_reports2022/alignment
bowtie2_index=/mnt/d/scATAC/sci_reports2022/genome/mm10
align_dir=/mnt/d/scATAC/sci_reports2022/alignment

cd /mnt/d/scATAC/sci_reports2022/trimmed
bowtie2  -p 7 -x  $bowtie2_index --very-sensitive -X 2000 -1  1.fq.gz -2 2.fq.gz \
  2>$align_dir/sample.summary \
  -S $align_dir/sample.sam 

docker run --rm -v /mnt/d/scATAC/sci_reports2022/sequence/:/data/ \
-v /mnt/d/scATAC/sci_reports2022/genome/bwa_index/:/genome/ \
bioraddbg/atac-seq-trim-reads \
-i /data/trimmed_reads/ -o /data/alignments/ \
-r /genome/
```
* [基因组网站](https://hgdownload.soe.ucsc.edu/downloads.html)









For generation of the fragments fle, which contain the start and end 
genomic coordinates of all aligned sequenced fragments, sorted bam fles were further process with “bap-frag” 
module of BAP (https://github.com/caleblareau/bap, v0.6.0).  


sci-ATAC-Seq数据预处理后得到的是二进制BAM文件，10X Genomics Chromium Single Cell ATAC数据预处理后的得到的是Fragments文本文档。10X Genomics Chromium Single Cell ATAC配套有“御用”预处理工具箱CellRanger，可以完成所有Read的预处理步骤。
