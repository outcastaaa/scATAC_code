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
ln -s  ~/data/sra/SRR13783625_?.fastq.gz  /mnt/d/scATAC/sci_reports2022/sequence/
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

echo "<=== quality control ===>"
mkdir -p /mnt/d/scATAC/sci_reports2022/fastqc
fastqc -t 8 -o /mnt/d/scATAC/sci_reports2022/fastqc /mnt/d/scATAC/sci_reports2022/sequence/*.gz 
cd /mnt/d/scATAC/sci_reports2022/fastqc 
multiqc .


echo "<=== debarcoding ===>"
cd ~/data/sra
source venv3/bin/activate
mkdir -p /mnt/d/scATAC/sci_reports2022/debarcode
cd /mnt/d/scATAC/sci_reports2022/debarcode/
bap-barcode v2.1 -a /mnt/d/scATAC/sci_reports2022/sequence/SRR13783625_1.fastq.gz -b /mnt/d/scATAC/sci_reports2022/sequence/SRR13783625_2.fastq.gz --nmismatches 1 


echo "<=== read trimming ===>"
mkdir -p /mnt/d/scATAC/sci_reports2022/trimmed
cd mkdir -p /mnt/d/scATAC/sci_reports2022/debarcode
trim_galore --phred33 --length 35 -e 0.1 --stringency 3 --paired -o /mnt/d/scATAC/sci_reports2022/trimmed/  R1.fastq.gz R2.fastq.gz  


echo "<=== quality control again ===>"
mkdir -p /mnt/d/scATAC/sci_reports2022/fastqc/again
fastqc -t 8 -o /mnt/d/scATAC/sci_reports2022/fastqc/again /mnt/d/scATAC/sci_reports2022/trimmed/*.gz 
cd /mnt/d/scATAC/sci_reports2022/fastqc/again
multiqc .
```

For generation of the fragments fle, which contain the start and end 
genomic coordinates of all aligned sequenced fragments, sorted bam fles were further process with “bap-frag” 
module of BAP (https://github.com/caleblareau/bap, v0.6.0).  


sci-ATAC-Seq数据预处理后得到的是二进制BAM文件，10X Genomics Chromium Single Cell ATAC数据预处理后的得到的是Fragments文本文档。10X Genomics Chromium Single Cell ATAC配套有“御用”预处理工具箱CellRanger，可以完成所有Read的预处理步骤。
