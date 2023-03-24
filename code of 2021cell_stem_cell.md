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


Upstream analysis of scATAC-Seq data
We performed the upstream analysis using the samtools (Li et al., 2009) (v1.9), bedtools (Quinlan, 2014) (v2.27.1), Picard tools (v2.9.0;
http://broadinstitute.github.io/picard/) and BWA (Li and Durbin, 2009) (v0.7.17). First, we aligned fastq files to the GRCh38 reference 
genome (average 473,886 reads per cell), followed by marking duplicates with MarkDuplicates function from Picard tools and
removing duplicates using samtools view with -F 1804 parameter per each cell. Overall with average duplicates rate 77% we obtained 91,554 reads per cell after removing duplicates. Next, we transformed bam files to bed files using bamtobed bedtools function
in bedpe mode and kept only fragments that are no bigger than 1000 bp using a custom script. We called peaks (for the clusters with
more than 50 cells) using the SnapATAC approach (Fang et al., 2019) with macs2 (Zhang et al., 2008) parameters ‘‘–nomodel–shift
100–ext 200–qval 5e-2 -B’’ and obtained 152,283 peaks. Importantly, for the downstream analysis in R, we binarized counts per cell
using Signac (Stuart et al., 2019; https://github.com/timoast/signac/) BinarizeCounts function, resulting in 32,217 fragments per cell
on average.

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

### SnapATAC.R.Rmd

```r
# title: "SnapATAC_merged"
# output: html_notebook

library(SnapATAC)
library(viridisLite);
library(ggplot2);

x.sp = createSnap(
    file="merged.snap",
    sample="merged234",
    num.cores=4
  );
# at this stage might be good to remove cells with low quality 

showBinSizes("merged.snap");

x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=4);

x.sp = makeBinary(x.sp, mat="bmat")


library(GenomicRanges);
black_list = read.table("hg38_blacklist.bed");
black_list.gr = GRanges(
    black_list[,1], 
    IRanges(black_list[,2], black_list[,3])
  );
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr));
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
x.sp

chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};
 x.sp

bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
hist(
    bin.cov[bin.cov > 0], 
    xlab="log10(bin cov)", 
    main="log10(Bin Cov)", 
    col="lightblue", 
    xlim=c(0, 5)
  );
 bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
 idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
 x.sp = x.sp[, idy, mat="bmat"];
 x.sp

x.sp = runDiffusionMaps(
    obj=x.sp,
    input.mat="bmat", 
    num.eigs=50
  );

plotDimReductPW(
    obj=x.sp, 
    eigs.dims=1:50,
    point.size=0.3,
    point.color="grey",
    point.shape=19,
    point.alpha=0.6,
    down.sample=5000,
    pdf.file.name=NULL, 
    pdf.height=7, 
    pdf.width=7
  );

# x.sp@metaData = cbind(x.sp@metaData, coldata.final)

x.sp = runKNN(
    obj=x.after.sp,
    eigs.dims= 1:50,
    k=10
  );

x.sprunCluster(
    obj=x.after.sp,
    tmp.folder=tempdir(),
    louvain.lib="R-igraph",
    seed.use=10
  );

x.sp@metaData$cluster = x.sp@cluster;

x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=2:50, 
    method="umap",
    seed.use=10
  );

 plotViz(
    obj=x.sp,
    method="umap",
    point.color=x.sp@sample, 
    point.size=0.1, 
    text.add= FALSE,
    down.sample=10000,
    legend.add=TRUE
  );



# calculate the ensemble signals for each cluster
ensemble.ls = lapply(split(seq(length(x.sp@cluster)), x.sp@cluster), function(x){
	SnapATAC::colMeans(x.sp[x,], mat="bmat");
	})
# cluster using 1-cor as distance  
hc = hclust(as.dist(1 - cor(t(do.call(rbind, ensemble.ls)))), method="ward.D2");
plotViz(
    obj=x.sp,
    method="umap", 
    main="batch1",
    point.color=x.sp@cluster, 
    point.size=1, 
    point.shape=19, 
    point.alpha=0.8, 
    text.add=TRUE,
    text.size=1.5,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    down.sample=10000,
    legend.add=FALSE
    );
plot(hc, hang=-1, xlab="");


clusters.sel = names(table(x.sp@cluster))[which(table(x.sp@cluster) > 50)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
    print(clusters.sel[i]);
    runMACS(
        obj=x.sp[which(x.sp@cluster==clusters.sel[i]),], 
        output.prefix=paste0("atac_v1_adult_brain_fresh_5k.", gsub(" ", "_", clusters.sel)[i]),
        path.to.snaptools="/home/berest/miniconda2/bin/snaptools",
        path.to.macs="/home/berest/miniconda2/bin/macs2",
        gsize="hs", # mm, hs, etc
        buffer.size=500, 
        num.cores=1,
        macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B ",
        tmp.folder=tempdir()
   );
 }, mc.cores=5);
# assuming all .narrowPeak files in the current folder are generated from the clusters
peaks.names = system("ls | grep narrowPeak", intern=TRUE);
peak.gr.ls = lapply(peaks.names, function(x){
    peak.df = read.table(x)
    GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
  })
peak.gr = GenomicRanges::reduce(Reduce(c, peak.gr.ls));
peak.gr
peak.df = as.data.frame(peak.gr)


peak.df$seqnames =  substr(substring(peak.df$seqnames,3),1,nchar(substring(peak.df$seqnames,3))-1) 

allowchr = c(paste0("chr",seq(1:23)),"chrY","chrX")
peak.filt.df = peak.df[which(peak.df$seqnames %in% allowchr),]

library(tidyverse)
write_tsv(peak.filt.df, path = "snapatac.clust.peaks.bed",col_names = F)
```


### CreateSeuratObjandQCplots.Rmd
```r
# title: "Create Seurat object + QC plots"
# output: html_notebook

set.seed(2017)

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
* Load consensus peakset  

```r
# get the couns matrix 
peakfile <- "/g/scb2/zaugg/berest/Projects/collaborations/SCRNAATAC/scratch/CVEJIC/test/snapatac.clust.peaks.bed"

peaks.df = as.data.frame(read_tsv(peakfile, col_names = F))
peaks.df$ID = paste0(peaks.df$X1,":",peaks.df$X2,"-",peaks.df$X3)
rownames(peaks.df) = peaks.df$ID

peaks <- getPeaks(peakfile, sort_peaks = TRUE)

# genome coordinates for binning 
genome <- BSgenome.Hsapiens.UCSC.hg38
seq.gen  = seqlengths(genome)


# external data for TSS enrichment
# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86)
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]

tss.ranges <- GRanges(
  seqnames = seqnames(gene.ranges),
  ranges = IRanges(start = start(gene.ranges), width = 2),
  strand = strand(gene.ranges)
)

seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')
```
* Load metadata
[数据所在](https://gitlab.com/cvejic-group/integrative-scrna-scatac-human-foetal/-/tree/master/Data/scATAC_CSV_file_for_Scanpy)  

```r
# load first metadata
coldata = as.data.frame(read_tsv("/datadisk/Desktop/Projects/SCRNATAC/scATAC/meta.lane234.txt"))
# load statistics from alignment and removing duplicates
stats.df = read_tsv("/datadisk/Desktop/Projects/SCRNATAC/scATAC/stats.lane234.txt")

# merge into on object
coldata.x = merge(coldata, stats.df, by = "ID")
# load experimental metadata
lane.meta = read_csv("/datadisk/Desktop/Projects/SCRNATAC/scATAC/MetaData_scATAC-Seq/Lane234.csv")
# final merging 
coldata.final = merge(coldata.x, lane.meta, by.x = "ID", by.y = "SANGER SAMPLE ID")

rm(coldata,coldata.x,stats.df,lane.meta)
rownames(coldata.final) = paste0(coldata.final$ID,".clean.bam")

# split the column 
coldata.final$origin = str_split(coldata.final$`SUPPLIER SAMPLE NAME`, " ", 4, simplify = TRUE)[,4]
```

* Create Signac object  
```r
# load or create peak and bin matrices 
fragment.path =  "/g/scb2/zaugg/berest/Projects/collaborations/SCRNAATAC/scratch/CVEJIC/test/merged.l234.sort.CO.bed.gz"


peak_matrix_name = "/g/scb2/zaugg/berest/Projects/collaborations/SCRNAATAC/scratch/CVEJIC/test/merged234.SNAPpeaksnotfilt.Rds"

if (file.exists(peak_matrix_name)) {
    peak_matrix = readRDS(file = peak_matrix_name)
} else {
    peak_matrix <- FeatureMatrix(fragments = fragment.path, features = peaks)
    saveRDS(peak_matrix, file = peak_matrix_name)
  }

bin_matrix_name = "/g/scb2/zaugg/berest/Projects/collaborations/SCRNAATAC/scratch/CVEJIC/test/merged234.SNAPpeaksnotfilt.BM.Rds"

if (file.exists(bin_matrix_name)) {
    bin_matrix = readRDS(file = bin_matrix_name)
} else {
    bin_matrix <- FeatureMatrix(fragments = fragment.path, genome = seq.gen,
                                 binsize = 5000, chunk = 100)
    saveRDS(bin_matrix, file = bin_matrix_name)
  }


# Create object 

merged <- CreateSeuratObject(
  counts = peak_matrix,
  assay = 'peaks',
  project = 'merged234',
  min.cells = 1,
  meta.data = coldata.final
)

# add bins matrix 
merged[["bins"]] = CreateAssayObject(counts = bin_matrix)
# add fragments file to the object 
merged <- SetFragments(
  object = merged,
  file = fragment.path
)

# Binarize bith bins and peaks matrix 
merged = BinarizeCounts(object = merged, assay = c("peaks","bins"))

rm(peak_matrix, bin_matrix)
```

* Quality checks  

```r
### QC

merged <- FRiP(
  object = merged,
  peak.assay = "peaks",
  bin.assay = "bins"
)

merged@meta.data$blacklist_ratio <- FractionCountsInRegion(
  object = merged,
  assay = 'peaks',
  regions = blacklist_hg38
)


merged <- NucleosomeSignal(object = merged)

merged$nucleosome_group <- ifelse(merged$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
p1 = FragmentHistogram(object = merged)

p1


p2 = VlnPlot(
  object = merged,
  features = c('FRiP', 'blacklist_ratio', 'nucleosome_signal', 'nFeature_peaks'),
  pt.size = 0.1,
  ncol = 4) + NoLegend()

p2 


# to save time use the first 2000 TSSs
merged <- TSSEnrichment(object = merged, tss.positions = tss.ranges[1:2000]) # 

merged$high.tss <- ifelse(merged$TSS.enrichment > 2, 'High', 'Low')
p3 = TSSPlot(merged, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend()

p3
```

* QC plots for paper  
```r
# 
library(ggthemes)
library(viridis)
library(reshape2)
theme_ms <- function(base_size=12, base_family="Helvetica") {
  library(grid)
  (theme_bw(base_size = base_size, base_family = base_family)+
      theme(text=element_text(color="black"),
            axis.title=element_text(face="bold", size = rel(1.3)),
            axis.text=element_text(size = rel(1), color = "black"),
            legend.title=element_text(face="bold"),
            legend.text=element_text(face="bold"),
            legend.background=element_rect(fill="transparent"),
            legend.key.size = unit(0.8, 'lines'),
            panel.border=element_rect(color="black",size=1),
            panel.grid=element_blank()
    ))
}

pFL = ggplot() + geom_histogram(aes( x = p1$data$length, y = (..count..)/sum(..count..)), binwidth = 5, colour = "black", fill = "black") + 
  xlab("Length of the fragments, bp") + ylab("Percent of Fragments") + theme_ms()
pFL
# ggsave("/datadisk/Desktop/Projects/SCRNATAC/scATAC/plots/QC/FL.pdf", plot = pFL, width = 6, height = 3)


pHEX = ggplot(merged@meta.data, aes(x = log10(nFeature_peaks), y = TSS.enrichment)) +
  geom_hex(bins = 100) +
  theme_bw() + scale_fill_viridis() +
  xlab("log10 Fragments in peaks") +
  ylab("TSS Enrichment") +
  ylim(0,10) + 
  geom_hline(yintercept = 2, lty = "dashed") +
  geom_vline(xintercept = log10(1000), lty = "dashed") +
  ggtitle(paste0("Total: ", nrow(merged@meta.data)," Passed: 3611")) + 
  theme_ms()
pHEX
# ggsave("/datadisk/Desktop/Projects/SCRNATAC/scATAC/plots/QC/TSSvsFIP.pdf", plot = pHEX, width = 5, height = 4)

data1 = merged@meta.data[,c("initial")]
data1 = melt(data1)
data1$variable = "initial"
pI = ggplot(data1, aes(x = variable, y = log10(value))) + geom_jitter(alpha = 0.2) + 
  geom_violin(colour = "darkred", fill = "darkred",alpha = 0.7)  +  
  xlab("") + ylab("log10 Initial \n amount of reads") + theme_ms() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())

data2 = merged@meta.data[,c("rmDup")]
data2 = melt(data2)
data2$variable = "initial"
pR = ggplot(data2, aes(x = variable, y = log10(value))) + geom_jitter(alpha = 0.2) + 
  geom_violin(colour = "darkred", fill = "darkred",alpha = 0.7)  +  
  xlab("") + ylab("log10 Reads after \n duplicates removal") + theme_ms() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())

data3 = merged@meta.data[,c("duplication_rate")]
data3 = melt(data3)
data3$variable = "initial"
pD = ggplot(data3, aes(x = variable, y = value*100)) + geom_jitter(alpha = 0.2) + 
  geom_violin(colour = "darkred", fill = "darkred",alpha = 0.7)  +  
  xlab("") + ylab("Duplication rate, %") + theme_ms() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())

data4 = merged@meta.data[,c("nCount_bins")]
data4 = melt(data4)
data4$variable = "initial"
pF = ggplot(data4, aes(x = variable, y = log10(value))) + geom_jitter(alpha = 0.2) + 
  geom_violin(colour = "darkred", fill = "darkred",alpha = 0.7)  +  
  xlab("") + ylab("log10 Fragments per cell")+ theme_ms()  + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                      axis.ticks.x=element_blank()) 

data5 = merged@meta.data[,c("FRiP")]
data5 = melt(data5)
data5$variable = "initial"
pFr = ggplot(data5, aes(x = variable, y = value * 100)) + geom_jitter(alpha = 0.2) + 
  geom_violin(colour = "darkred", fill = "darkred",alpha = 0.7)  +  
  xlab("") + ylab("Fragments in peaks, %") + theme_ms() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())

data6 = merged@meta.data[,c("blacklist_ratio")]
data6 = melt(data6)
data6$variable = "initial"
pB = ggplot(data6, aes(x = variable, y = value * 100)) + geom_jitter(alpha = 0.2) + 
  geom_violin(colour = "darkred", fill = "darkred",alpha = 0.7)  +  
  xlab("") + ylab("Reads in \n blacklisted regions, %") + theme_ms() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())

data7 = merged@meta.data[,c("nucleosome_signal")]
data7 = melt(data7)
data7$variable = "initial"
pN = ggplot(data7, aes(x = variable, y = value )) + geom_jitter(alpha = 0.2) + 
  geom_violin(colour = "darkred", fill = "darkred",alpha = 0.7)  +  
  xlab("") + ylab("Nuclesome signal") + theme_ms() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())

plotC = (pI|pR | pD) / (pF | pFr | pB | pN)

plotC

# ggsave("/datadisk/Desktop/Projects/SCRNATAC/scATAC/plots/QC/QCboxes.pdf", plot = plotC, width = 8, height = 6)

rm(data1,data2,data3,data4,data5,data6,data7,
   pI,pR,pD,pF ,pFr ,pB, pN)


data = merged@meta.data[,c("initial","rmDup","duplication_rate","nCount_bins",
                           "FRiP","blacklist_ratio","nucleosome_signal","TSS.enrichment")]

colnames(data) = c("Initial_reads","Reads_after_rmDUP","duplication_rate",
                   "Fragments_cell","FRIP","Reads_blacklist","Nucleosome_signal","TSS_enrichment")

# write_tsv(data, path = "/datadisk/Desktop/Projects/SCRNATAC/scATAC/plots/QC/QCdata.tsv", col_names = T)
data.add = p1$data[,c("chr","start","end","length")]
# write_tsv(data.add, path = "/datadisk/Desktop/Projects/SCRNATAC/scATAC/plots/QC/QCFragmentLengthdistr.tsv", col_names = T)

rm(data, data.add)
```  

* Filter merged object  
```r
merged.filt = subset(merged, subset = nFeature_peaks > 1000 & 
                        # nFeature_peaks < 20000 & 
                        FRiP > 0.5 & 
                        blacklist_ratio < 0.05 & 
                        nucleosome_signal < 10 & 
                        TSS.enrichment > 2)
merged.filt

table(merged.filt$origin)
table(merged.filt$batch)
table(merged.filt$SAMPLE)

saveRDS(merged.filt, file = "/datadisk/Desktop/Projects/SCRNATAC/scATAC/files/merged234.Snap.filt.Rds")
```  

!!! 此处即得到文件：merged234.Snap.filt.Rds (983.9 MB)    



### IndividualSamples.Rmd  该步骤看情况

```r
# title: "scATAC for individual samples"
# output: html_notebook

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

merged.filt = readRDS( file = "/datadisk/Desktop/Projects/SCRNATAC/scATAC/files/merged234.Snap.filt.Rds")
```

* sample 2 
```r
sample2 = subset(merged.filt, subset = SAMPLE == "Sample2")

pca.dim = 1:50
sample2 <- sample2 %>% 
               RunTFIDF( method = 2) %>%
               FindTopFeatures( min.cutoff = 'q0') %>%
               RunSVD(
                assay = 'peaks',
                reduction.key = 'LSI_',
                reduction.name = 'lsi') %>%
               RunUMAP( reduction = 'lsi', dims = pca.dim, umap.method = "uwot",
                        a = 1, b = 0.9, n.neighbors = 10) %>%
               FindNeighbors( reduction = 'lsi', dims = pca.dim,annoy.metric = "cosine") %>%
               FindClusters( resolution = 1,verbose = FALSE)

ggplot() +
  geom_point(aes(y = sample2@reductions$lsi@stdev / sum(sample2@reductions$lsi@stdev), x = 1:50)) + 
  xlab("LSI dimensions") + ylab("Variance explaned")

ggplot() + geom_point(aes(x = pca.dim, y = abs(sapply(pca.dim, function(i) cor(sample2@meta.data$nCount_peaks, sample2@reductions$lsi@cell.embeddings[,i]))))) +
  xlab("LSI dimensions") + ylab("Correlation fragments in peaks with lsi dimension")

ggplot() + geom_point(aes(x = pca.dim, y= abs(sapply(pca.dim, function(i) cor(sample2@meta.data$duplication_rate, sample2@reductions$lsi@cell.embeddings[,i]))))) +
  xlab("LSI dimensions") + ylab("Correlation duplication rate with lsi dimension")



DimPlot(object = sample2, label = TRUE) + NoLegend()

DimPlot(object = sample2, label = F, group.by = "origin")

DimPlot(object = sample2, label = F, group.by = "SAMPLE")

DimPlot(object = sample2, label = F, group.by = "batch")
```
* sample 3  
```r
sample3 = subset(merged.filt, subset = SAMPLE == "Sample3")

pca.dim = 1:50
sample3 <- sample3 %>% 
               RunTFIDF( method = 2) %>%
               FindTopFeatures( min.cutoff = 'q0') %>%
               RunSVD(
                assay = 'peaks',
                reduction.key = 'LSI_',
                reduction.name = 'lsi') %>%
               RunUMAP( reduction = 'lsi', dims = pca.dim, umap.method = "uwot",
                        a = 1, b = 0.9, n.neighbors = 10) %>%
               FindNeighbors( reduction = 'lsi', dims = pca.dim,annoy.metric = "cosine") %>%
               FindClusters( resolution = 1,verbose = FALSE)

ggplot() +
  geom_point(aes(y = sample3@reductions$lsi@stdev / sum(sample3@reductions$lsi@stdev), x = 1:50)) + 
  xlab("LSI dimensions") + ylab("Variance explaned")

ggplot() + geom_point(aes(x = pca.dim, y = abs(sapply(pca.dim, function(i) cor(sample3@meta.data$nCount_peaks, sample3@reductions$lsi@cell.embeddings[,i]))))) +
  xlab("LSI dimensions") + ylab("Correlation fragments in peaks with lsi dimension")

ggplot() + geom_point(aes(x = pca.dim, y= abs(sapply(pca.dim, function(i) cor(sample3@meta.data$duplication_rate, sample3@reductions$lsi@cell.embeddings[,i]))))) +
  xlab("LSI dimensions") + ylab("Correlation duplication rate with lsi dimension")




DimPlot(object = sample3, label = TRUE) + NoLegend()

DimPlot(object = sample3, label = F, group.by = "origin")

DimPlot(object = sample3, label = F, group.by = "SAMPLE")

DimPlot(object = sample3, label = F, group.by = "batch")
```

* sample 4  
```r
sample4 = subset(merged.filt, subset = SAMPLE == "Sample4")

pca.dim = 1:50
sample4 <- sample4 %>% 
               RunTFIDF( method = 2) %>%
               FindTopFeatures( min.cutoff = 'q0') %>%
               RunSVD(
                assay = 'peaks',
                reduction.key = 'LSI_',
                reduction.name = 'lsi') %>%
               RunUMAP( reduction = 'lsi', dims = pca.dim, umap.method = "uwot",
                        a = 1, b = 0.9, n.neighbors = 10) %>%
               FindNeighbors( reduction = 'lsi', dims = pca.dim,annoy.metric = "cosine") 


    sample4 = FindClusters(sample4, resolution = 1,verbose = TRUE)

ggplot() +
  geom_point(aes(y = sample4@reductions$lsi@stdev / sum(sample4@reductions$lsi@stdev), x = 1:50)) + 
  xlab("LSI dimensions") + ylab("Variance explaned")

ggplot() + geom_point(aes(x = pca.dim, y = abs(sapply(pca.dim, function(i) cor(sample4@meta.data$nCount_peaks, sample4@reductions$lsi@cell.embeddings[,i]))))) +
  xlab("LSI dimensions") + ylab("Correlation fragments in peaks with lsi dimension")

ggplot() + geom_point(aes(x = pca.dim, y= abs(sapply(pca.dim, function(i) cor(sample4@meta.data$duplication_rate, sample4@reductions$lsi@cell.embeddings[,i]))))) +
  xlab("LSI dimensions") + ylab("Correlation duplication rate with lsi dimension")




DimPlot(object = sample4, label = TRUE) + NoLegend()

DimPlot(object = sample4, label = F, group.by = "origin")

DimPlot(object = sample4, label = F, group.by = "SAMPLE")

DimPlot(object = sample4, label = F, group.by = "batch")
```

### MergedSamplesHarmony.Rmd  

```r
# title: "Merged ATAC-seq analysis + Harmony"
# output: html_notebook

library(Matrix)
library(tidyverse)
library(chromVAR)
set.seed(1978)
library(data.table)
library(Seurat)
library(Signac)
library(harmony)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)

merged.filt = readRDS( file = "/datadisk/Desktop/Projects/SCRNATAC/scATAC/files/merged234.Snap.filt.Rds")
```
* Continue with merged data  
```r
set.seed(2017)
merged.filt <- merged.filt %>% 
               RunTFIDF( method = 2) %>%
               FindTopFeatures( min.cutoff = 'q0') %>%
               RunSVD(
                assay = 'peaks',
                reduction.key = 'LSI_',
                reduction.name = 'lsi') 

FeatureScatter(merged.filt, 'LSI_1', 'nCount_peaks')
FeatureScatter(merged.filt, 'LSI_2', 'nCount_peaks')
pca.dim = 2:50
merged.filt = merged.filt %>%
               RunUMAP( reduction = 'lsi', dims = pca.dim, umap.method = "uwot",
                        a = 1, b = 0.9, n.neighbors = 10) %>%
               FindNeighbors( reduction = 'lsi', dims = pca.dim,annoy.metric = "cosine") %>%
               FindClusters( resolution = 1,verbose = FALSE)

ggplot() +
  geom_point(aes(y = merged.filt@reductions$lsi@stdev / sum(merged.filt@reductions$lsi@stdev), x = 1:50)) + 
  xlab("LSI dimensions") + ylab("Variance explaned")

ggplot() + geom_point(aes(x = pca.dim, y = abs(sapply(pca.dim, function(i) cor(merged.filt@meta.data$nCount_peaks, merged.filt@reductions$lsi@cell.embeddings[,i]))))) +
  xlab("LSI dimensions") + ylab("Correlation fragments in peaks with lsi dimension")

ggplot() + geom_point(aes(x = pca.dim, y= abs(sapply(pca.dim, function(i) cor(merged.filt@meta.data$duplication_rate, merged.filt@reductions$lsi@cell.embeddings[,i]))))) +
  xlab("LSI dimensions") + ylab("Correlation duplication rate with lsi dimension")




DimPlot(object = merged.filt, label = TRUE) + NoLegend()

DimPlot(object = merged.filt, label = F, group.by = "origin")

DimPlot(object = merged.filt, label = F, group.by = "SAMPLE")

DimPlot(object = merged.filt, label = F, group.by = "batch")
```

* Harmony  
```r

library(harmony)
set.seed(2017)
pca.dim.cor = 2:50


test <- harmony::RunHarmony(
  object = merged.filt,
  group.by.vars = c("SAMPLE","batch","origin"),
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE,
  dims.use = pca.dim.cor,
  max.iter.harmony = 20,
   plot_convergence = T,
  # tau = 5,
  max.iter.cluster = 200,
  sigma = 0.25,
  theta = c(2,4,4)
)




# re-compute the UMAP using corrected LSI embeddings
test <- test %>% RunUMAP(dims = 1:50, reduction = 'harmony',
                          umap.method = "uwot", n.neighbors = 10, n.components = 2
                        # a = 1, b = 0.9
                        ) %>%
                 FindNeighbors( reduction = 'harmony', dims = 1:50,annoy.metric = "cosine") %>%
                 FindClusters( resolution = 0.45, verbose = FALSE)


DimPlot(test, pt.size = 0.1) + ggplot2::ggtitle("Harmony integration")

DimPlot(test, pt.size = 0.1, group.by = "SAMPLE") 
DimPlot(test, pt.size = 0.1, group.by = "origin") 
DimPlot(test, pt.size = 0.1, group.by = "batch") 

saveRDS(test,  file = "/datadisk/Desktop/Projects/SCRNATAC/scATAC/files/test.merged234.Snap.HA.Rds")
```    

### chromVAR.Rmd
```r
# title: "chromVAR analysis"
# output: html_notebook

library(Matrix)
library(tidyverse)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2017)
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

* Load peaks and Seurat object  

```r
# overlap all peaks 
# get the couns matrix 
peakfile <- "/datadisk/Desktop/Projects/SCRNATAC/scATAC/peaks/snapatac.clust.peaks.bed"

peaks.df = as.data.frame(read_tsv(peakfile, col_names = F))
peaks.df$ID = paste0(peaks.df$X1,":",peaks.df$X2,"-",peaks.df$X3)
rownames(peaks.df) = peaks.df$ID
# peakset = "consensus"
peaks <- getPeaks(peakfile, sort_peaks = TRUE)



# genome coordinates for binning 
genome <- BSgenome.Hsapiens.UCSC.hg38
seq.gen  = seqlengths(genome)

 test = readRDS(file = "/datadisk/Desktop/Projects/SCRNATAC/scATAC/files/test.merged234.Snap.HA.Rds")
```
* Run ChromVAR  
```r

DefaultAssay(test) = "peaks"
# run chromVAR 

# prepare motifs 
 motifs_name = "/datadisk/Desktop/Projects/SCRNATAC/scATAC/HOCO11HUMAN.pfmlist.rds"
if (file.exists(motifs_name)) {
  motifs <- readRDS(motifs_name)
} else {
  library(universalmotif)
  library(TFBSTools)
  library(rlist)
  # data taken from https://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_pwm_HUMAN_mono.tar.gz
  # for file in pwm/*pwm; do echo "" >> "$file"; done
  pathToPWMs = "/scratch/berest/SC/CVEJIC/pwm/"
  pwmfiles = list.files(pathToPWMs)

  empty.list <- list()

  for (x in 1:length(pwmfiles)) {
  
    test1 = read_matrix(paste0(pathToPWMs,pwmfiles[x]), 
                      headers = ">", sep = "", positions = "rows")
    lel = convert_motifs(test1, class = "TFBSTools-PFMatrix")
  
    empty.list <- c(empty.list, lel)
  
  }

  pfm.list <- do.call(PFMatrixList, empty.list)

  list.save(pfm.list, file = motifs_name)
  
  motifs <- readRDS(motifs_name)
}

 
 

 # # Scan the DNA sequence of each peak for the presence of each motif
 
 motif.matrix_name = "/g/scb2/zaugg/berest/Projects/collaborations/SCRNAATAC/scratch/CVEJIC/test/motifmatrix.merged234snap.Rds"
 
 if (file.exists(motif.matrix_name)) {
   motif.matrix = readRDS(file  = motif.matrix_name)
 } else {
   
    motif.matrix <- CreateMotifMatrix(
      features = StringToGRanges(rownames(test), sep = c("-", "-")),
      pwm = motifs,
      genome = 'hg38',
      sep = c("-", "-"),
      use.counts = FALSE
    )
    colnames(motif.matrix) = TFBSTools::name(motifs)
    saveRDS(motif.matrix, file  = motif.matrix_name)

    motif.matrix = readRDS(file  = motif.matrix_name)

 }

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = motifs
)

# Add the Motif object to the assay
test[['peaks']] <- AddMotifObject(
  object = test[['peaks']],
  motif.object = motif
)



test <- RegionStats(
  object = test,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  sep = c("-", "-")
)


## actual running chromVar 
test <- RunChromVAR(
  object = test,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "peaks"
)

library(ggrepel)
library(rlist)
library(plotly)
deviations = as.matrix(GetAssayData(test, assay = "chromvar"))

deviations.s = deviations[order(rowSds(deviations), decreasing = T),]


nx = 20
ggplot() + geom_point(aes(y = rowSds(deviations.s), x = 1:nrow(deviations.s), label = rownames(deviations.s))) + 
  geom_label_repel(aes(y = rowSds(deviations.s[1:nx,]), x = 1:nrow(deviations.s[1:nx,]), 
                       label = unlist(str_split(rownames(deviations.s)[1:nx],"-",simplify =  F) %>% list.map(.[1])) ))


plot_ly(y = rowSds(deviations.s), x = 1:nrow(deviations.s),color = "black", name = unlist(str_split(rownames(deviations.s),"-",simplify =  F) %>% list.map(.[1]))
) %>% layout(showlegend = FALSE)

 deviations.s.plot = deviations.s[rowSds(deviations.s) > 1.2,]
 nrow(deviations.s.plot)

 saveRDS(test, file = "/datadisk/Desktop/Projects/SCRNATAC/scATAC/files/test.merged234.Snap.HA.ChromVar.Rds")
```


### IntegrateATACRNA3kb.Rmd

```r
# title: "Integrate scRNA and scATAC"
# output: html_notebook

library(Matrix)
library(tidyverse)
library(chromVAR)
set.seed(1978)
library(data.table)
library(Seurat)
library(Signac)
library(harmony)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(pals)
library(ggrepel)
library(patchwork)
```
* Load scRNA object  
```r
load("/datadisk/Desktop/Projects/SCRNATAC/scRNA/ScanpyToSeurat-20200329T150351Z-001/ScanpyToSeurat/RNA_all_samples.robj")
```
* load ATAC data    
```r
test = readRDS("/datadisk/Desktop/Projects/SCRNATAC/scATAC/files/test.merged234.Snap.HA.ChromVar.Rds")

peakfile <- "/datadisk/Desktop/Projects/SCRNATAC/scATAC/peaks/snapatac.clust.peaks.bed"

peaks.df = as.data.frame(read_tsv(peakfile, col_names = F))
peaks.df$ID = paste0(peaks.df$X1,":",peaks.df$X2,"-",peaks.df$X3)
rownames(peaks.df) = peaks.df$ID

peakset = "consensus"
peaks <- getPeaks(peakfile, sort_peaks = TRUE)


# size of extension upstream and downstream in kb 

extensions = c(3)
 # extract gene coordinates from Ensembl, and ensure name formatting is consistent with Seurat object 
gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')


for (ext in extensions) {
  
  
 
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 1000*ext, downstream = 1000*ext)


if (file.exists(paste0("/datadisk/Desktop/Projects/SCRNATAC/scATAC/geneact/GENEact.",ext,
                                       "kdownup.merged2l234.rds"))) {
  
   gene.activities = readRDS(file = paste0("/datadisk/Desktop/Projects/SCRNATAC/scATAC/geneact/GENEact.",ext,
                                         "kdownup.merged2l234.rds"))
} else {
  gene.activities <- FeatureMatrix(
    fragments = "/datadisk/Desktop/Projects/SCRNATAC/scATAC/fragments/merged.l234.sort.CO.bed.gz",
    features = genebodyandpromoter.coords,
    chunk = 20,
    cells = rownames(test@meta.data)
  )
  saveRDS(gene.activities, file = paste0("/datadisk/Desktop/Projects/SCRNATAC/scATAC/geneact/GENEact.",ext,
                                       "kdownup.merged2l234.rds"))


  gene.activities = readRDS(file = paste0("/datadisk/Desktop/Projects/SCRNATAC/scATAC/geneact/GENEact.",ext,
                                         "kdownup.merged2l234.rds"))

  
  }


 gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)


rownames(gene.activities) <- gene.key[rownames(gene.activities)]

gene.activities = gene.activities[unique(rownames(gene.activities)),]

# add the gene activity matrix to the Seurat object as a new assay, and normalize it
test[[paste0("RNA_geneprom_",ext,"kb")]] <- CreateAssayObject(counts = gene.activities)
# test = BinarizeCounts(object = test, assay = c("RNA"))

test <- NormalizeData(
  object = test,
  assay = paste0("RNA_geneprom_",ext,"kb"),
  normalization.method = 'LogNormalize',
  scale.factor = median(test@meta.data[,paste0("nFeature_RNA_geneprom_",ext,"kb")])
)

rm(gene.activities)
  
  
}



RNAcatsign = names(table(merged@meta.data[which(merged@meta.data$gate == "38-"),"annotation_merged"]))[table(merged@meta.data[which(merged@meta.data$gate == "38-"),"annotation_merged"]) > 20]


for (ext in extensions) {
  #  ext = 3
  transfer.anchors <- FindTransferAnchors(
    reference = subset(merged, subset = gate == "38-" & annotation_merged %in% RNAcatsign),
    query = test,
    reduction = 'cca',
    query.assay = paste0("RNA_geneprom_",ext,"kb"),
    features = rownames(test@assays[[paste0("RNA_geneprom_",ext,"kb")]]@counts),
    verbose = T
  )


  kweight = 6
    cutoff = 0.4 
    predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = subset(merged, subset = gate == "38-" & annotation_merged %in% RNAcatsign)$annotation_merged,  
    weight.reduction = test[['harmony']],
    dims = 1:50,
    k.weight = kweight
    
  )
  
  
  datatable.test = data.table(predicted.labels, keep.rownames = T)


  cutoff = 0.4

  datatable.test$predicted.id = ifelse(datatable.test$prediction.score.max < cutoff, "notsign", datatable.test$predicted.id)

  datatable.test = as.data.frame(datatable.test)
  rownames(datatable.test) = datatable.test$rn
  
  test <- AddMetaData(object = test, metadata = as.data.frame(datatable.test))
}
    
saveRDS(test, file = "/datadisk/Desktop/Projects/SCRNATAC/scATAC/files/test.merged234.Snap.HA.ChromVar.Integr3kb.Rds")

```   






