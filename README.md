
SeqWins Vignette
=============================
An R package allowing flexible base trimming and complete Fastq analysis on Windows System

> 💡 **Tip: [Please find **👉 MY BLOG** for an introduction and complete view of the project behind the code in this repository.](https://myhugoblog)**

### Contributors

<table>
  <tr>
    <td align="center">
      <a href="https://github.com/jliu678">
        <img src="https://avatars.githubusercontent.com/u/53794392?v=4" width="50px;" alt="jliu678"/><br />
        <sub><b>jliu678</b></sub>
      </a>
      <br />🙌😎🔬👨‍💻
    </td>
  </tr>
</table>

### Introduction

Words prevail that Fastq data cannot be elegantly processed in Windows. However, the fundamental low-level R package ShortRead and Rsubread have been available to provide memory-efficient, chunk-wise processing of FASTQ files and offer alignment performance that is competitive with or faster than many Linux-based aligners. Building on these strengths, I developed SeqWins (fastq **Seq**uence analysis on **Win**dows system), an R package for flexible base trimming and comprehensive FASTQ analysis on Windows. It achieves on **pure Windows system**:
 - memory-efficient, chunk-wise processing of FASTQ files
 - alignment performance competitive with or faster than many Linux-based aligners
 - flexible base-level (ATCG) quality control and trimming
 - convenient high-level whole-process analysis of fastq data-- spanning quality control report, trimming bases and reads accordingly, alignment and feature count



### Keywords

  + [SeqWins](https://myhugoblog)
  + Fastq
  + Windows system
  + R
  + Trim base quality
  + Filter reads
  + RNAseq
  + Quality control report
  + DNAseq
  + CHIPseq
  
***


### Installation

```{r init}
# Install dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ShortRead")

BiocManager::install("Rsubread")

# Install and load SeqWins
devtools::install_github("jliu678/SeqWins")
library(SeqWins)
```

***



### Prepare genomic Fastq sequence and the corresonding GTF annotation

These gemonic data of the same species as with your fastq data are required to build index when aligning reads and get feature counts, and can be download from NCBI assembly website, for example [the hg38 files](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/) by clicking the top right "**Download Assembly**" button and selecting "**Genomic FASTA (.fna)**" and "**Genomic GTF (.gtf)**" sequentially on the resultant drop-down options.

If necessary untar the downloaded files as blow

```{r}
untar("full/path/genome_assemblies_genome_fasta.tar",exdir = ".")
untar("full/path/genome_assemblies_genome_gtf.tar",exdir = ".")
```


### QC report

set working directory, it will be where the QC result folder (ShortRead 1.46.0 name it as "ShortRead Quality Assessment_files"),index files and the folder named "bam" storing aligned data are located

```{r}
setwd("/your/working/directory")
```

and get QC report

```{r}
# call ShortRead::qa to generate QC reports of all fastq.gz files
qa <- qa("full/path/FastqFolder", "fastq.gz")
browseURL(report(qa))
```
the above will generate report of the fastq files retrieved by
```{r}
list.files(path = "full/path/FastqFolder",pattern = "fastq.gz")
```

***

### Trim bases,filter reads,align and count feature

According to the QC report, customize the parameters used for trimming bases, filtering reads. Although most can be left default,please specify file paths totally decided by yourself; the  "subReadThreads","shortreadRAM" decided by your computer; and the sequence tech type decided by your project,hopefully not by money.

#### *specify file path*
If your fastq were generated from single-end sequencing (old-fashioned you are!), only specify `fileList1` as blow
```{r}
fl<-list.files(path = "full/path/FastqFolder",pattern = "fastq.gz")
seqW(fileList1=fl,genomeRefFile="./GCF_000001405.26_GRCh38_genomic.fna.gz",
                genomeAnnotFile="./GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")#RNAseq
```

If paired-end,specify both `fileList1` i.e. the "read 1" files and `fileList2` i.e. "read 2" files
```{r}
fl_1<-list.files(path = "full/path/FastqFolder",pattern = "^.*_1\\.fastq\\.gz$")
fl_2<-list.files(path = "full/path/FastqFolder",pattern = "^.*_2\\.fastq\\.gz$")

seqW(fileList1=fl_1,fileList1=fl_2,genomeRefFile="./GCF_000001405.26_GRCh38_genomic.fna.gz",
                genomeAnnotFile="./GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")#RNAseq
```

Please don't forget to specify `genomeRefFile` and `genomeAnnotFile` i.e. the path of "**Genomic FASTA (.fna)**" and "**Genomic GTF (.gtf)**" as any of the above examples for the first time when you run `seqW`.

If you already have index files located in the working dir, you can speed it up by setting `indexBasename="my_index"` to avoid regenerating index files like below. And this will make `seqW` function ignore whatever is assgned to `genomeRefFile`. 

```{r}
seqW(fileList1=fl_1,indexBasename="my_index",
                genomeAnnotFile="./GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")#RNAseq
```

`alignPairedOutput` only works for paired-end input,default is `gsub(basename(fileList1),pattern ="_1.*\\.fastq\\.gz",replacement = "\\.bam")`, optimize it please if it happens to cause overwriting of output files. For example for fastq files named as "a_1_sample1.fastq.gz",  "a_2_sample1.fastq.gz", "a_1_sample2.fastq.gz", as "a_2_sample2.fastq.gz", sample1 and sample2 files will both result in "a.bam" and either overwriting or error will occur.


#### *specify trim and filter*
Below shows a complete list of trim and filter parameters. Probably most can be left as default except the adaptor sequence specific to your sequence platform

1. `endTrimThrs`

Phred score threshold of the end base below which the end base will be trimmed,default "?"

2. `endTrimThrsend`

mean Phred score threshold of five bases in ends,below which the five bases will be trimmed,default "4"

3. `adpter1Seq`,`adpter2Seq` 

adapter sequence to be trimmed from end and also inner part of the read

+ default `adpter1Seq` "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" i.e. illumina HT4000 adapter for read 1

+ default `adpter2Seq` "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" i.e. illumina HT4000 adapter for read 2

4. `with.Lindels`,`with.Rindels`

True if your reads contain indels on left end (5',corresponding to `with.Lindels`) or right end (3',corresponding to `with.Rindels`), default False

5. `widthThrs`

width threshold of reads to be filtered out,default 14L

6. `cmplxThrs`

complexity threshold of reads to be filtered out,default 0.5 i.e. half of mean complexity of human genome

7. `innerNThrs`

number of N inside the read below which the read will be removed,default 2L

+ `innerS`

  which base should be considered the start of inside part of the read after the above trimming,default 4L

+ `innerE`

  which base should be considered the end of inside part of the read after the above trimming. Negative integer X for width(read)-abs(X) default -4L; positive integer X for Xth base

#### *summary of trimming and filtration*

A list containing dataframes of which each reports the result of trimming and filtration will be returned by the `seqw` function. You can assign the list to an object for further observation.
```{r}
reportList<-seqW(fileList1=fl_1,fileList1=fl_2,genomeRefFile="./GCF_000001405.26_GRCh38_genomic.fna.gz",
                genomeAnnotFile="./GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")#RNAseq
```

#### *set alignType*

Specifically for RNAseq leave `alignType` as default (simply don't mention it). By the way, unit of "shortreadRAM" is byte,don't panic.

```{r}
seqW(fileList1=fl,subReadThreads=3L,shortreadRAM=1e8,
                genomeRefFile="./GCF_000001405.26_GRCh38_genomic.fna.gz",
                genomeAnnotFile="./GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")
```

For DNAseq, longDNAseq, RNAseq aim at exon-exon junction and microRNAseq, set `alignType` as '`dna`','`dnaLong`','`rnaExon_Exon`'and '`microRNA`' respectively--
```{r}
seqW(fileList1=fl,alignType ='dna')#DNAseq
seqW(fileList1=fl,alignType ='dnaLong')#longDNAseq
seqW(fileList1=fl,alignType ='rnaExon_Exon')#RNAseq aim at exon-exon junction
seqW(fileList1=fl,alignType ='microRNA')#microRNAseq
```

#### *set featureCount*

If only gene counts are needed to be output, leave `useMetaFeatures` as default (simply don't mention it). Please note each row in the input "**Genomic GTF (.gtf)**" is a feature, and many features can belong to one gene while some features have no corresponding gene ID. You can see all fields that define a feature [here](https://mblab.wustl.edu/GTF22.html). Thus, if feature counts containing information of all the fields rather than gene counts are needed, specify `useMetaFeatures=F` as below, so that the output will give information like the number of reads mapped to "+" strand of chromosome "II" from "380" to "401".


```{r}
seqW(fileList1=fl,subReadThreads=3L,shortreadRAM=1e8,useMetaFeatures=F,
                genomeRefFile="./GCF_000001405.26_GRCh38_genomic.fna.gz",
                genomeAnnotFile="./GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")
```

***
### Even more flexibility is faciliated!

Hope the package can realize easy but deep appreciation of your fastq data, and you probably have found every parameters you wanna touch through above introduction. But there is even more! `...` arguments enable `seqW` function to access all the arguments in the classic functions wrapped inside it, including 

1. `Rsubread::align`, called when `alignType="dna"` ,`"rna"` ,`"microRNA"`
2. `Rsubread::sublong`, called when `alignType="longDNA"`
3. `Rsubread::subjunc`, called when `alignType="rnaExon_Exon"`
4. `Rsubread::featureCounts`, called every time `seqW` runs
(
thus simply add their parameters (see what they have by `?Rsubread::featureCounts` or [Rsubread vignettes](https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf)) when calling `seqW`, e.g.

```{r}
# minFragLength by default is 50L in Rsubread, but you can modify it easily as blow
seqW(minFragLength=40L,fileList1=fl,genomeRefFile="./GCF_000001405.26_GRCh38_genomic.fna.gz",
                genomeAnnotFile="./GCF_000001405.39_GRCh38.p13_genomic.gtf.gz")#RNAseq
```
#### *all are exported to facilitate your own pipeline*
You might even wanna add/alter the trimming and filtration pipeline; No problem! I have exported every units wrapped under the top-level function, including 

1. `trimEnds`
2. `trimTailw`
3. `trimAdapter`
4. `filterWidth`
5. `filterLowQuality`
6. `ShortRead::nFilter`
7. `filterLowComplexity`
8. `filterInnerN`

simply `?trimEnds` for example, to see their function and how to use it as bricks to make your own trimming and filter flow. And then use aforementioned `Rsubread` functions to complete alignment and feature count.


#### *run from trimmed fastq.gz or aligned bam*

`seqW` will sequentially generate

1. index files in working dir, whose names by default start with "my_index"; 
2. trimmed ".fastq.gz" files whose names contain "trimed" in the same folder as for you initially input ".fastq.gz" files;
3. aligned ".bam" files in the folder named "bam" located in working dir;
4. feature counts in ".rds" format in working dir.

The above steps are carried out by `Rsubread::buildindex`, `TrimAndFilter`, `Rsubread::align` or its variant and `Rsubread::featureCounts`, respectively. We have earlier discussed how to use `seqW` funciton avoid repeating generation of index files when there are already index files. Besides, if satisfactory trimmed fastq exist, simply sequentially call necessary functions listed above to achieve your goal, e.g. building index by `Rsubread::buildindex` then aligning and count features by  `Rsubread::align` and `Rsubread::featureCounts`. If you even have aligned bam, just run `Rsubread::buildindex` and `Rsubread::featureCounts`. If you have got feature counts at hand, well, thank you a lot for using or visiting this package.  

Also see [html vignetts](https://jliu678.github.io/SeqWins/)


