% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trimFilter.R
\name{TrimAndFilter}
\alias{TrimAndFilter}
\title{Trim low-quality bases out of each reads then remove low-quality reads}
\usage{
TrimAndFilter(
  fl1,
  fl2 = NULL,
  shortreadRAM = 1e+08,
  endTrimThrs = "?",
  endTrimMeanThrs = "4",
  widthThrs = 14L,
  meanQualThrs = 20L,
  totalNThrs = 4L,
  cmplxThrs = 0.5,
  innerNThrs = 2L,
  innerS = 4L,
  innerE = -4L,
  adpter1Seq = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
  adpter2Seq = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
  destination1 = paste0(gsub(fl1, pattern = "\\\\.fastq\\\\.gz$", replacement = ""),
    "_trimed.fastq.gz"),
  destination2 = paste0(gsub(fl2, pattern = "\\\\.fastq\\\\.gz$", replacement = ""),
    "_trimed.fastq.gz"),
  ...
)
}
\arguments{
\item{fl1}{full path of an fastq.gz file to be processed}

\item{fl2}{NULL for single-end, full path of the paired fastq.gz for paired-end}

\item{shortreadRAM}{RAM limit in bytes for ShortRead, default 1e8}

\item{endTrimThrs}{Phred score threshold of the end base below which the end base will be trimmed,default "?"}

\item{widthThrs}{width threshold of reads to be filtered out,default 14L}

\item{cmplxThrs}{complexity threshold of reads to be filtered out,default 0.5 i.e. half of mean complexity of human genome}

\item{innerNThrs}{number of N inside the read below which the read will be removed,default 2L}

\item{innerS}{which base should be considered the start of inside part of the read after the above trimming,default 4L}

\item{innerE}{which base should be considered the end of inside part of the read after the above trimming.negative integer X for width(read)-abs(X) default -4L;postive integer X for Xth base}

\item{adpter1Seq, adpter2Seq}{adapter sequence to be trimmed from end and also inner part of the read}

\item{destination1, destination2}{the output fastq.gz full path, default paste0(gsub(inputedFullPath,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz")}

\item{endTrimThrsend}{mean Phred score threshold of five bases in ends,below which the five bases will be trimmed,default "4"}
}
\value{
Date.frame summarizing trimming and filtration.
}
\description{
This function based on low-level package ShortRead provides high flexibility
to trim and filter according to QC report . For even higher extent of customization,
see vignette using sub-level functions wrapped in this function
}
\examples{
repDataFrame<-TrimAndFilter(fl1="data/a_1.fastq.gz",fl1="data/a_2.fastq.gz")
}
