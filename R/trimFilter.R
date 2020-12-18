#' Trim low-quality bases out of each reads then remove low-quality reads
#'
#' This function based on low-level package ShortRead provides high flexibility
#' to trim and filter according to QC report . For even higher extent of customization,
#' see vignette using sub-level functions wrapped in this function
#'
#' @param fl1 full path of an fastq.gz file to be processed
#' @param fl2 NULL for single-end, full path of the paired fastq.gz for paired-end
#' @param shortreadRAM RAM limit in bytes for ShortRead, default 1e8
#' @param endTrimThrs Phred score threshold of the end base below which the end base will be trimmed,default "?"
#' @param endTrimThrsend mean Phred score threshold of five bases in ends,below which the five bases will be trimmed,default "4"
#' @param adpter1Seq,adpter2Seq adapter sequence to be trimmed from end and also inner part of the read
#' @param widthThrs width threshold of reads to be filtered out,default 14L
#' @param cmplxThrs complexity threshold of reads to be filtered out,default 0.5 i.e. half of mean complexity of human genome
#' @param innerNThrs number of N inside the read below which the read will be removed,default 2L
#' @param innerS which base should be considered the start of inside part of the read after the above trimming,default 4L
#' @param innerE which base should be considered the end of inside part of the read after the above trimming.negative integer X for width(read)-abs(X) default -4L;postive integer X for Xth base
#' @param destination1,destination2 the output fastq.gz full path, default paste0(gsub(inputedFullPath,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz")
#' @return Date.frame summarizing trimming and filtration.
#' @examples
#' repDataFrame<-TrimAndFilter(fl1="data/a_1.fastq.gz",fl1="data/a_2.fastq.gz")

#' @export
TrimAndFilter <-  function(fl1,fl2=NULL,shortreadRAM=1e8,
                           endTrimThrs="?",endTrimMeanThrs="4",
                           widthThrs=14L,meanQualThrs=20L,
                           totalNThrs=4L,cmplxThrs=0.5,innerNThrs=2L,
                           innerS=4L,innerE=-4L,
                           adpter1Seq="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                           adpter2Seq="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
                           destination1=paste0(gsub(fl1,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),
                           destination2=paste0(gsub(fl2,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),...) {

  nthreads <- .Call(ShortRead:::.set_omp_threads, 1L)
  on.exit(.Call(ShortRead:::.set_omp_threads, nthreads))

  filterReport <- c(totalSequences = 0, endTrimed=0, matchTo5pAdapter = 0,
                    matchTo3pAdapter = 0, tooShort = 0, lowMeanQuality=0,
                    tooManyInternalN = 0, tooManyN = 0, lowComplexity = 0,
                    totalPassed = 0)

  isPaired<-!is.null(fl2)

  if (isPaired) {
    ## open input stream
    stream1 <- open(FastqStreamer(fl1,readerBlockSize = shortreadRAM,verbose = T))
    on.exit(close(stream1))
    stream2 <- open(FastqStreamer(fl2,readerBlockSize = shortreadRAM,verbose = T))
    on.exit(close(stream2))

    filterReport<-cbind(file1=filterReport,file2=filterReport,file1_2=filterReport)

    repeat {

      ## input chunk
      fq1 <- yield(stream1)
      if (length(fq1) == 0)
        break
      fq2 <- yield(stream2)

      filterReport["totalSequences",]<-filterReport["totalSequences",]+1e6L

      #trimEnd
      rng1a<-trimEnds(fq1, endTrimThrs,ranges = T)#30
      rng2a<-trimEnds(fq2, endTrimThrs,ranges = T)

      ## trim as soon as 2 of 5 nucleotides has quality encoding less
      ## than "4" (phred score 20)
      rng1b <- trimTailw(fq1, 2, endTrimMeanThrs, 2,ranges = T)
      rng2b <- trimTailw(fq2, 2, endTrimMeanThrs, 2,ranges = T)

      s1<-pmax(start(rng1a),start(rng1b))
      s2<-pmax(start(rng2a),start(rng2b))
      e1<-pmin(end(rng1a),end(rng1b))
      e2<-pmin(end(rng2a),end(rng2b))

      filterReport['endTrimed',1:2]<-  filterReport['endTrimed',1:2]+
        c(sum(s1!=1|e1!=width(fq1)),sum(s2!=1|e1!=width(fq2)))

      w1<-pmax(e1-s1+1,0)
      w2<-pmax(e2-s2+1,0)
      fq1<-narrow(fq1,start = s1,width = w1)
      fq2<-narrow(fq2,start = s2,width = w2)

      res1<-trimAdapter(shortQobj = fq1,Lpattern = adpter1Seq ,Rpattern =adpter1Seq)
      fq1<-res1$shortQobj
      res2<-trimAdapter(shortQobj = fq2,Lpattern = adpter2Seq ,Rpattern =adpter2Seq)
      fq2<-res2$shortQobj

      filterReport['matchTo5pAdapter',1:2]<-  filterReport['matchTo5pAdapter',1:2]+
        c(res1$matchL,res2$matchL)

      filterReport['matchTo3pAdapter',1:2]<-  filterReport['matchTo3pAdapter',1:2]+
        c(res1$matchR,res2$matchR)

      #filter too short
      l1<-filterWidth(threshold = widthThrs)(fq1)
      l2<-filterWidth(threshold = widthThrs)(fq2)
      filterReport['tooShort',]<-filterReport['tooShort',]+
        as.vector(c(stats(l1)['Input']-stats(l1)['Passing'],stats(l2)['Input']-stats(l2)['Passing'],
                    stats(l1)['Input']-sum(l1&l2)),mode="integer")
      fq1<-fq1[l1&l2]
      fq2<-fq2[l1&l2]

      #filter lowMeanQuality
      l1<-filterLowQuality(threshold = meanQualThrs)(fq1)
      l2<-filterLowQuality(threshold = meanQualThrs)(fq2)
      filterReport['lowMeanQuality',]<-filterReport['lowMeanQuality',]+
        as.vector(c(stats(l1)['Input']-stats(l1)['Passing'],stats(l2)['Input']-stats(l2)['Passing'],
                    stats(l1)['Input']-sum(l1&l2)),mode="integer")
      fq1<-fq1[l1&l2]
      fq2<-fq2[l1&l2]


      #filter tooManyN
      l1<-nFilter(threshold = totalNThrs)(fq1)
      l2<-nFilter(threshold = totalNThrs)(fq2)
      filterReport['tooManyN',]<-filterReport['tooManyN',]+
        as.vector(c(stats(l1)['Input']-stats(l1)['Passing'],stats(l2)['Input']-stats(l2)['Passing'],
                    stats(l1)['Input']-sum(l1&l2)),mode="integer")
      fq1<-fq1[l1&l2]
      fq2<-fq2[l1&l2]


      #filter lowComplexity
      l1<-filterLowComplexity(threshold =cmplxThrs)(fq1)
      l2<-filterLowComplexity(threshold =cmplxThrs)(fq2)
      filterReport['lowComplexity',]<-filterReport['lowComplexity',]+
        as.vector(c(stats(l1)['Input']-stats(l1)['Passing'],stats(l2)['Input']-stats(l2)['Passing'],
                    stats(l1)['Input']-sum(l1&l2)),mode="integer")
      fq1<-fq1[l1&l2]
      fq2<-fq2[l1&l2]

      #filter tooManyInternalN
      if (innerE<=0){
        l1<-filterInnerN(s = innerS,e = width(fq1)+innerE,thrs = innerNThrs)(fq1)
        l2<-filterInnerN(s = innerS,e = width(fq2)+innerE,thrs = innerNThrs)(fq2)
      }else{
        l1<-filterInnerN(s = innerS,e = innerE,thrs = innerNThrs)(fq1)
        l2<-filterInnerN(s = innerS,e = innerE,thrs = innerNThrs)(fq2)
      }


      filterReport['tooManyInternalN',]<-filterReport['tooManyInternalN',]+
        as.vector(c(stats(l1)['Input']-stats(l1)['Passing'],stats(l2)['Input']-stats(l2)['Passing'],
                    stats(l1)['Input']-sum(l1&l2)),mode="integer")
      fq1<-fq1[l1&l2]
      fq2<-fq2[l1&l2]

      filterReport['totalPassed',]<-filterReport['totalPassed',]+
        as.vector(c(stats(l1)['Passing'],stats(l2)['Passing'],sum(l1&l2)),mode="integer")


      #chunk output
      writeFastq(fq1, destination1, "a")
      writeFastq(fq2, destination2, "a")
    }
  }else{
    ## open input stream
    stream1 <- open(FastqStreamer(fl1,readerBlockSize = 1e8,verbose = T))
    on.exit(close(stream1))

    repeat {

      ## input chunk
      fq1 <- yield(stream1)
      if (length(fq1) == 0)
        break


      filterReport["totalSequences"]<-filterReport["totalSequences"]+1e6L

      #trimEnd
      rng1a<-trimEnds(fq1, endTrimThrs,ranges = T)#30

      ## trim as soon as 2 of 5 nucleotides has quality encoding less
      ## than "4" (phred score 20)
      rng1b <- trimTailw(fq1, 2, endTrimMeanThrs, 2,ranges = T)

      s1<-pmax(start(rng1a),start(rng1b))
      e1<-pmin(end(rng1a),end(rng1b))

      filterReport['endTrimed']<-  filterReport['endTrimed']+
        sum(s1!=1|e1!=width(fq1))

      w1<-pmax(e1-s1+1,0)

      fq1<-narrow(fq1,start = s1,width = w1)


      res1<-trimAdapter(shortQobj = fq1,Lpattern = adpter1Seq ,Rpattern =adpter1Seq,...)
      fq1<-res1$shortQobj

      filterReport['matchTo5pAdapter']<-  filterReport['matchTo5pAdapter']+
        res1$matchL

      filterReport['matchTo3pAdapter']<-  filterReport['matchTo3pAdapter']+
        res1$matchR

      #filter too short
      l1<-filterWidth(threshold = widthThrs)(fq1)
      filterReport['tooShort']<-filterReport['tooShort']+
        as.integer((stats(l1)['Input']-stats(l1)['Passing']))
      fq1<-fq1[l1]


      #filter lowMeanQuality
      l1<-filterLowQuality(threshold = meanQualThrs)(fq1)
      filterReport['lowMeanQuality']<-filterReport['lowMeanQuality']+
        as.integer((stats(l1)['Input']-stats(l1)['Passing']))
      fq1<-fq1[l1]



      #filter tooManyN
      l1<-nFilter(threshold = totalNThrs)(fq1)
      filterReport['tooManyN']<-filterReport['tooManyN']+
        as.integer((stats(l1)['Input']-stats(l1)['Passing']))
      fq1<-fq1[l1]



      #filter lowComplexity
      l1<-filterLowComplexity(threshold =cmplxThrs,...)(fq1)
      filterReport['lowComplexity']<-filterReport['lowComplexity']+
        as.integer((stats(l1)['Input']-stats(l1)['Passing']))
      fq1<-fq1[l1]


      #filter tooManyInternalN
      if (innerE<=0){
      l1<-filterInnerN(s = innerS,e = width(fq1)+innerE)(fq1)
      }else{
        l1<-filterInnerN(s = innerS,e = width(fq1)+innerE)(fq1)
        }
      filterReport['tooManyInternalN']<-filterReport['tooManyInternalN']+
        as.integer((stats(l1)['Input']-stats(l1)['Passing']))
      fq1<-fq1[l1]


      filterReport['totalPassed']<-filterReport['totalPassed']+
        as.integer((stats(l1)['Passing']))


      #chunk output
      writeFastq(fq1, destination1, "a")
    }
  }

  return(filterReport)
}

#' Trim adapter sequences out of reads
#' @export
trimAdapter<-function(shortQobj,Lpattern="",Rpattern="",
                      max.Lmismatch=rep(0:2, c(6,3,100)),
                      max.Rmismatch=rep(0:2, c(6,3,100)),
                      with.Lindels=F, with.Rindels=F){
  numNs <- 90
  if (nchar(Lpattern) > 0) {
    max.Lmismatch <- max.Lmismatch[1:nchar(Lpattern)]
    max.Lmismatch <- c(max.Lmismatch, 1:numNs + max(max.Lmismatch))
    Lpattern <- paste(c(rep("N", numNs), Lpattern), collapse = "")
  }
  if (nchar(Rpattern) > 0) {
    max.Rmismatch <- max.Rmismatch[1:nchar(Rpattern)]
    max.Rmismatch <- c(max.Rmismatch, 1:numNs + max(max.Rmismatch))
    Rpattern <- paste(c(Rpattern, rep("N", numNs)), collapse = "")
  }

  ranges <- trimLRPatterns(subject=shortQobj, Lpattern=Lpattern, Rpattern=Rpattern,
                           max.Lmismatch=max.Lmismatch, max.Rmismatch=max.Rmismatch,
                           with.Lindels=with.Lindels, with.Rindels=with.Rindels,
                           ranges=TRUE)

  list(shortQobj= narrow(x=shortQobj, start=start(ranges), end=end(ranges)),
       matchL= sum(start(ranges) != 1),
       matchR= sum(end(ranges) != width(shortQobj)))
}


#' Trim reads containing too few bases
#' @export
filterWidth <-
  function(threshold=14L, .name = "WidthFilter"){
    srFilter(function(x) {
      width(x) >= threshold
    }, name = .name)
  }

#' Trim reads containing too many N at inner part of the reads
#' @export
filterInnerN <-
  function(s=3L,e=12L,thrs=1L, .name = "InnerNFilter"){
    srFilter(function(x) {
      nFilter(threshold = thrs)(narrow(x,start = s,end = e))
    }, name = .name)
  }

#' Trim reads with low mean quality of bases
#' @export
filterLowQuality <-
  function(threshold=20, .name = "LowQualityFilter"){
    srFilter(function(x) {
      rowMeans(as(quality(x), "matrix"),na.rm = T) > threshold
    }, name = .name)
  }

#' Trim reads with low complexity
#' @export
filterLowComplexity <-
  function(threshold=0.5, referenceEntropy=3.908135, .name = "LowComplexityFilter"){
    srFilter(function(x) {
      # less than half the average entropy per dinucleotide of non-random chromosomes in hg18 (3.908135 bits):
      #   entropy(x)/3.908135 >= threshold
      diNucFreq <- dinucleotideFrequency(sread(x))
      if(is.null(dim(diNucFreq))) {
        diNucFreq <- diNucFreq/sum(diNucFreq)
        H <- -sum(diNucFreq * log2(diNucFreq), na.rm = TRUE)
      } else {
        diNucFreq <- diNucFreq/rowSums(diNucFreq)
        H <- -rowSums(diNucFreq * log2(diNucFreq), na.rm = TRUE)
      }
      H/referenceEntropy >= threshold
    }, name = .name)
  }


















