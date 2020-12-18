#' Master function wrapping trimming bases, filtering reads, alignment and feature counting sequentially
#'
#'
#' @param fileList1 full path of one or multiple fastq.gz files to be processed
#' @param fileList2 default NULL for single-end, full path of the paired fastq.gz files for paired-end
#' @param shortreadRAM RAM limit in bytes for ShortRead, default 1e8
#' @param subReadThreads number of threads to be used by SubRead, default 3L
#' @param alignType default "rna" for RNAseq and microRNAseq,see other option in vignettes, it will be passed to the "type" argument of SubRead::align
#' @param genomeRefFile,genomeAnnotFile full path of refseq sequence fastq file to generate alignment index, and its mate .gtf annotation file for feature counting.These files can be downloaded from NCBI,see vignettes
#' @param indexBasename if you never build index before, leave as defaul "NULL", otherwise set it as "my_index" and make sure index files are in working dir to speed up by avoiding regenerating index files
#' @param alignPairedOutput only works for paired-end input,default is gsub(basename(fileList1),pattern ="_1.*\\.fastq\\.gz",replacement = "\\.bam"), optimize it please if it happens to cause overwriting of output files
#'
#' @return A list contains Date.frame summarizing trimming and filtration of each input fastq.gz file. Generated .bam file will be save to ./bam; feature counts are saved as .rds file in working dir
#' @examples
#' trimFilterRes<-seqW(fileList1= fl,subReadThreads = 3L,shortreadRAM=1e8)

#' @export
seqW<- function(fileList1, fileList2=NULL,alignType="rna",subReadThreads=3L,shortreadRAM=1e8,
                genomeRefFile="./GCF_000001405.26_GRCh38_genomic.fna.gz",
                genomeAnnotFile="./GCF_000001405.39_GRCh38.p13_genomic.gtf.gz",
                indexBasename=NULL,
                alignPairedOutput=gsub(basename(fileList1),pattern ="_1.*\\.fastq\\.gz",replacement = "\\.bam"),
                ...){
  repList<-vector(mode="list",length=length(fileList1))
  for (i in 1:length(fileList1)){
    res<-TrimAndFilter(fl1 = fileList1[i],fl2 = fileList2[i],...)
    repList[i]<-list(res)}

  if (!alignType %in% c('dna','rna','dnaLong','rnaExon_Exon','microRNA')){
  stop("alignType must be one of 'dna' 'rna' 'dnaLong' 'rnaExon_Exon' 'microRNA'")
  }else{
    if (alignType %in% c("dna","rna")){
      if (is.null(indexBasename)){buildindex(basename="my_index",reference=genomeRefFile)
        indexBasename<-"my_index"}
      dir.create("./bam")
      if (is.null(fileList2)) {
        align(index=indexBasename,readfile1=paste0(gsub(fileList1,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),
              readfile2=NULL,
              output_file = file.path("./bam", gsub(basename(fileList1),pattern ="\\.fastq\\.gz$",replacement = "\\.bam")),
              type=alignType,nthreads=subReadThreads,...)

        saveRDS(object = featureCounts(files=file.path("./bam", gsub(basename(fileList1),pattern ="\\.fastq\\.gz$",replacement = "\\.bam")),
                                       isPairedEnd=F,nthreads=subReadThreads,
                                       annot.ext =genomeAnnotFile,
                                       isGTFAnnotationFile = T,
                                       verbose = T,...),
                file = paste0(alignType,"FeatureCount.rds"))

      }else{
        align(index=indexBasename,readfile1=paste0(gsub(fileList1,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),
              readfile2=paste0(gsub(fileList2,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),
              output_file = file.path("./bam",alignPairedOutput),
              type=alignType,nthreads=subReadThreads,...)

        saveRDS(object = featureCounts(files=file.path("./bam",alignPairedOutput),
                                       isPairedEnd=T,nthreads=subReadThreads,
                                       annot.ext =genomeAnnotFile,
                                       isGTFAnnotationFile = T,
                                       verbose = T,...),
                file = paste0(alignType,"FeatureCount.rds"))
      }
    }else{
      if (alignType=="dnaLong"){
        if (is.null(indexBasename)){buildindex(basename="my_index",reference=genomeRefFile)
          indexBasename<-"my_index"}
        dir.create("./bam")
        if (is.null(fileList2)) {
          sublong(index=indexBasename,readfile1=paste0(gsub(fileList1,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),
                readfile2=NULL,
                output_file = file.path("./bam", gsub(basename(fileList1),pattern ="\\.fastq\\.gz$",replacement = "\\.bam")),
                type="dna",nthreads=subReadThreads,...)

          saveRDS(object = featureCounts(files=file.path("./bam", gsub(basename(fileList1),pattern ="\\.fastq\\.gz$",replacement = "\\.bam")),
                                         isPairedEnd=F,nthreads=subReadThreads,
                                         annot.ext =genomeAnnotFile,
                                         isGTFAnnotationFile = T,
                                         verbose = T,...),
                  file = paste0(alignType,"FeatureCount.rds"))

        }else{
          sublong(index=indexBasename,readfile1=paste0(gsub(fileList1,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),
                readfile2=paste0(gsub(fileList2,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),
                output_file = file.path("./bam",alignPairedOutput),
                type="dna",nthreads=subReadThreads,...)

          saveRDS(object = featureCounts(files=file.path("./bam",alignPairedOutput),
                                         isPairedEnd=T,nthreads=subReadThreads,
                                         annot.ext =genomeAnnotFile,
                                         isGTFAnnotationFile = T,
                                         verbose = T,...),
                  file = paste0(alignType,"FeatureCount.rds"))
        }
      }
      if (alignType=="rnaExon_Exon"){
        if (is.null(indexBasename)){buildindex(basename="my_index",reference=genomeRefFile)
          indexBasename<-"my_index"}
        dir.create("./bam")
        if (is.null(fileList2)) {
          subjunc(index=indexBasename,readfile1=paste0(gsub(fileList1,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),
                  readfile2=NULL,
                  output_file = file.path("./bam", gsub(basename(fileList1),pattern ="\\.fastq\\.gz$",replacement = "\\.bam")),
                  type="rna",nthreads=subReadThreads,...)

          saveRDS(object = featureCounts(files=file.path("./bam", gsub(basename(fileList1),pattern ="\\.fastq\\.gz$",replacement = "\\.bam")),
                                         isPairedEnd=F,nthreads=subReadThreads,
                                         annot.ext =genomeAnnotFile,
                                         isGTFAnnotationFile = T,
                                         verbose = T,...),
                  file = paste0(alignType,"FeatureCount.rds"))

        }else{
          subjunc(index=indexBasename,readfile1=paste0(gsub(fileList1,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),
                  readfile2=paste0(gsub(fileList2,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),
                  output_file = file.path("./bam",alignPairedOutput),
                  type="rna",nthreads=subReadThreads,...)

          saveRDS(object = featureCounts(files=file.path("./bam",alignPairedOutput),
                                         isPairedEnd=T,nthreads=subReadThreads,
                                         annot.ext =genomeAnnotFile,
                                         isGTFAnnotationFile = T,
                                         verbose = T,...),
                  file = paste0(alignType,"FeatureCount.rds"))
        }
      }
      if (alignType=="microRNA"){
        if (is.null(indexBasename)){buildindex(basename="my_index",reference=genomeRefFile,gappedIndex =F,indexSplit = F)
          indexBasename<-"my_index"}
        dir.create("./bam")
        if (is.null(fileList2)) {
          subjunc(index=indexBasename,readfile1=paste0(gsub(fileList1,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),
                  readfile2=NULL,
                  output_file = file.path("./bam", gsub(basename(fileList1),pattern ="\\.fastq\\.gz$",replacement = "\\.bam")),
                  type="rna",nthreads=subReadThreads,...)

          saveRDS(object = featureCounts(files=file.path("./bam", gsub(basename(fileList1),pattern ="\\.fastq\\.gz$",replacement = "\\.bam")),
                                         isPairedEnd=F,nthreads=subReadThreads,
                                         annot.ext =genomeAnnotFile,
                                         isGTFAnnotationFile = T,
                                         verbose = T,...),
                  file = paste0(alignType,"FeatureCount.rds"))

        }else{
          subjunc(index=indexBasename,readfile1=paste0(gsub(fileList1,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),
                  readfile2=paste0(gsub(fileList2,pattern ="\\.fastq\\.gz$",replacement = ""),"_trimed.fastq.gz"),
                  output_file = file.path("./bam",alignPairedOutput),
                  type="rna",nthreads=subReadThreads,...)

          saveRDS(object = featureCounts(files=file.path("./bam",alignPairedOutput),
                                         isPairedEnd=T,nthreads=subReadThreads,
                                         annot.ext =genomeAnnotFile,
                                         isGTFAnnotationFile = T,
                                         verbose = T,...),
                  file = paste0(alignType,"FeatureCount.rds"))
        }
      }
    }
  }

  return(repList)
}
