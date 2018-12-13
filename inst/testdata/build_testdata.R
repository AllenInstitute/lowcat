library(lowcat)
library(GenomicRanges)
options(stringsAsFactors = FALSE)

inst_bams <- c(system.file("inst/testdata/bam/", "N702-N504_S13_L001_R1_001.rmd.srt.bam", package = "lowcat"),
               system.file("inst/testdata/bam/", "N702-N506_S15_L001_R1_001.rmd.srt.bam", package = "lowcat"),
               system.file("inst/testdata/bam/", "N702-N507_S16_L001_R1_001.rmd.srt.bam", package = "lowcat"))

bam_to_fragment_list(inst_bams)

bed1 <- data.frame(chr = c("chr1","chr1","chr2","chr6"),
                   start = c(10000L,20000L,5000L,2000L),
                   end = c(15000L,25000L,12000L,6000L),
                   name = c("site_1a","site_1b","site_1c","site_1d"),
                   score = c(0,20,50,3),
                   strand = c("+","-","+","-"))

bed2 <- data.frame(chr = c("chr1","chr1","chr4","chr6"),
                   start = c(12000L,2000L,5000L,1000L),
                   end = c(16000L,2500L,12000L,2000L),
                   name = c("site_2a","site_2b","site_2c","site_2d"),
                   score = c(0,20,50,3),
                   strand = c("-","+","+","-"))

gr1 <- GRanges(seqnames = bed1$chr,
               IRanges(start = bed1$start,
                       end = bed1$end),
               strand = bed1$strand)

gr2 <- GRanges(seqnames = bed2$chr,
               IRanges(start = bed2$start,
                       end = bed2$end),
               strand = bed2$strand)

ucsc1 <- c("chr1:10,000-15,000",
           "chr1:20,000-25,000",
           "chr2:5,000-12,000",
           "chr6:2,000-6,000")

ucsc2 <- c("chr1:12,000-16,000",
           "chr1:2,000-2,500",
           "chr4:5,000-12,000",
           "chr6:1,000-2,000")

save(bed1,bed2,gr1,gr2,ucsc1,ucsc2, file = "inst/testdata/rda/simple_regions.rda")
