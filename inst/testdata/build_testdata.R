library(lowcat)

inst_bams <- c(system.file("inst/testdata/bam/", "N702-N504_S13_L001_R1_001.rmd.srt.bam", package = "lowcat"),
               system.file("inst/testdata/bam/", "N702-N506_S15_L001_R1_001.rmd.srt.bam", package = "lowcat"),
               system.file("inst/testdata/bam/", "N702-N507_S16_L001_R1_001.rmd.srt.bam", package = "lowcat"))

frag_list <- bam_to_fragment_list(inst_bams)

