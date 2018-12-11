library(testthat)
library(lowcat)
options(stringsAsFactors = F)

context("Reading BAM files.")

inst_bams <- c(system.file("inst/testdata/bam/", "N702-N504_S13_L001_R1_001.rmd.srt.bam", package = "lowcat"),
               system.file("inst/testdata/bam/", "N702-N506_S15_L001_R1_001.rmd.srt.bam", package = "lowcat"),
               system.file("inst/testdata/bam/", "N702-N507_S16_L001_R1_001.rmd.srt.bam", package = "lowcat"))

test_that("bam_fragments_to_list() reads BAM files as a list of GenomicRanges",
          {
            frag_list <- bam_to_fragment_list(inst_bams)

            expect_is(frag_list, "list")
            expect_is(frag_list[[1]], "GenomicRanges")
            expect_equal(length(frag_list), length(inst_bams))
            expect_equal(names(frag_list), sub(".+/","",inst_bams))
          })
