# lowcat: Low Coverage Accessibility and Transcriptomics

<img src="https://upload.wikimedia.org/wikipedia/commons/f/f7/Rusty-spotted_cat_%28f._rubiginosa%29.JPG" alt="Rusty-spotted lowcat" width="300px"/>

Master: [![Master status](https://travis-ci.org/AllenInstitute/lowcat.svg?branch=master)](https://travis-ci.org/AllenInstitute/lowcat)
[![Master coverage](https://codecov.io/gh/AllenInstitute/lowcat/branch/master/graph/badge.svg)](https://codecov.io/github/AllenInstitute/lowcat?branch=master)


Dev: [![Dev status](https://travis-ci.org/AllenInstitute/lowcat.svg?branch=dev)](https://travis-ci.org/AllenInstitute/lowcat)
[![Dev coverage](https://codecov.io/gh/AllenInstitute/lowcat/branch/dev/graph/badge.svg)](https://codecov.io/github/AllenInstitute/lowcat?branch=dev)


## Description

lowcat is a collection of functions for reading, comparing, clustering, and identifying 
single-cell ATAC-seq or other similar, low-coverage chromatin accessibility data.

## Installation

lowcat depends on several Bioconductor packages. To start, install them with BiocManager:
```
if(!"BiocManager" %in% installed.packages()) {
  install.packages("BiocManager")
}

bioc_packages <- c("GenomicAlignments",
                   "GenomicRanges",
                   "Rsamtools",
                   "rtracklayer")

missing_bioc_packages <- setdiff(bioc_packages, installed.packages())

BiocManager::install(missing_bioc_packages)
```

Then, install lowcat using:
```
devtools::install_github("AllenInstitute/lowcat")
```

## Reading files

The main unit of operation for lowcat is a list of GRanges objects, one per cell/nucleus/sample.

Note that this is *not* a GRangesList object.

lowcat has helper functions for reading in many BAM files (one per sample) or multiplexed BAM files.

For multiple BAM files, the fastest option is run_pe_to_frag_parallel():
```
my_bam_dir <- "/path/to/bamfiles/"
bam_files <- list.files(my_bam_dir, 
                        pattern = ".bam$",
                        full_names = TRUE)

bam_names <- list.files(my_bam_dir, 
                        pattern = ".bam$",
                        full_names = FALSE)

bam_names <- sub(".bam$", "", bam_names)

fragment_list <- run_pe_to_frag_parallel(bam_files,
                                         sample_names = bam_names,
                                         n_cores = 4)
```

For multiplexed BAM data, like that from sci-ATAC-seq, use read_multiplexed_paired_bam().

Note 1: This requires that barcodes are present in the QNAME field of the BAM file.

Note 2: This function currently isn't very well optimized. 
It requires LOTS of available RAM - at least 4X the size of the target file itself, and can take a long time to run.
```
fragment_list <- read_multiplexed_paired_bam(bam_file,
                                             barcode_start = 1,
                                             barcode_end = 32,
                                             read_length = 50,
                                             min_frags = 100,
                                             remove_duplicates = TRUE)
```

## Downsampling and extending fragments

Once loaded, we usually evenly downsample our fragments for each cell for analysis:
```
downsample_n <- 1e4

downsampled_fragments <- downsampled_fragments(bam_fragments,
                                               downsample_n = downsample_n,
                                               discard_if_too_few = TRUE)

```

With an equal number of fragments, we next extend each fragment to an accessible region, and collapse any overlapping regions *within* each sample:
```
region_width <- 1e3

downsampled_regions <- expand_fragments(downsampled_fragments,
                                        width = region_width,
                                        collapse = TRUE)
```

## Counting overlaps between cells

After expanding to accessible regions per sample, we can compute the overlaps between every pair of cells. `lowcat` uses the `parallel` package to increase the speed of these calculations, but there are a few considerations for this algorithm.

1. This algorithm is not super well optimized. In my hands, each comparison takes ~13ms. You can estimate how long it will take using:
```
estimate_time <- function(n_samples, n_cores) {
  n_comparisons <- n_samples^2 - n_samples
  comparisons_per_core <- n_comparisons / n_cores
  total_s <- comparisons_per_core * 13 / 1000
  total_h <- round(total_s / 3600,4)
  
  cat(total_h,"hours\n")
}
```
Note that this does not include the time required to initialize paralellization and gather results.

2. If you're on Windows, you'll need to set the `cluster_type` parameter to `"PSOCK"`, and each core will need its own copy of the `downsampled_regions` list in RAM. You'll need to make sure you have enough RAM to allow these copies.

Most of the downstream analysis functions require only a matrix of Jaccard distances, which can be computed using `overlap_jaccard_mat()`.

On Windows, this will run the Jaccard overlaps in parallel:
```
jaccard_results <- overlap_jaccard_mat(downsampled_regions,
                                       n_cores = 8,
                                       cluster_type = "PSOCK")
```

On Linux or Mac, use `"FORK"` paralellization, which is more memory efficient:
```
jaccard_results <- overlap_jaccard_mat(downsampled_regions,
                                       n_cores = 8,
                                       cluster_type = "FORK")
```

If you would like more detailed results for each comparison, including the number of reads in each sample and in their intersection, use `overlap_jaccard_df()`, which returns more detailed statistics. This can later be coerced to a matrix with `jaccard_df_to_mat()`:
```
jaccard_results_df <- overlap_jaccard_df(downsampled_regions,
                                         n_cores = 8,
                                         cluster_type = "PSOCK")

jaccard_results <- jaccard_df_to_mat(jaccard_results_df)
```

## Level of Support

We are planning on occasional updating this tool with no fixed schedule. Community involvement is encouraged through both issues and pull requests.

## License

This code is released under the Allen Institute Software License. See file LICENSE for details.

## Contributions

Any contributions to this repository are subject to the Allen Institute Contribution Agreement. See file CONTRIBUTING.md for details.



