# lowcat: Low Coverage Accessibility and Transcriptomics

<img src="https://upload.wikimedia.org/wikipedia/commons/f/f7/Rusty-spotted_cat_%28f._rubiginosa%29.JPG" alt="Rusty-spotted lowcat" width="300px"/>

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
                   "rsamtools",
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
## Level of Support

We are planning on occasional updating this tool with no fixed schedule. Community involvement is encouraged through both issues and pull requests.

## License

This code is released under the Allen Institute Software License. See file LICENSE for details.

## Contributions

Any contributions to this repository are subject to the Allen Institute Contribution Agreement. See file CONTRIBUTING.md for details.



