#' Read a paired-end BAM file as GRanges fragments
#'
#' @param bamfile The BAM file to read
#'
#' @return a GenomicRanges object
#' @export
pe_to_frag <- function(bamfile) {
  bam <- GenomicAlignments::readGAlignmentPairs(bamfile)
  if(length(bam) > 0) {
    r1_bam <- GenomicAlignments::first(bam)
    st_r1_bam <- start(r1_bam)
    en_r1_bam <- end(r1_bam)
    r2_bam <- GenomicAlignments::last(bam)
    st_r2_bam <- start(r2_bam)
    en_r2_bam <- end(r2_bam)

    st_bam <- numeric(length(r1_bam))
    en_bam <- numeric(length(r1_bam))

    for(j in 1:length(r1_bam)) {
      st_bam[j] <- min(st_r1_bam[j],st_r2_bam[j])
      en_bam[j] <- max(en_r1_bam[j],en_r2_bam[j])
    }

    fr_bam <- GenomicRanges::GRanges(seqnames(r1_bam), IRanges::IRanges(st_bam,en_bam))
  } else {
    fr_bam <- GenomicRanges::GRanges(bam)
  }

  return(fr_bam)
}

#' Read a single-end BAM file as GRanges cut sites based on the 5' end of each read
#'
#' @param bamfile The BAM file to read
#'
#' @return a GenomicRanges object
#' @export
se_to_cuts <- function(bamfile) {
  bam <- GenomicAlignments::readGAlignments(bamfile)
  gr <- GenomicRanges::GRanges(bam)
  cuts <- GenomicRanges::resize(gr, 1, fix = "start")
  return(cuts)
}

#' Read multiple paired-end BAM files as a list of GRanges objects
#'
#' @param bamfiles A character vector of file locations for the BAM files to read
#'
#' @return a list object containing GenomicRanges objects. List names will BAM filenames without leading directories.
#' @export
bam_to_fragment_list <- function(bamfiles) {
  fragment_list <- vector("list",length(bamfiles))

  n_files <- length(bamfiles)

  for(i in 1:n_files) {

    if(i == 1 | i %% 100 == 0) {
      cat("\r","Reading ",i," of ",n_files)
      flush.console()
    }

    fragment_list[[i]] <- pe_to_frag(bamfiles[i])

    names(fragment_list)[i] <- sub(".+/","",bamfiles[i])
  }

  fragment_list
}

#' Read a multiplexed, paired-end BAM file and convert to a GRanges fragment list
#'
#' This requires that the barcodes for multiplexing are in the QNAME/read name field of the BAM file.
#'
#' @param bam A character vector of file locations for the BAM file
#' @param barcode_start A numeric value indicating where the barcode begins in QNAME (default = 1)
#' @param barcode_end A numeric value indicating where the barcode ends in QNAME (default = 32)
#' @param read_length A numeric value indicating the read length (default = 50)
#' @param min_frags A numeric value indicating the minimum number of fragments required to retain a barcode (default = 100)
#' @param remove_duplicates A logical value indicating whether or not to deduplicate the fragments after demultiplexing (default = FALSE)
#'
#' @return a list object containing GenomicRanges objects. List names will be barcodes.
#' @export
read_multiplexed_paired_bam <- function(bam,
                                        barcode_start = 1,
                                        barcode_end = 32,
                                        read_length = 50,
                                        min_frags = 100,
                                        remove_duplicates = TRUE) {

  param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isPaired = TRUE,
                                                                 isProperPair = TRUE,
                                                                 isFirstMateRead = TRUE),
                                   what = c("qname", "rname", "pos", "isize"))

  print("Reading multiplexed BAM file")

  bam_values <- data.frame(Rsamtools::scanBam(bam, param = param)[[1]])

  print("Identifying barcodes and correcting positions")



  if(min_frags > 1) {
    bam_values$barcode <- substr(bam_values$qname, barcode_start, barcode_end)
    bc_counts <- table(bam_values$barcode)
    keep_bc <- names(bc_counts)[which(bc_counts > min_frags)]
    keep_lgl <- bam_values$barcode %in% keep_bc
    bam_values <- bam_values[keep_lgl,]
    rm(keep_lgl)
    rm(keep_bc)
    rm(bc_counts)
    bam_values$barcode <- as.factor(bam_values$barcode)
  } else {
    bam_values$barcode <- as.factor(substr(bam_values$qname, barcode_start, barcode_end))
    bam_values <- bam_values[,-1]
  }

  bam_values$start <- 0

  bam_values$start[bam_values$isize >= 0] <- bam_values$pos[bam_values$isize >= 0]
  bam_values$start[bam_values$isize < 0] <- bam_values$pos[bam_values$isize < 0] + bam_values$isize[bam_values$isize < 0]

  bam_values$end <- 0

  bam_values$end[bam_values$isize >= 0] <- bam_values$pos[bam_values$isize >= 0] + bam_values$isize[bam_values$isize >= 0] + read_length
  bam_values$end[bam_values$isize < 0] <- bam_values$pos[bam_values$isize < 0] + read_length

  bam_values <- bam_values[,c("rname","start","end","barcode")]
  names(bam_values)[1] <- "chr"

  print(paste("Splitting based on", length(levels(bam_values$barcode)), "barcodes"))

  split_dfs <- split(bam_values, bam_values$barcode)

  names(split_dfs) <- map_chr(split_dfs,
                              function(x) {
                                as.character(x$barcode[1])
                              })

  print(head(names(split_dfs)))

  to_gr <- function(x) {
    GenomicRanges::GRanges(seqnames = x$chr,
                           IRanges::IRanges(start = x$start,
                                            end = x$end))
  }

  possible_GRanges <- purrr::possibly(to_gr, NULL)

  print("Converting to GenomicRanges")
  bam_regions <- purrr::map(1:length(split_dfs),
                            function(i) {
                              x <- split_dfs[[i]]
                              gr <- possible_GRanges(x)

                              if(remove_duplicates) {
                                GenomicRanges::unique(gr)
                              } else {
                                gr
                              }
                            })

  names(bam_regions) <- names(split_dfs)

  return(bam_regions)

}

#' Equally downsample a list of GenomicRanges objects
#'
#' @param fragment_list The list object containing GenomicRanges objects
#' @param downsample_n The number of fragments to retain. Default is 1e4.
#' @param discard_if_too_few Logical to determine if samples with fewer than downsample_n regions should be discarded.
#' If FALSE, samples with < downample_n regions will be retained without sampling. Default is TRUE.
#' @param seed A seed value for sampling for reproducibility. Default is 42.
#'
#' @return A list object containing downsampled GenomicRanges objects.
#' @export
downsample_fragments <- function(fragment_list,
                                 downsample_n = 1e4,
                                 discard_if_too_few = TRUE,
                                 seed = 42) {

  out_list <- vector("list",length(fragment_list))

  fragment_counts <- lapply(fragment_list, length)

  if(downsample_n == "auto") {
    downsample_n <- min(unlist(fragment_counts))
  }

  fragment_gt_downsample <- which(fragment_counts >= downsample_n)

  print(paste("Downsampling to",downsample_n,".",length(fragment_gt_downsample),"of",length(fragment_list),"have >",downsample_n,"."))

  set.seed(seed)

  for(i in fragment_gt_downsample) {
    out_list[[i]] <- sample(fragment_list[[i]], downsample_n)
  }

  names(out_list) <- names(fragment_list)

  if(discard_if_too_few) {
    out_list <- out_list[fragment_gt_downsample]
  }

  out_list

}

#' Expand the boundaries of regions in a list of GenomicRanges objects
#'
#' This function expands from the center of the GenomicRanges objects using GenomicRanges::resize(fix = "center")
#'
#' @param fragment_list The list object containing GenomicRanges objects
#' @param width The fragment size to output. Default is 1e4.
#' @param collapse Logical to determine if overlapping fragments should be merged after expansion. Default is TRUE.
#'
#' @return A list object containing resized GenomicRanges objects.
#' @export
expand_fragments <- function(fragment_list,
                             width = 1e4,
                             collapse = TRUE) {

  out_list <- vector("list", length(fragment_list))
  for(i in 1:length(fragment_list)) {

    expanded_fragments <- GenomicRanges::resize(fragment_list[[i]], width = width, fix = "center")

    if(collapse) {
      expanded_fragments <- GenomicRanges::reduce(expanded_fragments)
    }

    out_list[[i]] <- expanded_fragments
    names(out_list)[i] <- names(fragment_list)[i]
  }

  out_list

}

#' Convert a list of GenomicRanges objects to a list of genomic window ids
#'
#' This function converts genomic positions to integer window positions based on the 5' end of each GenomicRanges region.
#'
#' @param fragment_list The list object containing GenomicRanges objects
#' @param window_size The size of the genomic windows to use. Default is 1e4.
#' @param collapse Logical to determine if fragments found in the same window should be merged after conversion. Default is TRUE.
#'
#' @return A nested list object. Each sample in the fragment list will result in a list containing vectors with integer values for windows on each chromosome.
#' @export
fragments_to_windows <- function(fragment_list,
                                 window_size = 1e4,
                                 collapse = TRUE) {

  out_list <- list()
  for(i in 1:length(fragment_list)) {

    window_df <- data.frame(chr = as.character(GenomicRanges::seqnames(fragment_list[[i]])),
                            pos = as.integer(ceiling(start(fragment_list[[i]]) / window_size))) %>%
      dplyr::arrange(chr,pos)

    if(collapse) {
      window_df <- unique(window_df)
    }

    sample_list <- list()

    for(j in unique(window_df$chr)) {
      sample_list[[j]] <- window_df$pos[window_df$chr == j]
    }

    names(sample_list) <- unique(window_df$chr)

    out_list[[i]] <- sample_list
    names(out_list)[i] <- names(fragment_list)[i]
  }

  out_list

}

#' Filter a list of GenomicRanges objects against a single GenomicRanges object
#'
#' @param fragment_list The list object containing GenomicRanges objects
#' @param filter_GR The GenomicRanges object to use for filtering.
#' @param mode Whether to "remove" or "keep" overlaps. Default is "remove"
#' @param ignore_strand Logical, whether or not to ignore the strand of regions in the comparison
#'
#' @return A list object containing filtered GenomicRanges objects.
#' @export
filter_fragments <- function(fragment_list,
                             filter_GR,
                             mode = "remove",
                             ignore_strand = TRUE) {

  out_list <- list()
  for(i in 1:length(fragment_list)) {
    fragments <- fragment_list[[i]]

    overlapping_fragments <- unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(fragments,
                                                                                     filter_GR,
                                                                                     ignore.strand = ignore_strand)))

    if(mode == "remove") {
      filtered_fragments <- fragments[-overlapping_fragments]
    } else if(mode == "keep") {
      filtered_fragments <- fragments[overlapping_fragments]
    }

    out_list[[i]] <- filtered_fragments
    names(out_list)[i] <- names(fragment_list)[i]
  }

  out_list
}

#' Collapse a list of GenomicRanges objects to a single, merged GenomicRanges object
#'
#' @param fragment_list The list object containing GenomicRanges objects
#' @param width If not NULL, whill resize the regions in each GenomicRanges object to this width before merging. Default is 500.
#'
#' @return A GenomicRanges object containing merged regions.
#' @export
collapse_fragment_list <- function(fragment_list,
                                   width = 500) {

  if(!is.null(width)) {
    fragment_list <- expand_fragments(fragment_list,
                                      width = width,
                                      collapse = TRUE)
  } else {
    fragment_list <- lapply(fragment_list, reduce)
  }

  out_GRanges <- fragment_list[[1]]

  for(i in 2:length(fragment_list)) {

    merged_ranges <- c(out_GRanges, fragment_list[[i]])
    out_GRanges <- reduce(merged_ranges)

  }

  out_GRanges

}


#' Linking function for running pe_to_frag in parallel mode.
#'
#' @param N Index of the bam_file to read
#' @export
pe_to_frag_parallel <- function(N) {
  bam_file <- bam_files[N]

  pe_to_frag(bam_file)
}

#' Read paired-end BAM files to GenomicRanges objects in parallel
#'
#' @param bam_files a vector of bam file locations
#' @param sample_names Sample names. If NULL, will use BAM file names.
#' @param n_cores The number of cores to use in parallel. Use "auto" to detect and use all cores. Default is 6.
#'
#' @return a list of GenomicRanges objects
#' @export
run_pe_to_frag_parallel <- function(bam_files,
                                    sample_names = NULL,
                                    n_cores = 6) {
  # Set up parallelization
  if(n_cores == "auto") {
    n_cores <- parallel::detectCores()
  }

  print(paste("Starting",n_cores,"nodes"))

  cl <- parallel::makeCluster(n_cores)

  print("Exporting necessary objects to nodes")

  parallel::clusterEvalQ(cl, library(GenomicRanges))
  parallel::clusterEvalQ(cl, library(GenomicAlignments))
  parallel::clusterExport(cl, c("bam_files",
                                "pe_to_frag",
                                "pe_to_frag_parallel"),
                          # Use the function's local environment for export
                          envir = environment())

  N <- length(bam_files)

  res <- clusterApplyLB_chunks(N = N,
                               n_chunks = 20,
                               cl = cl,
                               FUN = pe_to_frag_parallel)

  parallel::stopCluster(cl)

  if(is.null(sample_names)) {
    names(res) <- bam_files
  } else {
    names(res) <- sample_names
  }

  res

}

#' Linking function for running se_to_cuts in parallel mode.
#'
#' @param N Index of the bam_file to read
#' @export
se_to_cuts_parallel <- function(N) {
  bam_file <- bam_files[N]

  se_to_cuts(bam_file)
}

#' Read single-end BAM files to GenomicRanges objects in parallel
#'
#' @param bam_files a vector of bam file locations
#' @param sample_names Sample names. If NULL, will use BAM file names.
#' @param n_cores The number of cores to use in parallel. Use "auto" to detect and use all cores. Default is 6.
#'
#' @return a list of GenomicRanges objects
#' @export
run_se_to_cuts_parallel <- function(bam_files,
                                    sample_names = NULL,
                                    n_cores = 6) {
  # Set up parallelization
  if(n_cores == "auto") {
    n_cores <- parallel::detectCores()
  }

  print(paste("Starting",n_cores,"nodes"))

  cl <- parallel::makeCluster(n_cores)

  print("Exporting necessary objects to nodes")

  parallel::clusterEvalQ(cl, library(GenomicRanges))
  parallel::clusterEvalQ(cl, library(GenomicAlignments))
  parallel::clusterExport(cl, c("bam_files",
                                "se_to_cuts",
                                "se_to_cuts_parallel"),
                          # Use the function's local environment for export
                          envir = environment())

  N <- length(bam_files)

  res <- clusterApplyLB_chunks(N = N,
                               n_chunks = 20,
                               cl = cl,
                               FUN = se_to_cuts_parallel)

  stopCluster(cl)

  if(is.null(sample_names)) {
    names(res) <- bam_files
  } else {
    names(res) <- sample_names
  }

  res

}

#' Merge multiple BAM files to a single output BAM file
#'
#' @param bam_files a character vector listing all of the BAM files to merge
#' @param out_file the target BAM file to write to
#' @param make_indexes whether or not to make indexes for the input files. Default is FALSE
#' @param sort_out_file whether or not to sort the output file. default is TRUE, which will force index_out_file to TRUE.
#' @param index_out_file whether or not to index the output file. Default is TRUE.
#' @param keep_unsorted If sort_out_file == TRUE, whether or not to keep the original, unsorted file. Default is FALSE.
#'
#' @export
merge_bam_files <- function(bam_files,
                            out_file,
                            make_indexes = FALSE,
                            sort_out_file = TRUE,
                            index_out_file = TRUE,
                            keep_unsorted = FALSE) {

  library(rbamtools)

  l_exist <- file.exists(bam_files)
  n_exist <- sum(l_exist)
  n_missing <- sum(!l_exist)

  l_indexed <- file.exists(paste0(bam_files,".bai"))
  n_indexed <- sum(l_indexed)
  n_no_index <- sum(!l_indexed)

  # Check for files and indexes
  if(n_exist < length(bam_files)) {
    stop(paste0("Can't find", n_missing ," bam files: ", paste(bam_files[!l_exist], collapse = ", ")))
  }
  if(n_indexed < length(bam_files)) {
    if(make_indexes) {
      print("Missing ", n_no_index, " indexes. These will be generated.")
    } else {
      stop(paste0("Can't find", n_no_index ," bam indexes (bam.bai): ", paste(paste0(bam_files[!l_indexed],".bai"), collapse = ", "), ". Use make_indexes = TRUE to generate them on the fly."))
    }
  }

  # Get header from first file
  bam1 <- bamReader(bam_files[1])
  header <- getHeader(bam1)
  bamClose(bam1)

  # Initialize the new file
  out <- bamWriter(header, out_file)

  for(i in 1:length(bam_files)) {
    # Open each bam
    bam <- bamReader(bam_files[i])

    # If an index doesn't exist, make one
    index_file <- paste0(bam_files[i],".bai")
    if(!l_indexed[i]) {
      createIndex(bam, index_file)
    }

    # Load Index
    loadIndex(bam, index_file)

    # Get all chromosomes
    refs <- getRefData(bam)
    for(j in 1:nrow(refs)) {
      # Get coords for each chromosome
      coords <- getRefCoords(bam, refs$SN[j])
      # Convert to bamRange
      range <- bamRange(bam, coords)
      # Write each chr to out.
      bamSave(out, range)
    }

    bamClose(bam)
  }

  bamClose(out)

  if(sort_out_file | index_out_file) {
    print(paste0("Sorting ",out_file,". (Required for indexing)"))
    out <- bamReader(out_file)
    sort_prefix <- sub(".bam$",".srt",out_file)
    sort_prefix <- sub(".+/","",sort_prefix)

    bamSort(out, prefix = sort_prefix)
    bamClose(out)

    if(!keep_unsorted) {
      file.remove(out_file)
      file.rename(sub(".bam$",".srt.bam",out_file), out_file)
    }

    if(index_out_file) {
      if(keep_unsorted) {
        out <- bamReader(sub(".bam$",".srt.bam",out_file))
        createIndex(out, paste0(sub(".bam$",".srt.bam",out_file),".bai"))
        bamClose(out)
      } else {
        out <- bamReader(out_file)
        createIndex(out, paste0(out_file,".bai"))
        bamClose(out)
      }

    }

  }

}

#' Downsample fragments within cluster sets to match the sample with the lowest number of fragments
#'
#' @param fragment_list The list object containing GenomicRanges objects.
#' @param clusters a vector with cluster assignments for each item in fragment_list
#'
#' @return a lis of GenomicRanges objects with all members of each cluster downsampled to the
#' minimum number of reads of all cluster members.
#' @export
#'
balance_fragment_clusters <- function(fragment_list,
                                      clusters) {
  unique_clusters <- unique(clusters)

  balanced_fragments <- map(clusters,
                            function(this_cluster) {
                              cluster_fragments <- fragment_list[clusters == this_cluster]
                              min_fragments <- min(map_int(cluster_fragments, length))
                              cluster_balanced_fragments <- downsample_fragments(cluster_fragments,
                                                                                 downsample_n = min_fragments)
                              cluster_balanced_fragments
                            })

  balanced_fragments <- unlist(balanced_fragments, recursive = F)
  balanced_fragments <- balanced_fragments[names(fragment_list)]

}


#' Convert data.frames in BED-like format to GRanges objects
#'
#' @param bed A data.frame with columns chr, start, end, name, strand, and score.
#'
#' @return A GenomicRanges object
#' @export
#'
bed_to_GRanges <- function(bed) {

  gr <- GenomicRanges::GRanges(seqnames=bed$chr,
                               IRanges::IRanges(start=bed$start,
                                                end=bed$end),
                               strand=bed$strand,
                               mcols=bed[,c("name","score")])

  return(gr)
}

#' Convert UCSC-style genomic locations to a BED-like data.frame
#'
#' @param gr A GenomicRanges object
#'
#' @return a data.frame with chr, start, end, and strand.
#' @export
#'
GRanges_to_bed <- function(gr) {
  bed <- data.frame(chr = GenomicRanges::seqnames(gr),
                    start = GenomicRanges::start(gr),
                    end = GenomicRanges::end(gr),
                    strand = GenomicRanges::strand(gr))

  bed$chr <- as.character(bed$chr)
  bed$strand <- as.character(bed$strand)

  bed
}

#' Convert UCSC-style genomic locations to a GenomicRegions object.
#'
#' @param ucsc_loc A vector of UCSC-like genomic locations. Example: "chr1:152,548,974-152,550,854"
#'
#' @return a GenomicRanges object for the UCSC locations
#' @export
#'
ucsc_loc_to_GRanges <- function(ucsc_loc) {
  chr <- sub(":.+","",ucsc_loc)
  start_pos <- sub(".+:","",sub("-.+","",ucsc_loc))
  start_pos <- as.numeric(gsub(",","",start_pos))
  end_pos <- sub(".+-","",ucsc_loc)
  end_pos <- as.numeric(gsub(",","",end_pos))
  GenomicRanges::GRanges(seqnames = chr,
                         IRanges::IRanges(start_pos,end_pos),
                         strand = "+")
}

#' Convert GenomicRanges objects to a vector of UCSC browser locations.
#'
#' @param gr A GenomicRanges object
#'
#' @return a character vector of UCSC browser locations
#' @export
#'
GRanges_to_ucsc_loc <- function(gr) {
  paste0(GenomicRanges::seqnames(gr), ":",
         prettyNum(GenomicRanges::start(gr), big.mark = ","), "-",
         prettyNum(GenomicRanges::end(gr), big.mark = ","))
}

#' Convert GenomicRanges to a bed-like data.frame
#'
#' @param ucsc_loc A vector of UCSC-like genomic locations. Example: "chr1:152,548,974-152,550,854"
#'
#' @return a BED-like data.frame with chr, start, end, and strand (+)
#' @export
#'
ucsc_loc_to_bed <- function(ucsc_loc) {
  chr <- sub(":.+","",ucsc_loc)
  start_pos <- sub(".+:","",sub("-.+","",ucsc_loc))
  start_pos <- as.numeric(gsub(",","",start_pos))
  end_pos <- sub(".+-","",ucsc_loc)
  end_pos <- as.numeric(gsub(",","",end_pos))
  data.frame(chr = chr,
             start = start_pos,
             end = end_pos,
             strand = "+")
}

#' Convert GenomicRanges objects to a vector of UCSC browser locations.
#'
#' @param bed A data.frame with columns chr, start, end, name, strand, and score.
#'
#' @return a character vector of UCSC browser locations
#' @export
#'
bed_to_ucsc_loc <- function(bed) {
  paste0(bed$chr, ":",
         prettyNum(start(bed$start), big.mark = ","), "-",
         prettyNum(end(bed$end), big.mark = ","))
}
