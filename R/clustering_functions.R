varibow <- function(n_colors) {
  sats <- rep_len(c(0.55,0.7,0.85,1),length.out = n_colors)
  vals <- rep_len(c(1,0.8,0.6),length.out = n_colors)
  sub("FF$","",rainbow(n_colors, s = sats, v = vals))
}

better_rainbow <- function(...) {
  varibow(...)
}

pe_to_frag <- function(bamfile) {
  bam <- readGAlignmentPairs(bamfile)

  r1_bam <- first(bam)
  st_r1_bam <- start(r1_bam)
  en_r1_bam <- end(r1_bam)
  r2_bam <- last(bam)
  st_r2_bam <- start(r2_bam)
  en_r2_bam <- end(r2_bam)

  st_bam <- numeric(length(r1_bam))
  en_bam <- numeric(length(r1_bam))

  for(j in 1:length(r1_bam)) {
    st_bam[j] <- min(st_r1_bam[j],st_r2_bam[j])
    en_bam[j] <- max(en_r1_bam[j],en_r2_bam[j])
  }

  fr_bam <- GRanges(seqnames(r1_bam),IRanges(st_bam,en_bam))
  return(fr_bam)
}

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

expand_fragments <- function(fragment_list,
                             width = 1e4,
                             collapse = TRUE) {

  out_list <- vector("list", length(fragment_list))
  for(i in 1:length(fragment_list)) {

    expanded_fragments <- resize(fragment_list[[i]], width = width, fix = "center")

    if(collapse) {
      expanded_fragments <- reduce(expanded_fragments)
    }

    out_list[[i]] <- expanded_fragments
    names(out_list)[i] <- names(fragment_list)[i]
  }

  out_list

}

fragments_to_windows <- function(fragment_list,
                                 window_size = 1e4,
                                 collapse = TRUE) {

  out_list <- list()
  for(i in 1:length(fragment_list)) {

    window_df <- data.frame(chr = as.character(seqnames(fragment_list[[i]])),
                            pos = as.integer(ceiling(start(fragment_list[[i]]) / window_size))) %>%
      arrange(chr,pos)

    if(collapse) {
      window_df <- unique(window_df)
    }

    sample_list <- list()

    for(j in unique(window_df$chr)) {
      sample_list[[j]] <- window_df$pos[window_df$chr == j]
    }

    out_list[[i]] <- sample_list
    names(out_list)[i] <- names(fragment_list)[i]
  }

  out_list

}

filter_fragments <- function(fragment_list,
                             filter_GR) {

  out_list <- list()
  for(i in 1:length(fragment_list)) {
    fragments <- fragment_list[[i]]

    overlapping_fragments <- unique(subjectHits(findOverlaps(fragments, filter_GR)))

    filtered_fragments <- fragments[-overlapping_fragments]

    out_list[[i]] <- filtered_fragments
    names(out_list)[i] <- names(fragment_list)[i]
  }

  out_list
}

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

count_fragment_overlaps <- function(fragment_list,
                                    target_GRanges,
                                    binarize = TRUE,
                                    aggregate = FALSE) {

  out_mat <- matrix(nrow=length(target_GRanges),
                    ncol=length(fragment_list))
  colnames(out_mat) <- names(fragment_list)

  for(i in 1:length(fragment_list)) {

    if(i %% 100 == 0) {
      cat("\r","Counting",i,"of",length(fragment_list))
      flush.console()
    }

    if(binarize) {
      out_mat[, i] <- countOverlaps(target_GRanges,fragment_list[[i]]) > 0
    } else {
      out_mat[, i] <- countOverlaps(target_GRanges,fragment_list[[i]])
    }

  }

  if(aggregate) {
    out_vals <- rowSums(out_mat)
    out_vals
  } else {
    out_mat
  }

}

percent_display <- function(i, n_chunks) {
  paste0(floor(i/n_chunks * 100),"%")
}

time_display <- function(i, n_chunks, start_time) {
  current_time <- Sys.time()
  time_diff <- current_time - start_time
  time_num <- as.numeric(time_diff)
  time_units <- units(time_diff)

  chunks_left <- n_chunks - i
  est_time_per_chunk <- time_num/i

  est_time_remain <- est_time_per_chunk * chunks_left

  time_num <- round(time_num, 2)
  est_time_remain <- round(est_time_remain, 2)

  paste0(percent_display(i, n_chunks), " took ",time_num," ",time_units,". Estimated ",est_time_remain," ",time_units," remaining.")
}

# Parallel Jaccard clustering functions
clusterApplyLB_chunks <- function(N, n_chunks, cl, FUN, ...) {
  chunk_size <- floor(N / n_chunks)
  chunk_starts <- ((1:n_chunks) - 1) * chunk_size + 1
  chunk_ends   <- c(1:(n_chunks - 1) * chunk_size, N)

  res <- list()

  start_time <- Sys.time()

  for(i in 1:n_chunks) {
    chunk_res <- clusterApplyLB(cl, chunk_starts[i]:chunk_ends[i], FUN, ...)
    res <- c(res, chunk_res)
    cat("\r",time_display(i, n_chunks, start_time))
    flush.console()
  }

  res
}



pe_to_frag_parallel <- function(N) {
  bam_file <- bam_files[N]

  pe_to_frag(bam_file)
}

run_pe_to_frag_parallel <- function(bam_files,
                                    sample_names = NULL,
                                    n_cores = 6) {
  # Set up parallelization
  if(n_cores == "auto") {
    n_cores <- detectCores()
  }

  print(paste("Starting",n_cores,"nodes"))

  cl <- makeCluster(n_cores)

  print("Exporting necessary objects to nodes")

  clusterEvalQ(cl, library(GenomicRanges))
  clusterEvalQ(cl, library(GenomicAlignments))
  clusterExport(cl, c("bam_files",
                      "pe_to_frag",
                      "pe_to_frag_parallel"),
                # Use the function's local environment for export
                envir = environment())

  N <- length(bam_files)

  res <- clusterApplyLB_chunks(N = N,
                               n_chunks = 20,
                               cl = cl,
                               FUN = pe_to_frag_parallel)

  stopCluster(cl)

  if(is.null(sample_names)) {
    names(res) <- bam_files
  } else {
    names(res) <- sample_names
  }

  res

}

fragment_overlap_jaccard_parallel <- function(N) {
  s1 <- index_pairs[N,1]
  s2 <- index_pairs[N,2]

  ol <- sum(countOverlaps(fragment_list[[s1]],
                          fragment_list[[s2]]) > 0)

  n1 <- length(fragment_list[[s1]])
  n2 <- length(fragment_list[[s2]])

  js <- ol/(n1+n2)
  jd <- 1-js

  out_df <- data.frame(s1 = s1,
                       s2 = s2,
                       s1_name = names(fragment_list)[s1],
                       s2_name = names(fragment_list)[s2],
                       overlaps = ol,
                       n1 = n1,
                       n2 = n2,
                       jaccard_similarity = js,
                       jaccard_distance = jd)

  out_df
}

run_fragment_overlap_jaccard_parallel <- function(fragment_list,
                                                  n_cores = 6) {

  library(data.table)

  index_pairs <- as.data.frame(t(combn(1:length(fragment_list),2)))

  # Set up parallelization
  if(n_cores == "auto") {
    n_cores <- detectCores()
  }

  print(paste("Starting",n_cores,"clusters"))

  cl <- makeCluster(n_cores)

  print("Exporting necessary objects to clusters")

  clusterEvalQ(cl, library(GenomicRanges))
  clusterExport(cl, c("index_pairs","fragment_list",
                      "fragment_overlap_jaccard_parallel"),
                # Use the function's local environment for export
                envir = environment())

  N <- nrow(index_pairs)

  res <- clusterApplyLB_chunks(N = N,
                               n_chunks = 20,
                               cl = cl,
                               FUN = fragment_overlap_jaccard_parallel)

  stopCluster(cl)

  results <- rbindlist(res)


}

window_overlap_jaccard_parallel <- function(N) {
  s1 <- index_pairs[N,1]
  s2 <- index_pairs[N,2]

  ol <- 0

  windows1 <- window_list[[s1]]
  windows2 <- window_list[[s2]]

  chrs <- intersect(names(windows1),names(windows2))

  if(length(chrs) > 0) {
    for(i in 1:length(chrs)) {

      chr <- chrs[i]

      ol <- ol + length(intersect(windows1[[chr]],windows2[[chr]]))

    }
  }

  n1 <- sum(unlist(lapply(windows1,length)))
  n2 <- sum(unlist(lapply(windows2,length)))

  js <- ol/(n1+n2)
  jd <- 1-js

  out_df <- data.frame(s1 = s1,
                       s2 = s2,
                       s1_name = names(window_list)[s1],
                       s2_name = names(window_list)[s2],
                       overlaps = ol,
                       n1 = n1,
                       n2 = n2,
                       jaccard_similarity = js,
                       jaccard_distance = jd)

  out_df
}

run_window_overlap_jaccard_parallel <- function(window_list,
                                                n_cores = 6) {

  library(data.frame)

  index_pairs <- as.data.frame(t(combn(1:length(window_list),2)))

  # Set up parallelization
  if(n_cores == "auto") {
    n_cores <- detectCores()
  }

  print(paste("Starting",n_cores,"nodes"))

  cl <- makeCluster(n_cores)

  print("Exporting necessary objects to nodes")

  clusterEvalQ(cl, library(GenomicRanges))
  clusterExport(cl, c("index_pairs","window_list",
                      "window_overlap_jaccard_parallel"),
                # Use the function's local environment for export
                envir = environment())

  N <- nrow(index_pairs)

  res <- clusterApplyLB_chunks(N = N,
                               n_chunks = 20,
                               cl = cl,
                               FUN = window_overlap_jaccard_parallel)

  stopCluster(cl)

  results <- rbindlist(res)


}

res_to_distance_matrix <- function(res) {

  n_samples <- max(res$s2)

  res_names <- c(res$s1_name[1],res$s2_name[1:(n_samples - 1)])

  res_matrix <- matrix(0,ncol = n_samples, nrow = n_samples)

  for(i in 1:nrow(res)) {
    res_matrix[res$s1[i], res$s2[i]] <- res$jaccard_distance[i]
    res_matrix[res$s2[i], res$s1[i]] <- res$jaccard_distance[i]
  }

  rownames(res_matrix) <- res_names
  colnames(res_matrix) <- res_names

  res_matrix
}
