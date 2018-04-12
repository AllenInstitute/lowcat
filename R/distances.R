varibow <- function(n_colors) {
  sats <- rep_len(c(0.55,0.7,0.85,1),length.out = n_colors)
  vals <- rep_len(c(1,0.8,0.6),length.out = n_colors)
  sub("FF$","",rainbow(n_colors, s = sats, v = vals))
}

better_rainbow <- function(...) {
  varibow(...)
}



#' Linking function for computing jaccard similarity scores based on GenomicRanges overlaps in parallel.
#'
#' @param N Index of the comparison to compute
#'
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

#' Compute jaccard similarity values based on GenomicRanges overlaps in parallel.
#'
#' @param fragment_list The list object containing GenomicRanges objects.
#' @param sample_names Sample names. If NULL, will use BAM file names.
#' @param n_cores The number of cores to use in parallel. Use "auto" to detect and use all cores. Default is 6.
#'
#' @return a list of GenomicRanges objects
#'
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

  print(paste("Running",N,"comparisons."))

  res <- clusterApplyLB_chunks(N = N,
                               n_chunks = 20,
                               cl = cl,
                               FUN = fragment_overlap_jaccard_parallel)

  stopCluster(cl)

  print("Collecting results.")

  results <- rbindlist(res)

  return(results)
}

#' Linking function for computing jaccard similarity scores based on window overlaps in parallel.
#'
#' @param N Index of the comparison to compute
#'
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

#' Compute jaccard similarity values based on window overlaps in parallel.
#'
#' @param fragment_list The list object containing genomic windows.
#' @param n_cores The number of cores to use in parallel. Use "auto" to detect and use all cores. Default is 6.
#'
#' @return a list of GenomicRanges objects
#'
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

#' Convert parallel jaccard similarity results from data.frame to a matrix
#'
#' @param res The results from run_fragment_overlap_jaccard_parallel() or run_window_overlap_jaccard_parallel()
#'
#' @return a jaccard distance matrix
#'
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
