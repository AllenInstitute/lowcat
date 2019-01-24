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
#' This version returns a data.frame with extended stats for every comparison.
#'
#' @param N Index of the comparison to compute
#'
fragment_overlap_jaccard_parallel_full <- function(N) {
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
#' This version returns a data.frame with extended stats for every comparison.
#'
#' @param fragment_list The list object containing GenomicRanges objects.
#' @param n_cores The number of cores to use in parallel. Use "auto" to detect and use all cores. Default is 6.
#' @param cluster_type EXPERIMENTAL: Either "PSOCK" (for Windows) or "FORK" (possible on Linux and Mac). FORK is more memory-efficient.
#'
#' @return a data.frame with results for overlap-based comparisons.
#'
run_fragment_overlap_jaccard_parallel_full <- function(fragment_list,
                                                       n_cores = 6,
                                                       cluster_type = "PSOCK") {

  index_pairs <- as.data.frame(t(combn(1:length(fragment_list),2)))

  # Set up parallelization
  if(n_cores == "auto") {
    n_cores <- parallel::detectCores()
  }

  print(paste("Starting",n_cores,"nodes"))

  cl <- parallel::makeCluster(n_cores, cluster_type)

  print("Exporting necessary objects to nodes")

  parallel::clusterEvalQ(cl, library(GenomicRanges))
  parallel::clusterExport(cl, c("index_pairs","fragment_list",
                      "fragment_overlap_jaccard_parallel"),
                # Use the function's local environment for export
                envir = environment())

  N <- nrow(index_pairs)

  print(paste("Running",N,"comparisons."))

  res <- clusterApplyLB_chunks(N = N,
                               n_chunks = 20,
                               cl = cl,
                               FUN = fragment_overlap_jaccard_parallel_full)

  parallel::stopCluster(cl)

  print("Collecting results.")

  results <- data.table::rbindlist(res)

  return(results)
}

#' Linking function for computing jaccard similarity scores based on GenomicRanges overlaps in parallel.
#'
#' @param N Index of the comparison to compute
#'
fragment_overlap_jaccard_parallel <- function(N) {
  s1 <- index_pairs[N,1]
  s2 <- index_pairs[N,2]

  ol <- sum(countOverlaps(fragment_list[[s1]],
                          fragment_list[[s2]]))# > 0)

  n1 <- fragment_lengths[s1]
  n2 <- fragment_lengths[s2]

  jd <- 1 - (ol / (n1 + n2))

  c(jd, N)
}

#' Compute jaccard similarity values based on GenomicRanges overlaps in parallel.
#'
#' @param fragment_list The list object containing GenomicRanges objects.
#' @param n_cores The number of cores to use in parallel. Use "auto" to detect and use all cores. Default is 6.
#' @param cluster_type EXPERIMENTAL: Either "PSOCK" (for Windows) or "FORK" (possible on Linux and Mac). FORK is more memory-efficient.
#'
#' @return a matrix of jaccard distances
#'
run_fragment_overlap_jaccard_parallel <- function(fragment_list,
                                                  n_cores = 6,
                                                  cluster_type = "PSOCK") {

  index_pairs <- t(combn(1:length(fragment_list),2))
  fragment_lengths <- purrr::map_int(fragment_list, length)

  # Set up parallelization
  if(n_cores == "auto") {
    n_cores <- parallel::detectCores()
  }

  print(paste("Starting",n_cores,"nodes"))

  cl <- parallel::makeCluster(n_cores, cluster_type)

  print("Exporting necessary objects to nodes")

  parallel::clusterEvalQ(cl, library(GenomicRanges))
  parallel::clusterExport(cl, c("index_pairs","fragment_list",
                                "fragment_lengths",
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

  res <- matrix(unlist(res), ncol = 2, byrow = TRUE)
  res <- res[,1][order(res[,2])]

  res_mat <- matrix(0, length(fragment_list), length(fragment_list))
  res_mat[lower.tri(res_mat, diag = FALSE)] <- res
  res_mat <- t(res_mat)
  res_mat[lower.tri(res_mat, diag = FALSE)] <- res

  rownames(res_mat) <- colnames(res_mat) <- names(fragment_list)

  return(res_mat)
}


#' # Returns only jd
#' This version is slower.
#'
#' fragment_overlap_jaccard_parallel3 <- function(N) {
#'   ol <- findOverlaps(fragment_list[[N]],
#'                       fragment_list[(N+1):length(fragment_list)])
#'
#'   ol <- table(subjectHits(ol))
#'
#'   n1 <- fragment_lengths[N]
#'   n2 <- fragment_lengths[(N+1):length(fragment_list)]
#'
#'   jd <- 1 - ol / (n1 + n2)
#'
#'   matrix(c(jd, rep(N, length(jd))), ncol = 2)
#' }
#'
#' #' Compute jaccard similarity values based on GenomicRanges overlaps in parallel.
#' #'
#' #' @param fragment_list The list object containing GenomicRanges objects.
#' #' @param sample_names Sample names. If NULL, will use BAM file names.
#' #' @param n_cores The number of cores to use in parallel. Use "auto" to detect and use all cores. Default is 6.
#' #' @param cluster_type EXPERIMENTAL: Either "PSOCK" (for Windows) or "FORK" (possible on Linux and Mac). FORK is more memory-efficient.
#' #'
#' #' @return a list of GenomicRanges objects
#' #'
#' run_fragment_overlap_jaccard_parallel3 <- function(fragment_list,
#'                                                    n_cores = 6,
#'                                                    cluster_type = "PSOCK") {
#'
#'   fragment_list <- GRangesList(fragment_list)
#'
#'   #index_pairs <- t(combn(1:length(fragment_list),2))
#'   fragment_lengths <- unlist(lapply(fragment_list, length))
#'
#'   # Set up parallelization
#'   if(n_cores == "auto") {
#'     n_cores <- parallel::detectCores()
#'   }
#'
#'   print(paste("Starting",n_cores,"nodes"))
#'
#'   cl <- parallel::makeCluster(n_cores, cluster_type)
#'
#'   print("Exporting necessary objects to nodes")
#'
#'   parallel::clusterEvalQ(cl, library(GenomicRanges))
#'   parallel::clusterExport(cl, c("fragment_list",
#'                                 "fragment_lengths",
#'                                 "fragment_overlap_jaccard_parallel"),
#'                           # Use the function's local environment for export
#'                           envir = environment())
#'
#'   N <- length(fragment_list) - 1
#'
#'   print(paste("Running",N,"comparisons."))
#'
#'   res <- clusterApplyLB_chunks(N = N,
#'                                n_chunks = 20,
#'                                cl = cl,
#'                                FUN = fragment_overlap_jaccard_parallel3)
#'
#'   stopCluster(cl)
#'
#'   print("Collecting results.")
#'
#'   res <- do.call("rbind", res)
#'   res <- res[,1][order(res[,2])]
#'
#'   res_mat <- matrix(0, length(fragment_list), length(fragment_list))
#'   res_mat[lower.tri(res_mat, diag = FALSE)] <- res
#'   res_mat <- t(res_mat)
#'   res_mat[lower.tri(res_mat, diag = FALSE)] <- res
#'
#'   return(res_mat)
#' }

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

  library(data.table)

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

  res <- res[order(res[,1], res[,2]),]

  res_names <- c(res$s1_name[1],res$s2_name[1:(n_samples - 1)])

  res_mat <- matrix(0, nrow(res) + 1, nrow(res) + 1)
  res_mat[lower.tri(res_mat, diag = FALSE)] <- res$jd
  res_mat <- t(res_mat)
  res_mat[lower.tri(res_mat, diag = FALSE)] <- res$jd

  rownames(res_matrix) <- res_names
  colnames(res_matrix) <- res_names

  res_matrix
}
