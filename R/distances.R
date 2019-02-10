#' Linking function for computing jaccard similarity scores based on GenomicRanges overlaps in parallel.
#'
#' This version returns a data.frame with extended stats for every comparison.
#'
#' @param N Index of the comparison to compute
#' @export
#'
overlap_jaccard_df_fun <- function(N) {
  s1 <- index_pairs[N,1]
  s2 <- index_pairs[N,2]

  ol <- sum(GenomicRanges::countOverlaps(fragment_list[[s1]],
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
#' For a more streamlined version that returns only a matrix of distances, see ?overlap_jaccard_mat.
#'
#' @param fragment_list The list object containing GenomicRanges objects.
#' @param n_cores The number of cores to use in parallel. Use "auto" to detect and use all cores. Default is "auto".
#' @param cluster_type Either "PSOCK" (for Windows) or "FORK" (possible on Linux and Mac). FORK is more memory-efficient.
#'
#' @return a data.frame with results for overlap-based comparisons:
#'
#' @export
overlap_jaccard_df <- function(fragment_list,
                               n_cores = "auto",
                               cluster_type = "PSOCK") {

  index_pairs <- as.data.frame(t(combn(1:length(fragment_list),2)))
  N <- nrow(index_pairs)

  # Set up parallelization
  if(n_cores == "auto") {
    n_cores <- parallel::detectCores()
  }

  if(n_cores == 1) {
    res <- purrr::map(1:N,
                      overlap_jaccard_df_fun)
  } else {
    print(paste("Starting",n_cores,"nodes"))

    cl <- parallel::makeCluster(n_cores, cluster_type)

    print("Exporting necessary objects to nodes")

    parallel::clusterEvalQ(cl, library(GenomicRanges))
    parallel::clusterExport(cl, c("index_pairs","fragment_list",
                                  "fragment_overlap_jaccard_parallel"),
                            # Use the function's local environment for export
                            envir = environment())


    print(paste("Running",N,"comparisons."))

    res <- clusterApplyLB_chunks(N = N,
                                 n_chunks = 20,
                                 cl = cl,
                                 FUN = overlap_jaccard_df_fun)

    parallel::stopCluster(cl)


  }
  print("Collecting results.")

  results <- data.table::rbindlist(res)

  return(results)
}

#' For compatibility with older scripts
#'
#' See ?overlap_jaccard_df for current function.
#' @export
run_fragment_overlap_jaccard_parallel_full <- overlap_jaccard_df

#' Linking function for computing jaccard similarity scores based on GenomicRanges overlaps in parallel.
#'
#' @param N Index of the comparison to compute
#' @export
overlap_jaccard_mat_fun <- function(N) {
  s1 <- index_pairs[N,1]
  s2 <- index_pairs[N,2]

  ol <- sum(GenomicRanges::countOverlaps(fragment_list[[s1]],
                                         fragment_list[[s2]]) > 0)

  n1 <- fragment_lengths[s1]
  n2 <- fragment_lengths[s2]

  jd <- 1 - (ol / (n1 + n2))

  c(jd, N)
}

#' Compute jaccard similarity values based on GenomicRanges overlaps in parallel.
#'
#' @param fragment_list The list object containing GenomicRanges objects.
#' @param n_cores The number of cores to use in parallel. Use "auto" to detect and use all cores. Default is 6.
#' @param cluster_type Either "PSOCK" (for Windows) or "FORK" (possible on Linux and Mac). FORK is more memory-efficient.
#'
#' @return a matrix of jaccard distances
#' @export
overlap_jaccard_mat <- function(fragment_list,
                                n_cores = "auto",
                                cluster_type = "PSOCK") {

  index_pairs <- t(combn(1:length(fragment_list),2))
  fragment_lengths <- purrr::map_int(fragment_list, length)
  N <- nrow(index_pairs)

  # Set up parallelization
  if(n_cores == "auto") {
    n_cores <- parallel::detectCores()
  }

  if(n_cores == 1) {
    res <- purrr::map(1:N,
                      overlap_jaccard_mat_fun)
  } else {

    print(paste("Starting",n_cores,"nodes"))

    cl <- parallel::makeCluster(n_cores, cluster_type)

    print("Exporting necessary objects to nodes")

    parallel::clusterEvalQ(cl, library(GenomicRanges))
    parallel::clusterExport(cl, c("index_pairs","fragment_list",
                                  "fragment_lengths",
                                  "fragment_overlap_jaccard_parallel"),
                            # Use the function's local environment for export
                            envir = environment())


    print(paste("Running",N,"comparisons."))

    res <- clusterApplyLB_chunks(N = N,
                                 n_chunks = 20,
                                 cl = cl,
                                 FUN = overlap_jaccard_mat_fun)

    stopCluster(cl)
  }


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

#' For compatibility with older scripts
#'
#' See ?overlap_jaccard_mat for current function.
#' @export
run_fragment_overlap_jaccard_parallel <- overlap_jaccard_mat

#' Linking function for computing jaccard similarity scores based on window overlaps in parallel.
#'
#' @param N Index of the comparison to compute
#' @export
window_jaccard_df_fun <- function(N) {
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
#' @param cluster_type Either "PSOCK" (for Windows) or "FORK" (possible on Linux and Mac). FORK is more memory-efficient.
#'
#' @return a list of GenomicRanges objects
#' @export
window_jaccard_df <- function(window_list,
                              n_cores = "auto",
                              cluster_type = "PSOCK") {

  index_pairs <- as.data.frame(t(combn(1:length(window_list),2)))
  N <- nrow(index_pairs)

  # Set up parallelization
  if(n_cores == "auto") {
    n_cores <- parallel::detectCores()
  }

  if(n_cores == 1) {
    res <- purrr::map(1:N,
                      window_jaccard_df_fun)
  } else {

    print(paste("Starting",n_cores,"nodes"))

    cl <- parallel::makeCluster(n_cores, cluster_type)

    print("Exporting necessary objects to nodes")

    parallel::clusterExport(cl,
                            c("index_pairs","window_list",
                              "window_overlap_jaccard_parallel"),
                            # Use the function's local environment for export
                            envir = environment())


    res <- clusterApplyLB_chunks(N = N,
                                 n_chunks = 20,
                                 cl = cl,
                                 FUN = window_jaccard_df_fun)

    parallel::stopCluster(cl)

  }

  results <- data.table::rbindlist(res)

}

#' For compatibility with older scripts
#'
#' See ?window_jaccard_df for current function.
#' @export
run_window_overlap_jaccard_parallel <- window_jaccard_df

#' Convert jaccard similarity results from data.frame-based functions to a matrix
#'
#' @param res The results from overlap_jaccard_df() or window_jaccard_df()
#'
#' @return a jaccard distance matrix
#' @export
res_to_distance_matrix <- function(res) {

  res <- res[order(res[,1], res[,2]),]

  n_samples <- max(res$s2[res$s1 == 1])
  res_names <- c(res$s1_name[1],res$s2_name[1:(n_samples - 1)])

  res_mat <- matrix(0, n_samples, n_samples)
  res_mat[lower.tri(res_mat, diag = FALSE)] <- res$jaccard_distance
  res_mat <- t(res_mat)
  res_mat[lower.tri(res_mat, diag = FALSE)] <- res$jaccard_distance

  rownames(res_mat) <- res_names
  colnames(res_mat) <- res_names

  res_mat
}
