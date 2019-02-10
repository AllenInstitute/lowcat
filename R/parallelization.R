#' Convert i and n_chunks into a percentage character object.
#'
#' @param i The current value
#' @param n_chunks number of chunks
#'
#' @return a character object.
#' @export
percent_display <- function(i, n_chunks) {
  paste0(floor(i/n_chunks * 100),"%")
}

#' Display time elapsed and projected completion time.
#'
#' @param i The current chunk under analysis
#' @param n_chunks number of total chunks
#' @param start_time The start time of the first chunk
#'
#' @return a character object.
#' @export
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

  paste0(percent_display(i, n_chunks), " took ",time_num," ",time_units,". Estimated ",est_time_remain," ",time_units," remaining.          ")
}

#' cat for replaceable messages
#'
#' @param ... Passed to cat()
#'
#' @export
cat_update <- function(...) {
  cat("\r", ...)
  flush.console()
}

#' Run clusterApplyLB across a single variable, N, using chunks with user feedback.
#'
#' @param N Variable passed to nodes running FUN
#' @param n_chunks Number of chunks to use. Default = 20 (5% of N per chunk).
#' @param cl A cluster object created by the parallel or snow package.
#' @param FUN The function to run using clusterApplyLB()
#' @param ... Additional parameters passed to clusterApplyLB()
#'
#' @return a list object with the outputs of FUN
#' @export
clusterApplyLB_chunks <- function(N,
                                  n_chunks = 20,
                                  cl, FUN, ...) {
  chunk_size <- floor(N / n_chunks)
  chunk_starts <- ((1:n_chunks) - 1) * chunk_size + 1
  chunk_ends   <- c(1:(n_chunks - 1) * chunk_size, N)

  res <- list()

  start_time <- Sys.time()

  for(i in 1:n_chunks) {
    chunk_res <- parallel::clusterApplyLB(cl, chunk_starts[i]:chunk_ends[i], FUN, ...)
    res <- c(res, chunk_res)
    cat("\r",time_display(i, n_chunks, start_time))
    flush.console()
  }

  res
}

