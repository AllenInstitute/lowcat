#' Filter a matrix based on standard deviation of variance
#'
#' This filters the rows of a matrix by computing the log of variance + 1 for each row,
#' computing the standard deviation of all log_vars, and then selecting rows
#' with a log_var greater than the value of mean(log_var) + sd(log_var) * sd_cut.
#'
#' @param x The matrix to filter
#' @param sd_cut The cutoff for sd(var()). Default = 2. If NULL, performs no filtering.
#' @param top The max number of rows to return. If NULL, returns all results.
#'
#' @return a matrix with filtered rows.
#' @export
#'
matrix_var_filter <- function(x, sd_cut = 2, top = NULL) {
  x_var <- apply(x, 1, var)
  log_var <- log(x_var + 1)
  if(!is.null(sd_cut)) {
    lv_cut <- mean(log_var) + sd_cut*sd(log_var)
    cutoff_filter <- which(log_var > lv_cut)
    x <- x[cutoff_filter,]
    x_var <- x_var[cutoff_filter]
  } else { cutoff_filter <- TRUE }

  if(!is.null(top)) {
    x_var_order <- rev(order(x_var))
    if(top <= length(x_var_order)) {
      x <- x[x_var_order[1:top],]
    } else {
      print(paste("Fewer than",top,"genes were available."))
      x <- x[x_var_order,]
    }
  }

  x

}

#' Compare matrices to find the maximum correlation of columns
#'
#' For each column in the query matrix, max_column_correlation will
#' find the column in the target matrix with the highest correlation score,
#' and will return a data.frame with the correlation score and target column name.
#'
#' @param query_mat The query matrix
#' @param target_mat The target matrix
#' @param method The name of the method to pass to cor(). Can be "pearson", "kendall", or "spearman". Default is "pearson".
#'
#' @return a data.frame with 3 columns: query, max_cor, and target.
#' @export
#'
max_column_correlation <- function(query_mat,
                                   target_mat,
                                   method = "pearson") {

  out_df <- data.frame(query = colnames(query_mat),
                       max_cor  = 0,
                       target = "",
                       stringsAsFactors = FALSE)

  for(i in 1:ncol(query_mat)) {
    query_vals <- query_mat[,i]
    cor_vals <- apply(target_mat, 2, function(x) cor(query_vals, x, method = method))
    max_cor <- max(cor_vals, na.rm = TRUE)
    max_name <- names(cor_vals)[which(cor_vals == max_cor)][1]

    out_df$max_cor[i] <- max_cor
    out_df$target[i] <- max_name
  }

  out_df
}

#' Compare matrices and return a matrix with the correlations of all columns
#'
#' For each column in the query matrix, max_column_correlation will
#' find the column in the target matrix with the highest correlation score,
#' and will return a data.frame with the correlation score and target column name.
#'
#' @param query_mat The query matrix
#' @param target_mat The target matrix
#' @param method The name of the method to pass to cor(). Can be "pearson", "kendall", or "spearman". Default is "pearson".
#'
#' @return a matrix with all correlation values, with the columns matching columns of target_mat,
#' and rows matching the columns of query_mat
#' @export
#'
all_column_correlation <- function(query_mat,
                                   target_mat,
                                   method = "pearson") {

  out_mat <- matrix(0, ncol = ncol(target_mat), nrow = ncol(query_mat))

  for(i in 1:ncol(query_mat)) {
    query_vals <- query_mat[,i]
    cor_vals <- apply(target_mat, 2, function(x) cor(query_vals, x, method = method))

    out_mat[i,] <- cor_vals
  }

  out_mat

}

