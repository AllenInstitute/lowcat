matrix_var_filter <- function(x, sd_cut = 2, top = NULL) {
  x_var <- apply(x, 1, var)
  log_var <- log(x_var + 1)
  if(!is.null(sd_cut)) {
    cutoff_filter <- which(log_var > mean(log_var) + sd_cut*sd(log_var))
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

region_fragment_count <- function(fragment_list,
                                  regions,
                                  name_col = "mcols.name") {
  out_matrix <- matrix(0,
                       nrow = length(regions),
                       ncol = length(fragment_list))
  rownames(out_matrix) <- mcols(regions)[[name_col]]
  colnames(out_matrix) <- names(fragment_list)

  for(i in 1:length(fragment_list)) {
    fragments <- fragment_list[[i]]

    region_overlaps <- countOverlaps(regions, fragments)

    out_matrix[,i] <- region_overlaps

  }

  out_matrix

}

region_fragment_count_clusters <- function(fragment_list,
                                           fragment_clusters,
                                           regions,
                                           name_col = "mcols.name") {

  clusters <- unique(fragment_clusters)

  out_matrix <- matrix(0,
                       nrow = length(regions),
                       ncol = length(clusters))
  rownames(out_matrix) <- mcols(regions)[[name_col]]
  colnames(out_matrix) <- clusters

  for(i in 1:length(clusters)) {
    cluster_fragment_list <- fragment_list[fragment_clusters == clusters[i]]

    cluster_matrix <- matrix(0,
                             nrow = length(regions),
                             ncol = length(cluster_fragment_list))
    rownames(cluster_matrix) <- mcols(regions)[[name_col]]
    colnames(cluster_matrix) <- names(cluster_fragment_list)

    for(j in 1:length(cluster_fragment_list)) {
      fragments <- cluster_fragment_list[[j]]
      region_overlaps <- countOverlaps(regions, fragments)
      cluster_matrix[,j] <- region_overlaps
    }

    out_matrix[,i] <- rowSums(cluster_matrix)

  }

  out_matrix

}

region_fragment_count_neighbors <- function(fragment_list,
                                            distance_matrix,
                                            n_neighbors = 4,
                                            max_dist = 3,
                                            regions,
                                            name_col = "mcols.name") {

  sample_names <- names(fragment_list)

  distance_matrix <- distance_matrix[sample_names,sample_names]

  out_matrix <- matrix(0,
                       nrow = length(regions),
                       ncol = length(sample_names))
  rownames(out_matrix) <- mcols(regions)[[name_col]]
  colnames(out_matrix) <- sample_names

  for(i in 1:length(sample_names)) {
    sample_name <- sample_names[i]

    sample_distances <- distance_matrix[,sample_names[i]]
    sample_distances <- sample_distances[order(sample_distances)]
    if(!is.null(max_dist)) {
      sample_distances <- sample_distances[sample_distances < max_dist]
    }

    if(length(sample_distances) > n_neighbors + 1) {
      neighbor_samples <- names(sample_distances[1:(n_neighbors + 1)])
    } else {
      neighbor_samples <- names(sample_distances)
    }


    neighbor_fragment_list <- fragment_list[neighbor_samples]

    neighbor_matrix <- matrix(0,
                              nrow = length(regions),
                              ncol = length(neighbor_fragment_list))
    rownames(neighbor_matrix) <- mcols(regions)[[name_col]]
    colnames(neighbor_matrix) <- names(neighbor_fragment_list)

    for(j in 1:length(neighbor_fragment_list)) {
      fragments <- neighbor_fragment_list[[j]]
      region_overlaps <- countOverlaps(regions, fragments)
      neighbor_matrix[,j] <- region_overlaps
    }

    out_matrix[,i] <- rowSums(neighbor_matrix)

  }

  out_matrix

}

max_column_correlation <- function(query_mat,
                                   target_mat,
                                   method = "pearson") {

  out_df <- data.frame(bam_file = colnames(query_mat),
                       max_cor  = 0,
                       cl = "",
                       stringsAsFactors = FALSE)

  for(i in 1:ncol(query_mat)) {
    query_vals <- query_mat[,i]
    cor_vals <- apply(target_mat, 2, function(x) cor(query_vals, x, method = method))
    max_cor <- max(cor_vals)
    max_name <- names(cor_vals)[cor_vals == max_cor][1]

    out_df$max_cor[i] <- max_cor
    out_df$cl[i] <- max_name
  }

  out_df
}

all_column_correlation <- function(query_mat,
                                   target_mat,
                                   method = "pearson") {

  out_mat <- matrix(0, ncol = ncol(target_mat), nrow = nrow(query_mat))

  for(i in 1:ncol(query_mat)) {
    query_vals <- query_mat[,i]
    cor_vals <- apply(target_mat, 2, function(x) cor(query_vals, x, method = method))

    out_mat[i,] <- cor_vals
  }

  out_mat

}

