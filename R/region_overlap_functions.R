#' Count the sum of overlaps for each member of a fragment list in a set of regions
#'
#' @param fragment_list The list of GenomicRanges objects to compare to regions
#' @param regions  A single GenomicRanges object as targets
#' @param name_col The column in regions that has the name of each region.
#'
#' @return A matrix with each fragment list as columns and each target region as rows containing
#' the counts of overlaps between each fragment list in each region.
#'
#' @export
#'
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

    region_overlaps <- GenomicRanges::countOverlaps(regions, fragments)

    out_matrix[,i] <- region_overlaps

  }

  out_matrix

}

#' Aggregate sample overlaps by cluster
#'
#' @param fragment_list The list of GenomicRanges objects to compare to regions
#' @param regions  A single GenomicRanges object as targets
#' @param name_col The column in regions that has the name of each region.
#'
#' @return A matrix with each fragment list as columns and each target region as rows containing
#' the counts of overlaps between each fragment list in each region.
#'
#' @export
#'
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
      region_overlaps <- GenomicRanges::countOverlaps(regions, fragments)
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
      region_overlaps <- GenomicRanges::countOverlaps(regions, fragments)
      neighbor_matrix[,j] <- region_overlaps
    }

    out_matrix[,i] <- rowSums(neighbor_matrix)

  }

  out_matrix

}
