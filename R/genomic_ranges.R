#' Combine individual GRanges from a list of GRanges objects based on group annotations.
#'
#' This is useful for making combined objects at the level of subclasses, clusters, or
#' other annotations for analysis within R.
#'
#' @param fragment_list a list of GRanges objects, each named by sample_id
#' @param samples a sample meatdata data.frame with a column called sample_id and other columns
#' to be used for grouping.
#' @param group_col The column of the samples data.frame with group assignments.
#'
#' @return A list of GRanges objects, one for each unique value in group_col.
#' @export
#'
combine_group_GRanges <- function(fragment_list,
                                  samples,
                                  group_col) {

  library(purrr)

  groups <- unique(samples[[group_col]])
  out_GRanges <- list(length(groups))

  for(i in seq_along(groups)) {
    group <- groups[i]
    group_samples <- samples$sample_id[samples[[group_col]] == group]
    group_fragments <- fragment_list[group_samples]
    if(length(group_fragments) > 1) {
      chrs <- unlist(map(group_fragments, function(x) as.character(seqnames(x))))
      starts <- unlist(map(group_fragments, start))
      ends <- unlist(map(group_fragments, end))
      strands <- unlist(map(group_fragments, function(x) as.character(strand(x))))
      group_mcols <- map(group_fragments, mcols)
      group_mcols <- do.call("rbind", group_mcols)

      group_GRanges <- GRanges(seqnames = chrs,
                               IRanges(start = starts,
                                       end = ends),
                               strand = strands,
                               mcols = group_mcols)
      out_GRanges[[i]] <- group_GRanges
    } else {
      out_GRanges[[i]] <- group_fragments
    }

  }

  names(out_GRanges) <- groups

  out_GRanges

}


#' Compare two lists of GRanges to get counts or overlap coverage.
#'
#' If mode = "counts", generates a matrix for the number of overlaps between
#' each set of query and target GRanges. Normalization parameters specify whether
#' to divide by the number of regions in the query, target, or the sum of both (total).
#'
#' If mode = "coverage", generates a matrix with the total number of overlapping bases
#' between each set of query and target GRanges. Normalization parameters specify whether
#' to divide by the total coverage of regions in the query, target, or the sum of both (total).
#'
#' @param query_regions a list of GRanges objects
#' @param target_regions a second list of GRanges objects
#' @param mode Either "counts" or "coverage".
#' @param norm Normalization mode: "none", "query","target", or "total".
#'
#' @return a matrix with query regions as rows and target regions as columns.
#' @export
#'
compare_fragment_lists <- function(query_regions,
                                  target_regions,
                                  mode = c("counts","coverage"),
                                  norm = c("none","query","target","total")) {

  mode <- match.arg(mode, c("counts","coverage"))
  norm <- match.arg(norm, c("none","query","target","total"))

  queries <- names(query_regions)
  targets <- names(target_regions)

  scores <- matrix(0, nrow = length(queries), ncol = length(targets))
  rownames(scores) <- queries
  colnames(scores) <- targets

  for(t in seq_along(targets)) {
    for(q in seq_along(queries)) {
      if(mode == "counts") {
        scores[q,t] <- sum(countOverlaps(query_regions[[q]], target_regions[[t]]))
        if(norm == "query") {
          scores[q,t] <- scores[q,t] / length(query_regions[[q]])
        } else if(norm == "target") {
          scores[q,t] <- scores[q,t] / length(target_regions[[t]])
        } else if(norm == "total") {
          scores[q,t] <- scores[q,t] / (length(query_regions[[q]] + length(target_regions[[t]])))
        }
      } else if(mode == "coverage") {
        scores[q,t] <- sum(width(intersect(query_regions[[q]], target_regions[[t]], ignore.strand = TRUE)))
        if(norm == "query") {
          scores[q,t] <- scores[q,t] / sum(width(query_regions[[q]]))
        } else if(norm == "target") {
          scores[q,t] <- scores[q,t] / sum(width(length(target_regions[[t]])))
        } else if(norm == "total") {
          scores[q,t] <- scores[q,t] / (sum(width(query_regions[[q]])) + sum(width(target_regions[[t]])))
        }
      }


    }
  }

  scores

}
