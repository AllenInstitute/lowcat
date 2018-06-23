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


sum_GRanges_coverage <- function(bam_fragments,
                                 target_regions) {



}
