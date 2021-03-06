% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{build_pile_plot}
\alias{build_pile_plot}
\title{Build a multi-track pileup plot}
\usage{
build_pile_plot(gr_list, ucsc_loc, highlight_loc = NULL,
  padding = c(1e+05, 1e+05), gr_groups = NULL, group_colors = NULL,
  norm = c("PM", "max"), max_val = NULL, window_size = NULL,
  window_mode = c("max", "mean", "median"), target_color = "#B7B7B7",
  highlight_color = "#F9ED32")
}
\arguments{
\item{gr_list}{A list of GRanges objects}

\item{ucsc_loc}{A target location, in UCSC format (e.g. "chr1:533,235-552,687)}

\item{highlight_loc}{A location to use for highlights in UCSC format}

\item{padding}{A 2-element numeric vector with upstream and downstream padding around the ucsc_loc to extend the plotting window.}

\item{gr_groups}{A vector indicating which group each GRanges object in gr_list belongs to. If NULL, will group all samples into a single track.}

\item{group_colors}{A named vector, one per group, with colors for each group. Names should match values in gr_groups.}

\item{norm}{Normalization for each group. Currently support is per-million ("PM") and "max".}

\item{max_val}{A maximum value to use for scaling the y-values in each track.}

\item{window_size}{Window bin size to use for region down-sampling. If NULL, will not downsample.}

\item{window_mode}{If using window_size, what value to use for each window. Options are "max","mean", and "median".}

\item{target_color}{The color of the background rectangle to highlight the region in ucsc_loc.}

\item{highlight_color}{The color the background rectangle to highlight the region in highlight_loc.}
}
\value{
A ggplot2 plot object
}
\description{
Build a multi-track pileup plot
}
