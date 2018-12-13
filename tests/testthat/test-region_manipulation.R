library(lowcat)
options(stringsAsFactors = FALSE)

context("region type conversion")
test_that(
  "GRanges_to_ucsc_loc converts correctly",
  {
    regions_file <- system.file("testdata/rda",
                                "simple_regions.rda",
                                package = "lowcat")

    load(regions_file)

    ucsc_conv1 <- GRanges_to_ucsc_loc(gr1)
    ucsc_conv2 <- GRanges_to_ucsc_loc(gr2)

    expect_identical(ucsc_conv1, ucsc1)
    expect_identical(ucsc_conv2, ucsc2)

  }
)

test_that(
  "ucsc_loc_to_GRanges converts correctly",
  {
    regions_file <- system.file("testdata/rda",
                                "simple_regions.rda",
                                package = "lowcat")

    load(regions_file)

    gr_conv1 <- ucsc_loc_to_GRanges(ucsc1)
    gr_conv2 <- ucsc_loc_to_GRanges(ucsc2)

    gr_comp1 <- gr1
    strand(gr_comp1) <- rep("+",4)
    gr_comp2 <- gr2
    strand(gr_comp2) <- rep("+",4)

    expect_identical(gr_conv1, gr_comp1)
    expect_identical(gr_conv2, gr_comp2)

  }
)

test_that(
  "GRanges_to_bed converts correctly",
  {
    regions_file <- system.file("testdata/rda",
                                "simple_regions.rda",
                                package = "lowcat")

    load(regions_file)

    bed_conv1 <- GRanges_to_bed(gr1)
    bed_conv2 <- GRanges_to_bed(gr2)

    bed_comp1 <- bed1[,c("chr","start","end","strand")]
    bed_comp2 <- bed2[,c("chr","start","end","strand")]

    expect_identical(bed_conv1, bed_comp1)
    expect_identical(bed_conv2, bed_comp2)

  }
)

# test_that(
#   "bed_to_GRanges converts correctly",
#   {
#     regions_file <- system.file("testdata/rda",
#                                 "simple_regions.rda",
#                                 package = "lowcat")
#
#     load(regions_file)
#
#     gr_conv1 <- bed_to_GRanges(bed1)
#     gr_conv2 <- bed_to_GRanges(bed2)
#
#     gr_comp1 <- gr1
#     strand(gr_comp1) <- rep("+",4)
#     gr_comp2 <- gr2
#     strand(gr_comp2) <- rep("+",4)
#
#     expect_identical(gr_conv1, gr1)
#     expect_identical(gr_conv2, gr2)
#
#   }
# )


context("filtering regions")

test_that(
  "filter_fragments() correctly removes regions",
  {
    regions_file <- system.file("testdata/rda",
                                "simple_regions.rda",
                                package = "lowcat")

    load(regions_file)

    filtered <- suppressWarnings(
      filter_fragments(fragment_list = list(gr1),
                                 filter_GR = gr2,
                                 mode = "remove",
                                 ignore_strand = TRUE)
    )

    expect_is(filtered, "list")
    expect_is(filtered[[1]], "GRanges")
    expect_equal(length(filtered[[1]]), 2)

    filtered <- suppressWarnings(
      filter_fragments(fragment_list = list(gr1),
                                 filter_GR = gr2,
                                 mode = "remove",
                                 ignore_strand = FALSE)
    )

    expect_is(filtered, "list")
    expect_is(filtered[[1]], "GRanges")
    expect_equal(length(filtered[[1]]), 3)
  }
)

test_that(
  "filter_fragments() correctly retains regions",
  {
    regions_file <- system.file("testdata/rda",
                                "simple_regions.rda",
                                package = "lowcat")

    load(regions_file)

    filtered <- suppressWarnings(
      filter_fragments(fragment_list = list(gr1),
                                 filter_GR = gr2,
                                 mode = "keep",
                                 ignore_strand = TRUE)
    )

    expect_is(filtered, "list")
    expect_is(filtered[[1]], "GRanges")
    expect_equal(length(filtered[[1]]), 2)

    filtered <- suppressWarnings(
      filter_fragments(fragment_list = list(gr1),
                                 filter_GR = gr2,
                                 mode = "keep",
                                 ignore_strand = FALSE)
    )

    expect_is(filtered, "list")
    expect_is(filtered[[1]], "GRanges")
    expect_equal(length(filtered[[1]]), 1)
  }
)
