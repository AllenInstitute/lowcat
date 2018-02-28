merge_bam_files <- function(bam_files,
                            out_file,
                            make_indexes = FALSE,
                            sort_out_file = TRUE,
                            index_out_file = TRUE) {

  library(rbamtools)

  l_exist <- file.exists(bam_files)
  n_exist <- sum(l_exist)
  n_missing <- sum(!l_exist)

  l_indexed <- file.exists(paste0(bam_files,".bai"))
  n_indexed <- sum(l_indexed)
  n_no_index <- sum(!l_indexed)

  # Check for files and indexes
  if(n_exist < length(bam_files)) {
    stop(paste0("Can't find", n_missing ," bam files: ", paste(bam_files[!l_exist], collapse = ", ")))
  }
  if(n_indexed < length(bam_files)) {
    if(make_indexes) {
      print("Missing ", n_no_index, " indexes. These will be generated.")
    } else {
      stop(paste0("Can't find", n_no_index ," bam indexes (bam.bai): ", paste(paste0(bam_files[!l_indexed],".bai"), collapse = ", "), ". Use make_indexes = TRUE to generate them on the fly."))
    }
  }

  # Get header from first file
  bam1 <- bamReader(bam_files[1])
  header <- getHeader(bam1)
  bamClose(bam1)

  # Initialize the new file
  out <- bamWriter(header, out_file)

  for(i in 1:length(bam_files)) {
    # Open each bam
    bam <- bamReader(bam_files[i])

    # If an index doesn't exist, make one
    index_file <- paste0(bam_files[i],".bai")
    if(!l_indexed[i]) {
      createIndex(bam, index_file)
    }

    # Load Index
    loadIndex(bam, index_file)

    # Get all chromosomes
    refs <- getRefData(bam)
    for(j in 1:nrow(refs)) {
      # Get coords for each chromosome
      coords <- getRefCoords(bam, refs$SN[j])
      # Convert to bamRange
      range <- bamRange(bam, coords)
      # Write each chr to out.
      bamSave(out, range)
    }

    bamClose(bam)
  }

  bamClose(out)

  if(sort_out_file | index_out_file) {
    print(paste0("Sorting ",out_file,". (Required for indexing)"))
    out <- bamReader(out_file)
    bamSort(out, prefix = sub(".bam$","",out_file))

    if(index_out_file) {
      createIndex(out, paste0(out_file,".bai"))
    }

    bamClose(out)
  }

}
