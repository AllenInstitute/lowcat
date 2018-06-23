#' Read TSS regions from the UCSC RefGenes table
#'
#' This function requires an internet connection.
#'
#' @param symbols A character object with gene symbols to retain. If NULL, will return all genes. Default is NULL.
#' @param expand Distance from TSS to expand. If a single value is provided, will use this values in both directions.
#' If two values are provided, they will be used to expand regions in the 5' and 3' direction relative to transcription direction. Default is 5000.
#' @param genome The genome to use. Default is "mm10".
#'
#' @return a data.frame with columns matching BED format (chr, start, end, name, score, and strand).
get_tss_regions <- function(symbols = NULL,
                            expand = 5000,
                            genome = "mm10") {
  library(dplyr)
  library(rtracklayer)

  session <- browserSession("UCSC")
  genome(session) <- genome
  refgene <- getTable(ucscTableQuery(session,table="refGene"))

  if(length(expand) == 1) {
    expand5 <- expand
    expand3 <- expand
  } else if(length(expand) == 2) {
    expand5 <- expand[1]
    expand3 <- expand[2]
  }

  if(!is.null(symbols)) {

    refgene <- refgene %>%
      filter(name2 %in% symbols)

  }

  tss_windows <- refgene %>%
    mutate(tss = ifelse(strand == "+", txStart, txEnd)) %>%
    mutate(tss_up = ifelse(strand == "+", tss - expand5, tss - expand3),
           tss_dn = ifelse(strand == "+", tss + expand3, tss + expand5)) %>%
    select(chrom,tss_up,tss_dn,name2,score,strand)

  names(tss_windows) <- c("chr","start","end","name","score","strand")

  return(tss_windows)

}

#' Read gene body locations from the UCSC RefGenes table
#'
#' @param symbols A character object with gene symbols to retain. If NULL, will return all genes. Default is NULL.
#' @param expand Distance from the ends of the gene bodies to expand. If a single value is provided, will use this values in both directions.
#' If two values are provided, they will be used to expand regions in the 5' and 3' direction relative to transcription direction. Default is 0.
#' @param genome The genome to use. Default is "mm10".
#'
#' @return a data.frame with columns matching BED format (chr, start, end, name, score, and strand).
get_gene_bodies <- function(symbols = NULL,
                            expand = 0,
                            genome = "mm10") {
  library(dplyr)
  library(rtracklayer)

  session <- browserSession("UCSC")
  genome(session) <- genome
  refgene <- getTable(ucscTableQuery(session,table="refGene"))

  if(length(expand) == 1) {
    expand5 <- expand
    expand3 <- expand
  } else if(length(expand) == 2) {
    expand5 <- expand[1]
    expand3 <- expand[2]
  }

  if(!is.null(symbols)) {

    refgene <- refgene %>%
      filter(name2 %in% symbols)

  }

  gene_bodies <- refgene %>%
    mutate(tss_up = ifelse(strand == "+", txStart - expand5, txStart - expand3),
           end_dn = ifelse(strand == "+", txEnd + expand3, txEnd + expand5)) %>%
    select(chrom,tss_up,end_dn,name2,score,strand) %>%
    unique()

  names(gene_bodies) <- c("chr","start","end","name","score","strand")

  return(gene_bodies)

}

#' Generate regions similar to GREAT based on the UCSC RefGenes table
#'
#' GREAT is the [Genomic Regions Enrichment of Annotations Tool](http://great.stanford.edu).
#'
#' GREAT regions have a minimum size around the TSS. In the real deal, this is an asymmetrical region (more upstream than downstream), but
#' in this implementation, the minimum region specified by minexpand is symmetrical. Regions then expand up to the maximum distance or until
#' they hit the minimum region of an adjacent gene, whichever is smaller.
#'
#' @param symbols A character object with gene symbols to retain. If NULL, will return all genes. Default is NULL.
#' @param minexpand Minimum distance from the TSS to expand. This region will be retained even if it overlaps an adjacent gene. Default is 5000.
#' @param maxexpand Maximum distance from the TSS to expand. This region will be bounded by adjacent minimum regions of the nearest gene.
#' @param genome The genome to use. Default is "mm10".
#'
#' @return a data.frame with columns matching BED format (chr, start, end, name, score, and strand).
get_great_regions <- function(symbols = NULL,
                              minexpand = 5000,
                              maxexpand = 1e6,
                              genome = "mm10") {
  library(dplyr)
  library(rtracklayer)

  session <- browserSession("UCSC")
  genome(session) <- genome
  refgene <- getTable(ucscTableQuery(session,table="refGene"))
  chromsizes <- read.table(paste0("http://hgdownload.cse.ucsc.edu/goldenPath/",genome,"/bigZips/",genome,".chrom.sizes"))
  names(chromsizes) <- c("chrom","max_chr_size")

  if(!is.null(symbols)) {

    refgene <- refgene %>%
      filter(name2 %in% symbols)

  }

  tss_windows <- refgene %>%
    left_join(chromsizes) %>%
    mutate(tss = ifelse(strand == "+", txStart, txEnd)) %>%
    arrange(chrom,tss) %>%
    mutate(min_up = tss - minexpand,
           min_dn = tss + minexpand) %>%
    group_by(chrom) %>%
    mutate(max_up = ifelse(lag(min_dn, default = 0) < tss - maxexpand, # If the nearest minimal region is > maxexpand away
                           tss - maxexpand, # Set the window boundary to tss - maxexpand
                           ifelse(lag(min_dn, default = 0) > min_up, # If the nearest minimal region is < minexpand away
                                  min_up, # Set the window boundary to tss - minexpand
                                  lag(min_dn, default = 0) + 1)) # Otherwise, set the boundary to be just outside the minexpand of the nearby gene.
    ) %>%
    mutate(max_dn = ifelse(lead(min_up, default = 0) > tss + maxexpand, # If the nearest minimal region is > maxexpand away
                           tss + maxexpand, # Set the window boundary to tss + maxexpand
                           ifelse(lead(min_up, default = 0) < min_dn, # If the nearest minimal region is < minexpand away
                                  min_dn, # Set the window boundary to tss + minexpand
                                  lead(min_up, default = 0) - 1))) %>%
    mutate(max_up = ifelse(max_up < 0, 0, max_up),
           max_dn = ifelse(max_dn > max_chr_size, max_chr_size, max_dn)) %>%
    ungroup() %>%
    select(chrom,max_up,max_dn,name2,score,strand)

  names(tss_windows) <- c("chr","start","end","name","score","strand")

  return(tss_windows)

}


#' Convert data.frames in BED-like format to GRanges objects
#'
#' @param bed A data.frame with columns chr, start, end, name, strand, and score.
#'
#' @return A GenomicRanges object
#'
bed_to_GRanges <- function(bed) {
  library(rtracklayer)

  gr <- GRanges(seqnames=bed$chr,
                IRanges(start=bed$start,
                        end=bed$end),
                strand=bed$strand,
                mcols=bed[,c("name","score")])

  return(gr)
}

#' Convert UCSC-style genomic locations to a GenomicRegions object.
#'
#' @param ucsc_loc A UCSC-like genomic location. Example: "chr1:152,548,974-152,550,854"
#'
#' @return a GenomicRanges object for the UCSC location
#'
ucsc_loc_to_gr <- function(ucsc_loc) {
  chr <- sub(":.+","",ucsc_loc)
  start <- sub(".+:","",sub("-.+","",ucsc_loc))
  start <- as.numeric(gsub(",","",start))
  end <- sub(".+-","",ucsc_loc)
  end <- as.numeric(gsub(",","",end))
  GRanges(seqnames = chr,
          IRanges(start,end),
          strand = "+")
}

