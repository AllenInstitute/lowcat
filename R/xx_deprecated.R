# ### These are historical functions
# # Still here for reference, but not exported.
#
# makePeakList <- function(peakDir,peakNames) {
#
#   peakFiles <- list.files(peakDir)
#
#   allPeaks <- list()
#
#   for(i in 1:length(peakFiles)) {
#
#     info <- file.info(paste0(peakDir,peakFiles[i]))
#     no_lines <- info$size == 0
#
#     if(no_lines) {
#       bedPeaks <- data.frame(chr="chr",start=0,end=0,stringsAsFactors=F)
#     } else {
#       bedPeaks <- read.table(paste0(peakDir,peakFiles[i]),header=F)
#     }
#     GRPeaks <- GRanges(     seqnames = bedPeaks[,1],
#                             IRanges(start = bedPeaks[,2],
#                                     end = bedPeaks[,3]))
#     allPeaks <- c(allPeaks,GRPeaks)
#     names(allPeaks)[i] <- peakNames[i]
#
#   }
#
#   return(allPeaks)
# }
#
#
# readCounts <- function(bamFiles) {
#   counts <- numeric()
#
#   for(b in 1:length(bamFiles)) {
#
#     bam <- import(bamFiles[b],format="bam")
#     counts <- c(counts,length(bam))
#
#   }
#   return(counts)
# }
#
# calculateCounts <- function(peakList,bamFiles) {
#
#   counts <- matrix(nrow=length(peakList),ncol=length(bamFiles))
#
#   for(b in 1:length(bamFiles)) {
#     print(paste0("Running ",b," of ",length(bamFiles)))
#
#     bam <- import(bamFiles[b],format="bam")
#
#     for(p in 1:length(peakList)) {
#       counts[p,b] <- sum(countOverlaps(bam,peakList[[p]]))
#
#     }
#   }
#
#   counts <- as.data.frame(counts)
#   names(counts) <- bamFiles
#   row.names(counts) <- names(peakList)
#
#   return(counts)
#
# }
#
# calculateBHSUM <- function(peakSet,bamFiles,downsample=NULL,bigmode=F,verbose=F) {
#
#   if(bigmode) {
#     peakSet_size <- length(peakSet)
#     chunk_size <- 100000
#     n_chunks <- floor(peakSet_size / chunk_size)
#
#     bhsum <- matrix(nrow=0,ncol=length(bamFiles))
#
#     for(i in 1:n_chunks) {
#       peakSet_chunk <- peakSet[((i - 1) * chunk_size + 1):(i * chunk_size),]
#
#       chunk_bhsum <- matrix(nrow=chunk_size,ncol=length(bamFiles))
#
#       for(b in 1:length(bamFiles)) {
#
#         bam <- import(bamFiles[b],format="bam")
#
#         if(!is.null(downsample)) {
#           bam <- sample(bam,downsample)
#         }
#
#         bamCounts <- numeric()
#
#         bamCounts <- as.numeric(c(bamCounts,countOverlaps(peakSet_chunk,bam)) > 0)
#
#         chunk_bhsum[,b] <- bamCounts
#
#       }
#
#       chunk_bhsum <- chunk_bhsum[rowSums(chunk_bhsum > 0), ]
#
#       bhsum <- rbind(bhsum,chunk_bhsum)
#     }
#
#     # Remaining sites after chunking
#     peakSet_chunk <- peakSet[(n_chunks * chunk_size + 1):peakSet_size,]
#
#     chunk_bhsum <- matrix(nrow=chunk_size,ncol=length(bamFiles))
#
#     for(b in 1:length(bamFiles)) {
#
#       bam <- import(bamFiles[b],format="bam")
#
#       if(!is.null(downsample)) {
#         ds_n <- downsample[1]
#         ds_type <- downsample[2]
#
#         if(ds_type == "reads") {
#
#           bam <- sample(bam,ds_n)
#
#         } else if(ds_type == "fragments") {
#
#           bam_p <- bam[strand(bam) == "+",]
#           bam_n <- bam[strand(bam) == "-",]
#           bam_f <- GRanges(seqname = seqnames(bam_p),IRanges(start = start(bam_p),end = end(bam_p)))
#           bam <- sample(bam,ds_n)
#
#         }
#       }
#
#       bamCounts <- numeric()
#
#       bamCounts <- as.numeric(c(bamCounts,countOverlaps(peakSet_chunk,bam)) > 0)
#
#       chunk_bhsum[,b] <- bamCounts
#
#     }
#
#     chunk_bhsum <- chunk_bhsum[rowSums(chunk_bhsum > 0), ]
#
#     bhsum <- rbind(bhsum,chunk_bhsum)
#
#   }  else {
#
#     bhsum <- matrix(nrow=length(peakSet),ncol=length(bamFiles))
#
#     for(b in 1:length(bamFiles)) {
#
#       if(verbose) {
#         print(paste("Analyzing",bamFiles[b]))
#       }
#
#       bam <- import(bamFiles[b],format="bam")
#
#       if(!is.null(downsample)) {
#         bam <- sample(bam,downsample)
#       }
#
#       bamCounts <- numeric()
#
#       bamCounts <- as.numeric(c(bamCounts,countOverlaps(peakSet,bam)) > 0)
#
#       bhsum[,b] <- bamCounts
#
#     }
#
#   }
#   names(bhsum) <- bamFiles
#
#   return(bhsum)
#
# }
#
#
# calculateRPKMs <- function(peakList,bamFiles) {
#
#   rpkm <- matrix(nrow=length(peakList),ncol=length(bamFiles))
#
#   for(b in 1:length(bamFiles)) {
#
#     bam <- import(bamFiles[b],format="bam")
#     m <- length(bam)/1e6
#
#     for(p in 1:length(peakList)) {
#
#       k <- sum(width(peakList[[p]]))/1e3
#
#       rpkm[p,b] <- sum(countOverlaps(bam,peakList[[p]]))/k/m
#
#     }
#   }
#
#   rpkm <- as.data.frame(rpkm)
#   names(rpkm) <- bamFiles
#   row.names(rpkm) <- names(peakList)
#
#   return(rpkm)
# }
#
#
# subsampleRPKMs <- function(peakList,bamFiles,subcount=100000) {
#
#   rpkm <- matrix(nrow=length(peakList),ncol=length(bamFiles))
#
#   for(b in 1:length(bamFiles)) {
#
#     bam <- import(bamFiles[b],format="bam")
#
#     bam <- sample(bam,subcount)
#
#     m <- length(bam)/1e6
#
#     for(p in 1:length(peakList)) {
#
#       k <- sum(width(peakList[[p]]))/1e3
#
#       rpkm[p,b] <- sum(countOverlaps(bam,peakList[[p]]))/k/m
#
#     }
#   }
#
#   rpkm <- as.data.frame(rpkm)
#   names(rpkm) <- bamFiles
#   row.names(rpkm) <- names(peakList)
#
#   return(rpkm)
# }
#
# calculateBinomTests <- function(peakList,bamFiles,bamFilter=".+",correction=T,min_count=0) {
#
#   mm10Size <- 2652783500
#
#   binomPValues <- matrix(nrow=length(peakList),ncol=length(bamFiles))
#
#   for(b in 1:length(bamFiles)) {
#
#     bam <- import(bamFiles[b],format="bam")
#
#     if(length(bam) >= min_count) {
#
#       for(i in 1:length(peakList)) {
#
#         p <- sum(width(peakList[[i]]))/mm10Size
#         n <- length(bam)
#         x <- sum(countOverlaps(bam,peakList[[i]]))
#
#         binomPValues[i,b] <- binom.test(x,n,p)$p.value
#
#       }
#     } else {
#       binomPValues[,b] <- rep(1,length(peakList))
#     }
#   }
#
#   if(correction) {
#     binomPValues <- binomPValues * ( ncol(binomPValues) * nrow(binomPValues) )
#     binomPValues[binomPValues > 1] <- 1
#   }
#
#   binomPValues <- as.data.frame(binomPValues)
#   names(binomPValues) <- bamFiles
#   row.names(binomPValues) <- names(peakList)
#
#   return(binomPValues)
#
# }
#
#
# plotRPKMs <- function(rpkm,lowCut=0,highCut=5) {
#
#
#   rpkmPlotData <- data.frame(val=unlist(data.frame(rpkm)),
#                              x=rep(1:ncol(rpkm),each=nrow(rpkm)),
#                              yname=rep(rownames(rpkm),ncol(rpkm)),
#                              color="white")
#
#   yorder <- data.frame(yname=rownames(rpkm),y=1:length(rownames(rpkm))) %>% arrange(y)
#
#   rpkmPlotData <- rpkmPlotData %>% left_join(yorder)
#
#   colorSet <- colorRampPalette(c("darkblue","orangered"))(1001)
#   lowColor <- "white"
#   highColor <- "red"
#
#   rpkmPlotData$color[rpkmPlotData$val < lowCut] <- lowColor
#   rpkmPlotData$color[rpkmPlotData$val > highCut] <- highColor
#   intermediates <- rpkmPlotData$val >= lowCut & rpkmPlotData$val <= highCut
#   rpkmPlotData$color[intermediates] <- colorSet[round((rpkmPlotData$val[intermediates]-lowCut)/(highCut-lowCut)*1000,0)+1]
#
#
#   xlab <- names(rpkm)
#
#   ylab <- yorder$yname
#
#   rpkmPlot <- ggplot(rpkmPlotData) +
#     geom_tile(aes(x=x,y=y,fill=color)) +
#     scale_fill_identity() +
#     scale_y_continuous(labels=ylab,breaks=1:nrow(rpkm),expand=c(0,0),name="peakSet") +
#     scale_x_continuous(labels=xlab,breaks=1:ncol(rpkm),expand=c(0,0),name="bamFile") +
#     theme_classic() +
#     theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.3))
#
#   return(rpkmPlot)
# }
#
# plotBinarized <- function(rpkm) {
#
#   binarized <- matrix(nrow=nrow(rpkm),ncol=ncol(rpkm))
#
#   for(b in 1:ncol(rpkm)) {
#
#     binarized[,b] <- as.numeric(rpkm[,b] == max(rpkm[,b]))
#
#   }
#
#   binarizedPlotData <- data.frame(val=unlist(data.frame(binarized)),
#                                   x=rep(1:ncol(rpkm),each=nrow(rpkm)),
#                                   yname=rep(rownames(rpkm),ncol(rpkm)))
#
#   yorder <- data.frame(yname=rownames(rpkm),y=1:length(rownames(rpkm))) %>% arrange(y)
#   binarizedPlotData <- binarizedPlotData %>% left_join(yorder)
#
#   xlab <- names(rpkm)
#   ylab <- yorder$yname
#
#   binarizedPlot <- ggplot(binarizedPlotData) +
#     geom_tile(aes(x=x,y=y,fill=val),color="black") +
#     scale_fill_gradient2(low="white",mid="white",high="orangered",midpoint=max(binarizedPlotData$val)/2) +
#     scale_y_continuous(labels=ylab,breaks=1:nrow(rpkm),expand=c(0,0),name="peakSet") +
#     scale_x_continuous(labels=xlab,breaks=1:ncol(rpkm),expand=c(0,0),name="bamFile") +
#     theme_classic() +
#     theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.3))
#
#   return(binarizedPlot)
# }
#
# plotBinomTests <- function(rpkm,counts=F,lowCut=3,highCut=50) {
#
#   rpkmPlotData <- data.frame(val=-log10(unlist(data.frame(rpkm))),
#                              x=rep(1:ncol(rpkm),each=nrow(rpkm)),
#                              yname=rep(rownames(rpkm),ncol(rpkm)),
#                              color="white")
#
#   yorder <- data.frame(yname=rownames(rpkm),y=1:length(rownames(rpkm))) %>% arrange(y)
#
#   rpkmPlotData <- rpkmPlotData %>% left_join(yorder)
#
#   colorSet <- colorRampPalette(c("darkblue","orangered"))(1001)
#   lowColor <- "white"
#   highColor <- "red"
#
#   rpkmPlotData$color[rpkmPlotData$val < lowCut] <- lowColor
#   rpkmPlotData$color[rpkmPlotData$val > highCut] <- highColor
#   intermediates <- rpkmPlotData$val >= lowCut & rpkmPlotData$val <= highCut
#   rpkmPlotData$color[intermediates] <- colorSet[round((rpkmPlotData$val[intermediates]-lowCut)/(highCut-lowCut)*1000,0)+1]
#
#   if(counts) {
#     xlab <- counts
#   } else {
#     xlab <- names(rpkm)
#   }
#   ylab <- yorder$yname
#
#   rpkmPlot <- ggplot(rpkmPlotData) +
#     geom_tile(aes(x=x,y=y,fill=color),color="black") +
#     scale_fill_identity() +
#     scale_y_continuous(labels=ylab,breaks=1:nrow(rpkm),expand=c(0,0),name="peakSet") +
#     scale_x_continuous(labels=xlab,breaks=1:ncol(rpkm),expand=c(0,0),name="bamFile") +
#     theme_classic() +
#     coord_flip() +
#     theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.3))
#
#   return(rpkmPlot)
# }
#
#
# calculateHyperTests <- function(peakList,bamDir,bamFilter=".+",correction=T) {
#
#   n_all_peaks <- sum(unlist(lapply(peakList,length)))
#
#   bamFiles <- list.files(bamDir,bamFilter)
#
#   hyperPValues <- matrix(nrow=length(peakList),ncol=length(bamFiles))
#
#   for(b in 1:length(bamFiles)) {
#
#     bam <- import(paste0(bamDir,bamFiles[b]),format="bam")
#
#     for(i in 1:length(peakList)) {
#
#       q <- sum(countOverlaps(bam,peakList[[i]]))
#       m <- length(peakList[[i]])
#       n <- n_all_peaks - m
#       k <- length(bam)
#
#       hyperPValues[i,b] <- 1 - phyper(q, m, n, k)
#
#     }
#   }
#
#   if(correction) {
#     hyperPValues <- hyperPValues * ( ncol(hyperPValues) * nrow(hyperPValues) )
#     hyperPValues[hyperPValues > 1] <- 1
#   }
#
#   hyperPValues <- as.data.frame(hyperPValues)
#   names(hyperPValues) <- bamFiles
#   row.names(hyperPValues) <- names(peakList)
#
#   return(hyperPValues)
#
# }
#
# plotHyperTests <- function(rpkm,counts=F) {
#
#   rpkmPlotData <- data.frame(val=-log10(unlist(data.frame(rpkm))),
#                              x=rep(1:ncol(rpkm),each=nrow(rpkm)),
#                              yname=rep(rownames(rpkm),ncol(rpkm)),
#                              color="white")
#
#   yorder <- data.frame(yname=rownames(rpkm),y=1:length(rownames(rpkm))) %>% arrange(y)
#
#   rpkmPlotData <- rpkmPlotData %>% left_join(yorder)
#
#   colorSet <- colorRampPalette(c("darkblue","orangered"))(1001)
#   lowCut <- 1
#   lowColor <- "white"
#   highCut <- 25
#   highColor <- "red"
#
#   rpkmPlotData$color[rpkmPlotData$val < lowCut] <- lowColor
#   rpkmPlotData$color[rpkmPlotData$val > highCut] <- highColor
#   intermediates <- rpkmPlotData$val >= lowCut & rpkmPlotData$val <= highCut
#   rpkmPlotData$color[intermediates] <- colorSet[round((rpkmPlotData$val[intermediates]-lowCut)/(highCut-lowCut)*1000,0)+1]
#
#   if(counts) {
#     xlab <- counts
#   } else {
#     xlab <- names(rpkm)
#   }
#   ylab <- yorder$yname
#
#   rpkmPlot <- ggplot(rpkmPlotData) +
#     geom_tile(aes(x=x,y=y,fill=color),color="black") +
#     scale_fill_identity() +
#     scale_y_continuous(labels=ylab,breaks=1:nrow(rpkm),expand=c(0,0),name="peakSet") +
#     scale_x_continuous(labels=xlab,breaks=1:ncol(rpkm),expand=c(0,0),name="bamFile") +
#     theme_classic() +
#     theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.3))
#
#   return(rpkmPlot)
# }
#
#
# # # Returns only jd
# # This version is slower.
# #
# # fragment_overlap_jaccard_parallel3 <- function(N) {
# #   ol <- findOverlaps(fragment_list[[N]],
# #                       fragment_list[(N+1):length(fragment_list)])
# #
# #   ol <- table(subjectHits(ol))
# #
# #   n1 <- fragment_lengths[N]
# #   n2 <- fragment_lengths[(N+1):length(fragment_list)]
# #
# #   jd <- 1 - ol / (n1 + n2)
# #
# #   matrix(c(jd, rep(N, length(jd))), ncol = 2)
# # }
# #
# # # Compute jaccard similarity values based on GenomicRanges overlaps in parallel.
# # #
# # # @param fragment_list The list object containing GenomicRanges objects.
# # # @param sample_names Sample names. If NULL, will use BAM file names.
# # # @param n_cores The number of cores to use in parallel. Use "auto" to detect and use all cores. Default is 6.
# # # @param cluster_type EXPERIMENTAL: Either "PSOCK" (for Windows) or "FORK" (possible on Linux and Mac). FORK is more memory-efficient.
# # #
# # # @return a list of GenomicRanges objects
# # #
# # run_fragment_overlap_jaccard_parallel3 <- function(fragment_list,
# #                                                    n_cores = 6,
# #                                                    cluster_type = "PSOCK") {
# #
# #   fragment_list <- GRangesList(fragment_list)
# #
# #   #index_pairs <- t(combn(1:length(fragment_list),2))
# #   fragment_lengths <- unlist(lapply(fragment_list, length))
# #
# #   # Set up parallelization
# #   if(n_cores == "auto") {
# #     n_cores <- parallel::detectCores()
# #   }
# #
# #   print(paste("Starting",n_cores,"nodes"))
# #
# #   cl <- parallel::makeCluster(n_cores, cluster_type)
# #
# #   print("Exporting necessary objects to nodes")
# #
# #   parallel::clusterEvalQ(cl, library(GenomicRanges))
# #   parallel::clusterExport(cl, c("fragment_list",
# #                                 "fragment_lengths",
# #                                 "fragment_overlap_jaccard_parallel"),
# #                           # Use the function's local environment for export
# #                           envir = environment())
# #
# #   N <- length(fragment_list) - 1
# #
# #   print(paste("Running",N,"comparisons."))
# #
# #   res <- clusterApplyLB_chunks(N = N,
# #                                n_chunks = 20,
# #                                cl = cl,
# #                                FUN = fragment_overlap_jaccard_parallel3)
# #
# #   stopCluster(cl)
# #
# #   print("Collecting results.")
# #
# #   res <- do.call("rbind", res)
# #   res <- res[,1][order(res[,2])]
# #
# #   res_mat <- matrix(0, length(fragment_list), length(fragment_list))
# #   res_mat[lower.tri(res_mat, diag = FALSE)] <- res
# #   res_mat <- t(res_mat)
# #   res_mat[lower.tri(res_mat, diag = FALSE)] <- res
# #
# #   return(res_mat)
# # }



# # Merge multiple BAM files to a single output BAM file
# #
# # @param bam_files a character vector listing all of the BAM files to merge
# # @param out_file the target BAM file to write to
# # @param make_indexes whether or not to make indexes for the input files. Default is FALSE
# # @param sort_out_file whether or not to sort the output file. default is TRUE, which will force index_out_file to TRUE.
# # @param index_out_file whether or not to index the output file. Default is TRUE.
# # @param keep_unsorted If sort_out_file == TRUE, whether or not to keep the original, unsorted file. Default is FALSE.
# #
# # @export
# merge_bam_files <- function(bam_files,
#                             out_file,
#                             make_indexes = FALSE,
#                             sort_out_file = TRUE,
#                             index_out_file = TRUE,
#                             keep_unsorted = FALSE) {
#
#   library(rbamtools)
#
#   l_exist <- file.exists(bam_files)
#   n_exist <- sum(l_exist)
#   n_missing <- sum(!l_exist)
#
#   l_indexed <- file.exists(paste0(bam_files,".bai"))
#   n_indexed <- sum(l_indexed)
#   n_no_index <- sum(!l_indexed)
#
#   # Check for files and indexes
#   if(n_exist < length(bam_files)) {
#     stop(paste0("Can't find", n_missing ," bam files: ", paste(bam_files[!l_exist], collapse = ", ")))
#   }
#   if(n_indexed < length(bam_files)) {
#     if(make_indexes) {
#       print("Missing ", n_no_index, " indexes. These will be generated.")
#     } else {
#       stop(paste0("Can't find", n_no_index ," bam indexes (bam.bai): ", paste(paste0(bam_files[!l_indexed],".bai"), collapse = ", "), ". Use make_indexes = TRUE to generate them on the fly."))
#     }
#   }
#
#   # Get header from first file
#   bam1 <- bamReader(bam_files[1])
#   header <- getHeader(bam1)
#   bamClose(bam1)
#
#   # Initialize the new file
#   out <- bamWriter(header, out_file)
#
#   for(i in 1:length(bam_files)) {
#     # Open each bam
#     bam <- bamReader(bam_files[i])
#
#     # If an index doesn't exist, make one
#     index_file <- paste0(bam_files[i],".bai")
#     if(!l_indexed[i]) {
#       createIndex(bam, index_file)
#     }
#
#     # Load Index
#     loadIndex(bam, index_file)
#
#     # Get all chromosomes
#     refs <- getRefData(bam)
#     for(j in 1:nrow(refs)) {
#       # Get coords for each chromosome
#       coords <- getRefCoords(bam, refs$SN[j])
#       # Convert to bamRange
#       range <- bamRange(bam, coords)
#       # Write each chr to out.
#       bamSave(out, range)
#     }
#
#     bamClose(bam)
#   }
#
#   bamClose(out)
#
#   if(sort_out_file | index_out_file) {
#     print(paste0("Sorting ",out_file,". (Required for indexing)"))
#     out <- bamReader(out_file)
#     sort_prefix <- sub(".bam$",".srt",out_file)
#     sort_prefix <- sub(".+/","",sort_prefix)
#
#     bamSort(out, prefix = sort_prefix)
#     bamClose(out)
#
#     if(!keep_unsorted) {
#       file.remove(out_file)
#       file.rename(sub(".bam$",".srt.bam",out_file), out_file)
#     }
#
#     if(index_out_file) {
#       if(keep_unsorted) {
#         out <- bamReader(sub(".bam$",".srt.bam",out_file))
#         createIndex(out, paste0(sub(".bam$",".srt.bam",out_file),".bai"))
#         bamClose(out)
#       } else {
#         out <- bamReader(out_file)
#         createIndex(out, paste0(out_file,".bai"))
#         bamClose(out)
#       }
#
#     }
#
#   }
#
# }
