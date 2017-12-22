makePeakList <- function(peakDir,peakNames) {

  peakFiles <- list.files(peakDir)

  allPeaks <- list()

  for(i in 1:length(peakFiles)) {

    info <- file.info(paste0(peakDir,peakFiles[i]))
    no_lines <- info$size == 0

    if(no_lines) {
      bedPeaks <- data.frame(chr="chr",start=0,end=0,stringsAsFactors=F)
    } else {
      bedPeaks <- read.table(paste0(peakDir,peakFiles[i]),header=F)
    }
    GRPeaks <- GRanges(     seqnames = bedPeaks[,1],
                            IRanges(start = bedPeaks[,2],
                                    end = bedPeaks[,3]))
    allPeaks <- c(allPeaks,GRPeaks)
    names(allPeaks)[i] <- peakNames[i]

  }

  return(allPeaks)
}


readCounts <- function(bamFiles) {
  counts <- numeric()

  for(b in 1:length(bamFiles)) {

    bam <- import(bamFiles[b],format="bam")
    counts <- c(counts,length(bam))

  }
  return(counts)
}

calculateCounts <- function(peakList,bamFiles) {

  counts <- matrix(nrow=length(peakList),ncol=length(bamFiles))

  for(b in 1:length(bamFiles)) {
    print(paste0("Running ",b," of ",length(bamFiles)))

    bam <- import(bamFiles[b],format="bam")

    for(p in 1:length(peakList)) {
      counts[p,b] <- sum(countOverlaps(bam,peakList[[p]]))

    }
  }

  counts <- as.data.frame(counts)
  names(counts) <- bamFiles
  row.names(counts) <- names(peakList)

  return(counts)

}

calculateBHSUM <- function(peakSet,bamFiles,downsample=NULL,bigmode=F,verbose=F) {

  if(bigmode) {
    peakSet_size <- length(peakSet)
    chunk_size <- 100000
    n_chunks <- floor(peakSet_size / chunk_size)

    bhsum <- matrix(nrow=0,ncol=length(bamFiles))

    for(i in 1:n_chunks) {
      peakSet_chunk <- peakSet[((i - 1) * chunk_size + 1):(i * chunk_size),]

      chunk_bhsum <- matrix(nrow=chunk_size,ncol=length(bamFiles))

      for(b in 1:length(bamFiles)) {

        bam <- import(bamFiles[b],format="bam")

        if(!is.null(downsample)) {
          bam <- sample(bam,downsample)
        }

        bamCounts <- numeric()

        bamCounts <- as.numeric(c(bamCounts,countOverlaps(peakSet_chunk,bam)) > 0)

        chunk_bhsum[,b] <- bamCounts

      }

      chunk_bhsum <- chunk_bhsum[rowSums(chunk_bhsum > 0), ]

      bhsum <- rbind(bhsum,chunk_bhsum)
    }

    # Remaining sites after chunking
    peakSet_chunk <- peakSet[(n_chunks * chunk_size + 1):peakSet_size,]

    chunk_bhsum <- matrix(nrow=chunk_size,ncol=length(bamFiles))

    for(b in 1:length(bamFiles)) {

      bam <- import(bamFiles[b],format="bam")

      if(!is.null(downsample)) {
        ds_n <- downsample[1]
        ds_type <- downsample[2]

        if(ds_type == "reads") {

          bam <- sample(bam,ds_n)

        } else if(ds_type == "fragments") {

          bam_p <- bam[strand(bam) == "+",]
          bam_n <- bam[strand(bam) == "-",]
          bam_f <- GRanges(seqname = seqnames(bam_p),IRanges(start = start(bam_p),end = end(bam_p)))
          bam <- sample(bam,ds_n)

        }
      }

      bamCounts <- numeric()

      bamCounts <- as.numeric(c(bamCounts,countOverlaps(peakSet_chunk,bam)) > 0)

      chunk_bhsum[,b] <- bamCounts

    }

    chunk_bhsum <- chunk_bhsum[rowSums(chunk_bhsum > 0), ]

    bhsum <- rbind(bhsum,chunk_bhsum)

  }  else {

    bhsum <- matrix(nrow=length(peakSet),ncol=length(bamFiles))

    for(b in 1:length(bamFiles)) {

      if(verbose) {
        print(paste("Analyzing",bamFiles[b]))
      }

      bam <- import(bamFiles[b],format="bam")

      if(!is.null(downsample)) {
        bam <- sample(bam,downsample)
      }

      bamCounts <- numeric()

      bamCounts <- as.numeric(c(bamCounts,countOverlaps(peakSet,bam)) > 0)

      bhsum[,b] <- bamCounts

    }

  }
  names(bhsum) <- bamFiles

  return(bhsum)

}


calculateRPKMs <- function(peakList,bamFiles) {

  rpkm <- matrix(nrow=length(peakList),ncol=length(bamFiles))

  for(b in 1:length(bamFiles)) {

    bam <- import(bamFiles[b],format="bam")
    m <- length(bam)/1e6

    for(p in 1:length(peakList)) {

      k <- sum(width(peakList[[p]]))/1e3

      rpkm[p,b] <- sum(countOverlaps(bam,peakList[[p]]))/k/m

    }
  }

  rpkm <- as.data.frame(rpkm)
  names(rpkm) <- bamFiles
  row.names(rpkm) <- names(peakList)

  return(rpkm)
}


subsampleRPKMs <- function(peakList,bamFiles,subcount=100000) {

  rpkm <- matrix(nrow=length(peakList),ncol=length(bamFiles))

  for(b in 1:length(bamFiles)) {

    bam <- import(bamFiles[b],format="bam")

    bam <- sample(bam,subcount)

    m <- length(bam)/1e6

    for(p in 1:length(peakList)) {

      k <- sum(width(peakList[[p]]))/1e3

      rpkm[p,b] <- sum(countOverlaps(bam,peakList[[p]]))/k/m

    }
  }

  rpkm <- as.data.frame(rpkm)
  names(rpkm) <- bamFiles
  row.names(rpkm) <- names(peakList)

  return(rpkm)
}

calculateBinomTests <- function(peakList,bamFiles,bamFilter=".+",correction=T,min_count=0) {

  mm10Size <- 2652783500

  binomPValues <- matrix(nrow=length(peakList),ncol=length(bamFiles))

  for(b in 1:length(bamFiles)) {

    bam <- import(bamFiles[b],format="bam")

    if(length(bam) >= min_count) {

      for(i in 1:length(peakList)) {

        p <- sum(width(peakList[[i]]))/mm10Size
        n <- length(bam)
        x <- sum(countOverlaps(bam,peakList[[i]]))

        binomPValues[i,b] <- binom.test(x,n,p)$p.value

      }
    } else {
      binomPValues[,b] <- rep(1,length(peakList))
    }
  }

  if(correction) {
    binomPValues <- binomPValues * ( ncol(binomPValues) * nrow(binomPValues) )
    binomPValues[binomPValues > 1] <- 1
  }

  binomPValues <- as.data.frame(binomPValues)
  names(binomPValues) <- bamFiles
  row.names(binomPValues) <- names(peakList)

  return(binomPValues)

}


plotRPKMs <- function(rpkm,lowCut=0,highCut=5) {


  rpkmPlotData <- data.frame(val=unlist(data.frame(rpkm)),
                             x=rep(1:ncol(rpkm),each=nrow(rpkm)),
                             yname=rep(rownames(rpkm),ncol(rpkm)),
                             color="white")

  yorder <- data.frame(yname=rownames(rpkm),y=1:length(rownames(rpkm))) %>% arrange(y)

  rpkmPlotData <- rpkmPlotData %>% left_join(yorder)

  colorSet <- colorRampPalette(c("darkblue","orangered"))(1001)
  lowColor <- "white"
  highColor <- "red"

  rpkmPlotData$color[rpkmPlotData$val < lowCut] <- lowColor
  rpkmPlotData$color[rpkmPlotData$val > highCut] <- highColor
  intermediates <- rpkmPlotData$val >= lowCut & rpkmPlotData$val <= highCut
  rpkmPlotData$color[intermediates] <- colorSet[round((rpkmPlotData$val[intermediates]-lowCut)/(highCut-lowCut)*1000,0)+1]


  xlab <- names(rpkm)

  ylab <- yorder$yname

  rpkmPlot <- ggplot(rpkmPlotData) +
    geom_tile(aes(x=x,y=y,fill=color)) +
    scale_fill_identity() +
    scale_y_continuous(labels=ylab,breaks=1:nrow(rpkm),expand=c(0,0),name="peakSet") +
    scale_x_continuous(labels=xlab,breaks=1:ncol(rpkm),expand=c(0,0),name="bamFile") +
    theme_classic() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.3))

  return(rpkmPlot)
}

plotBinarized <- function(rpkm) {

  binarized <- matrix(nrow=nrow(rpkm),ncol=ncol(rpkm))

  for(b in 1:ncol(rpkm)) {

    binarized[,b] <- as.numeric(rpkm[,b] == max(rpkm[,b]))

  }

  binarizedPlotData <- data.frame(val=unlist(data.frame(binarized)),
                                  x=rep(1:ncol(rpkm),each=nrow(rpkm)),
                                  yname=rep(rownames(rpkm),ncol(rpkm)))

  yorder <- data.frame(yname=rownames(rpkm),y=1:length(rownames(rpkm))) %>% arrange(y)
  binarizedPlotData <- binarizedPlotData %>% left_join(yorder)

  xlab <- names(rpkm)
  ylab <- yorder$yname

  binarizedPlot <- ggplot(binarizedPlotData) +
    geom_tile(aes(x=x,y=y,fill=val),color="black") +
    scale_fill_gradient2(low="white",mid="white",high="orangered",midpoint=max(binarizedPlotData$val)/2) +
    scale_y_continuous(labels=ylab,breaks=1:nrow(rpkm),expand=c(0,0),name="peakSet") +
    scale_x_continuous(labels=xlab,breaks=1:ncol(rpkm),expand=c(0,0),name="bamFile") +
    theme_classic() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.3))

  return(binarizedPlot)
}

plotBinomTests <- function(rpkm,counts=F,lowCut=3,highCut=50) {

  rpkmPlotData <- data.frame(val=-log10(unlist(data.frame(rpkm))),
                             x=rep(1:ncol(rpkm),each=nrow(rpkm)),
                             yname=rep(rownames(rpkm),ncol(rpkm)),
                             color="white")

  yorder <- data.frame(yname=rownames(rpkm),y=1:length(rownames(rpkm))) %>% arrange(y)

  rpkmPlotData <- rpkmPlotData %>% left_join(yorder)

  colorSet <- colorRampPalette(c("darkblue","orangered"))(1001)
  lowColor <- "white"
  highColor <- "red"

  rpkmPlotData$color[rpkmPlotData$val < lowCut] <- lowColor
  rpkmPlotData$color[rpkmPlotData$val > highCut] <- highColor
  intermediates <- rpkmPlotData$val >= lowCut & rpkmPlotData$val <= highCut
  rpkmPlotData$color[intermediates] <- colorSet[round((rpkmPlotData$val[intermediates]-lowCut)/(highCut-lowCut)*1000,0)+1]

  if(counts) {
    xlab <- counts
  } else {
    xlab <- names(rpkm)
  }
  ylab <- yorder$yname

  rpkmPlot <- ggplot(rpkmPlotData) +
    geom_tile(aes(x=x,y=y,fill=color),color="black") +
    scale_fill_identity() +
    scale_y_continuous(labels=ylab,breaks=1:nrow(rpkm),expand=c(0,0),name="peakSet") +
    scale_x_continuous(labels=xlab,breaks=1:ncol(rpkm),expand=c(0,0),name="bamFile") +
    theme_classic() +
    coord_flip() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.3))

  return(rpkmPlot)
}


calculateHyperTests <- function(peakList,bamDir,bamFilter=".+",correction=T) {

  n_all_peaks <- sum(unlist(lapply(peakList,length)))

  bamFiles <- list.files(bamDir,bamFilter)

  hyperPValues <- matrix(nrow=length(peakList),ncol=length(bamFiles))

  for(b in 1:length(bamFiles)) {

    bam <- import(paste0(bamDir,bamFiles[b]),format="bam")

    for(i in 1:length(peakList)) {

      q <- sum(countOverlaps(bam,peakList[[i]]))
      m <- length(peakList[[i]])
      n <- n_all_peaks - m
      k <- length(bam)

      hyperPValues[i,b] <- 1 - phyper(q, m, n, k)

    }
  }

  if(correction) {
    hyperPValues <- hyperPValues * ( ncol(hyperPValues) * nrow(hyperPValues) )
    hyperPValues[hyperPValues > 1] <- 1
  }

  hyperPValues <- as.data.frame(hyperPValues)
  names(hyperPValues) <- bamFiles
  row.names(hyperPValues) <- names(peakList)

  return(hyperPValues)

}

plotHyperTests <- function(rpkm,counts=F) {

  rpkmPlotData <- data.frame(val=-log10(unlist(data.frame(rpkm))),
                             x=rep(1:ncol(rpkm),each=nrow(rpkm)),
                             yname=rep(rownames(rpkm),ncol(rpkm)),
                             color="white")

  yorder <- data.frame(yname=rownames(rpkm),y=1:length(rownames(rpkm))) %>% arrange(y)

  rpkmPlotData <- rpkmPlotData %>% left_join(yorder)

  colorSet <- colorRampPalette(c("darkblue","orangered"))(1001)
  lowCut <- 1
  lowColor <- "white"
  highCut <- 25
  highColor <- "red"

  rpkmPlotData$color[rpkmPlotData$val < lowCut] <- lowColor
  rpkmPlotData$color[rpkmPlotData$val > highCut] <- highColor
  intermediates <- rpkmPlotData$val >= lowCut & rpkmPlotData$val <= highCut
  rpkmPlotData$color[intermediates] <- colorSet[round((rpkmPlotData$val[intermediates]-lowCut)/(highCut-lowCut)*1000,0)+1]

  if(counts) {
    xlab <- counts
  } else {
    xlab <- names(rpkm)
  }
  ylab <- yorder$yname

  rpkmPlot <- ggplot(rpkmPlotData) +
    geom_tile(aes(x=x,y=y,fill=color),color="black") +
    scale_fill_identity() +
    scale_y_continuous(labels=ylab,breaks=1:nrow(rpkm),expand=c(0,0),name="peakSet") +
    scale_x_continuous(labels=xlab,breaks=1:ncol(rpkm),expand=c(0,0),name="bamFile") +
    theme_classic() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.3))

  return(rpkmPlot)
}


