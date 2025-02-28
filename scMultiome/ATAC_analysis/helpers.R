#' Mapping & aggregating modalities with genomic coordinates to reference 
#' coordinates.
#'
#' Internal convencience function for mapping different modality scores with 
#' cell type and further labels such as tfs to reference coordinates. The resulting
#' table will have dimension ref coords x byCols (or ref coord x byCols[1] x byCols[2]).
#' ByCols can be for instance cell type labels and/or transcription factor names.
#'
#'@name .genomicRangesMapping
#'@author Emanuel Sonder
#'@param refRanges GRanges object with reference coordinates
#'@param assayTable modality table to be mapped to the reference coordinates. 
#'Needs to containing genomic coordinates (see args: seqNamesCol, startCol, endCol), and byCols.
#'@param byCols will be the columns /depths of the resulting matrix with dimension 
#' ref coords x byCols (or ref coord x byCols[1] x byCols[2]). 
#' ByCols can be for instance cell type labels and/or transcription factor names.
#' @param seqNamesCol name of the column in motif, atac and chIP-seq data.tables containing
#' the sequence information.
#' @param startCol name of the column in motif, atac and chIP-seq data.tables containing
#' the start coordinate.
#' @param endCol name of the column in motif, atac and chIP-seq data.tables containing
#' the end coordinate.
#' @param scoreCol name of the score column (e.g. motif matching scores, atac fragment counts)
#' @param aggregationFun function (e.g. mean, median, sum) used to aggregate
#' if multiple rows of the assayTable overlap a reference coordinate.
#' @param BPPARAM BiocParallel argument either SerialParam() or MulticoreParam(workers=n)
#' @export
.genomicRangesMapping <- function(refRanges, 
                                  assayTable,
                                  byCols=c("tf_uniprot_id",
                                           "cell_type_id"),
                                  seqNamesCol="chr",
                                  startCol="start",
                                  endCol="end",
                                  scoreCol=NULL,
                                  calledInCol=NULL,
                                  aggregationFun=NULL,
                                  minoverlap=1,
                                  BPPARAM=SerialParam()){
  # TODO: - add warning for integer overflows - data.table size
  assayTable <- copy(assayTable)  
  
  # attribute generic names to dimensionalities
  if(length(byCols)==2)
  {
    setnames(assayTable, byCols, c("col_depth", "col_width"))
    byCols <- c("col_depth", "col_width")
    multiTf <- TRUE
  }
  else
  {
    setnames(assayTable, byCols, c("col_width"))
    byCols <- c("col_width")
    multiTf <- FALSE
  }
  
  # get dimensions of tables
  nRefs <- length(refRanges)
  nColsWidth <- length(unique(assayTable$col_width))
  
  # convert to integer for speed-up
  levels <- unique(assayTable$col_width)
  assayTable[,col_width:=as.integer(factor(assayTable$col_width, levels=levels))]
  
  # convert to GRanges for faster overlap finding
  suppressWarnings(assayTable$width <- NULL)
  suppressWarnings(assayTable$strand <- NULL)
  assayRanges <- makeGRangesFromDataFrame(as.data.frame(assayTable),
                                          keep.extra.columns=TRUE,
                                          seqnames.field=seqNamesCol,
                                          start.field=startCol,
                                          end.field=endCol,
                                          ignore.strand=TRUE)
  
  # find overlaps with ref. coordinates
  overlapTable <- as.data.table(findOverlaps(refRanges, assayRanges,
                                             type="any",
                                             minoverlap=minoverlap,
                                             ignore.strand=TRUE))
  rm(refRanges, assayRanges)
  
  # retrieve tf and cell type ids
  overlapTable <- cbind(overlapTable$queryHits, 
                        assayTable[overlapTable$subjectHits, 
                                   c(byCols, scoreCol, calledInCol), 
                                   with=FALSE])
  
  if(multiTf)
  {
    setkey(overlapTable, V1, col_width)
    if(!is.null(calledInCol)) setnames(overlapTable, calledInCol, "calledInCol")
    if(!is.null(scoreCol)) setnames(overlapTable, scoreCol, "scoreCol")
    overlapTable <- split(overlapTable, by=c("col_depth"))
    
    #TODO:  BPPARAM= serialParam or multiCoreParam() as function arguments
    #TODO:  bplapply(fls[1:3], FUN, BPPARAM = MulticoreParam(), param = param)
    # overlap with ref. coordinates
    overlapTable <- BiocParallel::bplapply(overlapTable, function(table){
      
      table <- table[,.(value=aggregationFun(scoreCol)),
                     by=c("V1", "col_width")]
      
      # one would need to construct a second table here for the neg labels
      
      # convert to sparse matrix
      table <- sparseMatrix(i=table$V1, 
                            j=as.integer(table$col_width), # 11.07.2024 as.integer is not needed
                            dims=c(nRefs, nColsWidth),
                            x=table$value)
      colnames(table) <- levels
      return(table)},
      BPPARAM=BPPARAM)
    
  }
  else
  {
    # setkeys for speed-up
    overlapTable[,V1:=as.integer(V1)]
    setkey(overlapTable, col_width, V1)
    
    # overlap with ref. coordinates
    if(is.null(aggregationFun)) error("Aggregation function needs to be defined")
    
    # why an lapply 05.04.24
    #overlapTable <- overlapTable[,.("scoreCol"=lapply(.SD, aggregationFun, na.rm=TRUE)),
    #                              by=c("col_width", "V1"),
    #                             .SDcols=scoreCol]
    setnames(overlapTable, scoreCol, "scoreCol")
    overlapTable <- overlapTable[,.(scoreCol=aggregationFun(scoreCol)),
                                 by=c("col_width", "V1")]
    
    # convert to sparse matrix
    overlapTable <- Matrix::sparseMatrix(i=overlapTable$V1, 
                                         j=overlapTable$col_width,
                                         dims=c(nRefs, nColsWidth),
                                         x=overlapTable$scoreCol)
    
    colnames(overlapTable) <- levels
  }
  
  gc()
  return(overlapTable)
}

getNonRedundantMotifs <- function(format=c("PFMatrix","universal","PWMatrix"),
                                  species=c("Hsapiens","Mmusculus")){
  species <- match.arg(species)
  motifs <- MotifDb::query(MotifDb::MotifDb, c(species,"HOCOMOCO"))
  pat <- paste0("^",species,"-HOCOMOCOv1[0-1]-|_HUMAN.+|_MOUSE.+|core-[A-D]-|secondary-[A-D]-")
  modf <- data.frame(row.names=names(motifs),
                     TF=gsub(pat,"",names(motifs)),
                     grade=gsub(".+\\.","",names(motifs)))
  modf <- modf[order(modf$TF,-as.numeric(grepl("HOCOMOCOv11",row.names(modf))),modf$grade),]
  modf <- modf[!duplicated(modf$TF),]
  motifs <- motifs[row.names(modf)]
  switch(match.arg(format),
         universal=setNames(universalmotif::convert_motifs(motifs), modf$TF),
         PFMatrix=do.call(TFBSTools::PFMatrixList, setNames(
           universalmotif::convert_motifs(motifs, class="TFBSTools-PFMatrix"),
           modf$TF)),
         PWMatrix=do.call(TFBSTools::PWMatrixList, 
                          setNames(universalmotif::convert_motifs(motifs, 
                                                                  class="TFBSTools-PWMatrix"), modf$TF))
  )
}


FQnorm <- function(counts, type="mean"){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  if(type=="mean"){
    # refdist <- apply(counts.sort,1,mean)
    refdist <- base::rowMeans(counts.sort)
  } else if(type=="median"){
    #refdist <- apply(counts.sort,1,median)
    refdist <- matrixStats::rowMedians(counts.sort)
  }
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}

gcqn <- function(counts, gcGroups, summary='mean', round=TRUE){
  gcBinNormCounts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts), 
                            dimnames=list(rownames(counts),colnames(counts)))
  for(ii in 1:nlevels(gcGroups)){
    id <- which(gcGroups==levels(gcGroups)[ii])
    if(length(id) == 1){
      normCountBin <- counts[id,]
      if(round) normCountBin <- round(normCountBin)
      gcBinNormCounts[id,] <- normCountBin
      next
    }
    countBin <- counts[id,,drop=FALSE]
    normCountBin <- FQnorm(countBin, type=summary)
    if(round) normCountBin <- round(normCountBin)
    normCountBin[normCountBin<0] <- 0
    gcBinNormCounts[id,] <- normCountBin
  }
  return(gcBinNormCounts)
}

gcqn_qsmooth <- function(counts, gcGroups, bio, round=TRUE){
  gcBinNormCounts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts), 
                            dimnames=list(rownames(counts),colnames(counts)))
  for(ii in 1:nlevels(gcGroups)){
    id <- which(gcGroups==levels(gcGroups)[ii])
    countBin <- counts[id,]
    qs <- qsmooth::qsmooth(countBin, group_factor=bio)
    normCountBin <- qs@qsmoothData
    if(round) normCountBin <- round(normCountBin)
    normCountBin[normCountBin<0L] <- 0L
    gcBinNormCounts[id,] <- normCountBin
  }
  return(gcBinNormCounts)
}

#' wrapper function by Jiayi Wang
GCQuantileNorm <- function(se, genome=NULL, g=20, 
                           summary='median', round=TRUE) {
  
  stopifnot(!is.null(genome) || !is.null(rowData(se)$bias))
  # get GC content
  gr <- rowRanges(se)
  if(!is.null(rowData(se)$bias)){
    gcContent <- rowData(se)$bias
  }else{
    peakSeqs <- getSeq(x = genome, gr)
    gcContent <- letterFrequency(peakSeqs, "GC", as.prob = TRUE)[,1]
  }
  gcGroups <- Hmisc::cut2(gcContent, g = g)
  
  # run QC quantile normalization
  countsGCQN <- gcqn(counts(se), gcGroups, summary = summary, round = round)
}

#' wrapper function by Jiayi Wang
GCSmoothQuantile <- function(se, genome=NULL, g=20, bio) {
  library(qsmooth)
  # get GC content
  gr <- rowRanges(se)
  se <- as(se,"SummarizedExperiment")
  if(!is.null(rowData(se)$bias)){
    gcContent <- rowData(se)$bias
  }else{
    peakSeqs <- getSeq(x = genome, gr)
    gcContent <- letterFrequency(peakSeqs, "GC", as.prob = TRUE)[,1]
  }
  gcGroups <- Hmisc::cut2(gcContent, g = g)
  
  # run QC smooth quantile normalization
  # colData(se)$group_id <- factor(substr(colnames(se), 1, 
  #     nchar(colnames(se)) - 5))
  
  # countsGCSQ <- gcqn_qsmooth(counts(se), gcGroups,
  #        bio=droplevels(colData(se)$group_id))
  countsGCSQ <- gcqn_qsmooth(se, gcGroups,
                             bio=droplevels(as.factor(colData(se)[,bio])))
}