setSechmOption("hmcols", value=colorspace::diverging_hcl(101,"Blue-Red 3"))

#' Perform staged-FDR for DEA interaction terms using the stageR package
#'
#' @param se A SummarizedExperiment object
#' @param screen A global DEA data.frame containing the screening p-values in
#'   the `PValue` column
#' @param fullFit An edgeR `glmFit` object with all the interaction terms
#' @param coefs The interaction coefficients to test (default all)
#'
#' @return The `se` object with additional rowData columns containing the
#'   adjusted p-values
stagedInteractions <- function(se, screen, fullFit, coefs=NULL){
  pScreen <- setNames(screen$PValue, row.names(screen))
  if(is.null(coefs))
    coefs <- grep(":",colnames(fullFit$coefficients),value=TRUE)
  if(is.null(names(coefs))) names(coefs) <- make.names(coefs)
  deas <- lapply(coefs, FUN=function(co){
    as.data.frame(topTags(glmLRT(fullFit, co), Inf))
  })
  pConf <- as.data.frame(lapply(deas, FUN=function(x) x[names(pScreen),"PValue"]))
  row.names(pConf) <- names(pScreen)
  stageRObj <- stageR(pScreen=pScreen, pConfirmation=as.matrix(pConf), pScreenAdjusted=FALSE)
  stageRObj <- stageWiseAdjustment(object=stageRObj, method="none", alpha=0.05)
  res <- as.data.frame(getAdjustedPValues(stageRObj, onlySignificantGenes=FALSE, order=FALSE))
  res <- res[,-1,drop=FALSE]
  for(f in colnames(res)){
    deas[[f]][row.names(res),"FDR"] <- res[,f]
    rowData(se)[[paste0("DEA.inter.",f)]] <- deas[[f]][row.names(se),]
  }
  if(ncol(res)>1){
    rowData(se)$int.minFDR <- 1
    rowData(se)[row.names(res),"int.minFDR"] <- matrixStats::rowMins(as.matrix(res))
  }
  se
}

lfcComp <- function(se, compared="Sex", splitBy=NULL){
  fcs <- scuttle::sumCountsAcrossCells(assays(se)$log2FC, colData(se)[,c("Sex","TimePoint")], average=TRUE)
  colnames(fcs) <- paste(fcs$Sex,fcs$TimePoint)
  fcs <- fcs[adegs,]
  fcs <- reshape2:::melt(assay(fcs), value.name="log2FC", factorsAsStrings=FALSE)
  fcs$Sex <- as.factor(gsub(" .+","",fcs$Var2))
  fcs$TimePoint <- factor(gsub("^.+ ","",fcs$Var2), levels(se$TimePoint))
  fcs <- fcs[fcs$TimePoint!="0min",]
  fcs2 <- do.call(rbind, lapply(split(fcs, droplevels(fcs$TimePoint)), FUN=function(x){
    dea <- as.character(x$TimePoint[1])
    dea <- intersect(c(dea,
                       paste0("inter.Sexmale.TimePoint",dea)), names(deas))
    degs <- lapply(deas[dea],FUN=function(dea){
      row.names(dea)[which(abs(dea$logFC)>th.lfc & dea$FDR<0.05)]
    })
    int.degs <- c()
    if(length(dea)>1) int.degs <- degs[[2]]
    x <- x[which(x$Var1 %in% unique(unlist(degs))),]
    y <- x[x$Sex=="male",]
    y <- y[order(y$log2FC),]
    g2r <- setNames(seq_along(y$Var1), y$Var1)
    x$rank <- g2r[as.character(x$Var1)]
    x$sig.int <- x$Var1 %in% int.degs
    x
  }))
}

agg4Ribb <- function(d, value.var, by, sem=TRUE){
  dag <- aggregate(d[,value.var,drop=FALSE], d[,by,drop=FALSE], FUN=mean)
  ds <- aggregate(d[,value.var,drop=FALSE], d[,by,drop=FALSE], FUN=function(x){
    if(sem) return(sd(x)/sqrt(length(x)))
    sd(x)
  })
  dag$lower <- dag[[value.var]]-ds[[value.var]]
  dag$upper <- dag[[value.var]]+ds[[value.var]]
  dag
}

#' getMsigSets
#'
#' Retrieves MsigDB genesets.
#'
#' @param species Default "Homo sapiens"; see `msigdbr::msigdbr_show_species()` for available species.
#' @param collections Collections to fetch, default H, C2 and C5.
#' @param subcats Sub categories to fetch
#'
#' @return A names list of genesets
#' @export
getMsigSets <- function(species="Mus musculus", collections=c("H","C2","C5"),
                        subcats=c("CP:KEGG", "GO:BP", "GO:CC","CP:WIKIPATHWAYS")){
  library(msigdbr)
  m <- msigdbr(species)
  m <- m[which(m$gs_cat %in% collections),]
  m <- m[which(m$gs_subcat %in% subcats),]
  m$name2 <- paste0(m$gs_subcat, ":", m$gs_name)
  split(m$gene_symbol, m$name2)
}



#' breakStrings
#'
#' breaks a string of words (or vector thereof) into two lines
#'
#' @param x a character vector
#' @param minSizeForBreak the minimum number of characters to break on two lines (default 20)
#' @param lb the line break character (default "\n")
#'
#' @return a character vector of length=length(x)
#'
#' @export
breakStrings <- function(x, minSizeForBreak=20, lb="\n"){
  sapply(x,minSizeForBreak=minSizeForBreak,lb=lb,FUN=function(x,minSizeForBreak,lb){
    if(nchar(x)<=minSizeForBreak)	return(x)
    g <- gregexpr(" ", x)[[1]]
    if(length(g)==0) return(x)
    if(length(g)==1 & all(g==-1)) return(x)
    mid <- nchar(x)/2
    mid <- g[order(abs(g-mid))[1]]
    substr(x, mid, mid) <- lb
    return(x)
  })
}

#' Prepares the collectTRI object into a regulon
#'
#' @param reg The collectTRI object
#' @param e The expression dataset or list of features included
#' @param cor.thres Correlation threshold for removal of colinear features
#' @param minsize Minimum regulon size
#' @param discardUndirected Logical; whether to discard targets where it's
#'   unclear whether they're activated or repressed by the TF
#'
#' @return A filtered network (data.frame of edges)
prepareRegulon <- function(reg, e, cor.thres=0.95, minsize=5, discardUndirected=FALSE){
  if(!is.character(e)) e <- row.names(e)
  if(discardUndirected){
    reg <- reg[which(reg$mor!=0)]
  }else{
    reg$mor[which(reg$mor==0)] <- 0.25
  }
  if(!is.null(reg$likelihood)){
    reg$mor <- reg$mor * reg$likelihood
    reg$likelihood <- NULL
  }
  # remove colinearity
  m <- reshape2::dcast(reg, target~source, value.var = "mor", fill = 0)
  row.names(m) <- m[,1]
  m <- as.matrix(m[,-1])
  m <- m[intersect(row.names(m), e),]
  m <- m[,colSums(m!=0)>=minsize]
  cc <- cor(m)
  cc[upper.tri(cc, diag=TRUE)] <- NA_real_
  removed <- vector("character")
  while(nrow(w <- which(cc>cor.thres, arr.ind=TRUE))>0){
    removed <- c(removed, row.names(w))
    keep <- setdiff(row.names(cc),row.names(w))
    cc <- cc[keep,][,keep]
  }
  message(paste0(ncol(cc),"/",length(unique(reg$source))," regulons kept."))
  if(length(removed)>0)
    message("The following factors were removed due to collinearity with other factors:\n",
            paste(removed, collapse=", "))
  reg[reg$source %in% row.names(cc),]
}

chipAtlasProcTargets <- function(d, patterns=c("Neuron","Hippocam","Cortex","Cerebellum")){
  suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(Matrix)
  })
  if(is.character(d)){
    if(is.null(names(d))) names(d) <- gsub("\\.[0-9]+\\.tsv$","",basename(d))
    d <- lapply(d, row.names=1, FUN=read.delim)
  }
  wp <- lapply(setNames(names(d),names(d)), FUN=function(n){
    d <- d[[n]]
    aw <- lapply(setNames(patterns, patterns), FUN=function(p){
      grep(p, colnames(d), ignore.case=TRUE)
    })
    aw <- aw[which(lengths(aw)>0)]
    wt <- "unspecific"
    if(length(aw)==0) return(list("unspecific",ncol(d)-2,d[,1,drop=FALSE]))
    d <- d[,aw[[1]],drop=FALSE]
    d2 <- as.data.frame(rowMeans(d))
    colnames(d2) <- n
    return(list(names(aw)[1],ncol(d),d2))
  })
  cd <- dplyr::bind_rows(lapply(wp, FUN=function(x) data.frame(type=x[[1]], nDatasets=x[[2]])))
  row.names(cd) <- names(wp)
  g <- unique(sort(unlist(lapply(wp,FUN=function(x) row.names(x[[3]])))))
  e <- sapply(wp, FUN=function(x) round(x[[3]][g,]))
  e[is.na(e)] <- 0L
  row.names(e) <- g
  se <- SummarizedExperiment(list(score=as(e, "sparseMatrix")), colData=cd[colnames(e),])
  assays(se)$logscore <- log1p(assay(se))
  assays(se)$norm <- sqrt(t(t(assay(se))/colMaxs(assay(se))))
  se
}

decoupler2SE <- function(res, colData=DataFrame(), metadata=list()){
  TFs <- sort(unique(res$source))
  al <- lapply(split(res, res$statistic), FUN=function(x){
    a <- reshape2::dcast(x, formula=source~condition, value.var="score")
    row.names(a) <- a$source
    as.matrix(a[TFs,-1])
  })
  SummarizedExperiment(al, colData=colData[colnames(al[[1]]),], metadata=metadata)
}

dampingPlot <- function(dd, goi, bins=100, position=position_dodge(width=0.2), ...){
  p <- ggplot(dd, aes(rank, damping, colour=category, label=gene)) + geom_point() +
    theme_classic() + theme(legend.position = "bottom") + labs(colour="Damping category") +
    ggrepel::geom_text_repel(data=dd[goi,], max.overlaps = 100, colour="black", position=position,
                             min.segment.length=0, ...) + coord_cartesian(xlim=c(0,1))
  ggExtra::ggMarginal(p, type="histogram", margins="y", groupFill=TRUE, bins=bins, colour=NA, alpha=1)
}


getDampingCoef <- function(normal, tampered, fdr.normal=0.05, fdr.tampered=0.5){
  d <- merge(normal, tampered, by="row.names", suffix=c("",".tampered"))
  d$damping <- ifelse(d$logFC>0,
                      1-(exp(d$logFC+d$logFC.tampered)-1)/(exp(d$logFC)-1),
                      1-(exp(-(d$logFC+d$logFC.tampered))-1)/(exp(-d$logFC)-1))
  d$damping <- pmax(pmin(d$damping,1),0)
  d$damping[which(d$FDR>fdr.normal)] <- NA
  d$damping[which(d$FDR.tampered>fdr.tampered & d$damping>0.5)] <- NA
  d
}


# ref=which(d$Experiment=="RestraintTimeline")
mdsTimeProj <- function(se, degs, ref=seq_len(ncol(se))){
  a <- assays(se)$corrected[intersect(degs,row.names(se)),]
  mds <- limma::plotMDS(a, top=length(degs), plot=FALSE)
  d <- cbind(x=mds$x, y=mds$y, as.data.frame(colData(se)), isRef=FALSE)
  d$isRef[ref] <- TRUE
  d$mds1 <- cmdscale(dist(t(a)),k=1)
  # center
  agg <- aggregate(d[ref,1:2], by=d[ref,c("TimePoint"), drop=FALSE], FUN=median)
  center <- optim(c(x=0,y=0), function(p){
    sum((agg[,2]-p[1])^2+(agg[,3]-p[2])^2)
  })$par
  d$x <- d$x-center[1]
  d$y <- d$y-center[2]
  d$rad <- atan2(d$y,d$x)
  fit <- princurve::principal_curve(as.matrix(d[,1:2]))
  d$lambda <- fit$lambda
  agg <- aggregate(d[ref,c("x","y","rad","lambda")], by=d[ref,c("TimePoint"), drop=FALSE], FUN=median)
  ctl <- ifelse(any(d$Treatment=="Handling"),"Handling","Handling +ARS")
  p.rad <- p.adjust(sapply(split(d, d$TimePoint), FUN=function(d){
    wilcox.test(d$rad[which(d$Treatment==ctl)],
                d$rad[which(d$Treatment!=ctl)])$p.value
    }))
  p.lambda <- p.adjust(sapply(split(d, d$TimePoint), FUN=function(d){
    wilcox.test(d$lambda[which(d$Treatment==ctl)],
                d$lambda[which(d$Treatment!=ctl)])$p.value
  }))
  list(d=d, agg=agg, padj.rad=p.rad, padj.lambda=p.lambda)
}

mdsTimePlot <- function(se, degs, hull=TRUE, arrows=TRUE, labels=TRUE, concavity=3){
  dds <- calcNormFactors(DGEList(assay(se), group = se$TimePoint))
  mds <- limma::plotMDS(dds[intersect(degs,row.names(dds)),], top=length(degs), plot=FALSE)
  d <- cbind(x=mds$x, y=mds$y, as.data.frame(colData(se)))
  tpcols <- metadata(se)$anno_colors$TimePoint
  agg <- aggregate(d[,1:2], by=d[,c("TimePoint"), drop=FALSE], FUN=median)
  agg2 <- agg[rep(1:nrow(agg), each=2)[c(-1,-2*nrow(agg))],]
  agg2$TimePoint <- c(agg2$TimePoint[-1], agg2$TimePoint[nrow(agg2)])
  agg2$color <- tpcols[agg2$TimePoint]
  p <- ggplot(d, aes(x,y))
  if(hull) p <- p + ggforce::geom_mark_hull(aes(fill=TimePoint), colour=NA, alpha=0.25,
                      concavity=concavity, radius=unit(5,"mm"), expand=unit(5, "mm"))
  if(arrows) p <- p + geom_path(data=agg2, aes(group=TimePoint), colour=agg2$color,
                                size=2, alpha=0.5, arrow=grid::arrow())
  p <- p + geom_point(aes(colour=TimePoint, shape=Sex), size=2.5) + theme_bw()
  if(labels) p <- p +
    ggrepel::geom_text_repel(data=agg, aes(label=TimePoint), force=4, min.segment.length=0,
                             show.legend=FALSE, nudge_y=c(-0.12,rep(0,3),0.06,0.1))
  p + scale_colour_manual(values=tpcols) + scale_fill_manual(values=tpcols) +
    labs(x=paste0("logFC dim 1 (",round(100*mds$var.explained[1]),"%)"),
         y=paste0("logFC dim 2 (",round(100*mds$var.explained[2]),"%)"))
}


plotGeneProf <- function(se, g, ...){
  d <- meltSE(se, g)
  d$feature <- factor(d$feature, unique(d$feature))
  dag <- aggregate(d[,"log2FC",drop=FALSE], by=d[,c("feature","Treatment","TimePoint")], FUN=median)
  dag <- agg4Ribb(d, "log2FC", by=c("Treatment","feature","TimePoint"), sem=TRUE)
  dag <- dag[which(!(dag$TimePoint %in% c("24h"))),]
  dag$Time <- c(0,45,90,180,240,5.5*60,24*60)[as.integer(dag$TimePoint)]
  if(any(dag$Treatment=="Handling")){
    crscols <- c(Handling = "grey", `10days CRS` = "blue3", `20days CRS` = "midnightblue")
  }else{
    crscols <- c("Handling +ARS" = "lightgray", `10days CRS +ARS` = "slateblue", `20days CRS +ARS` = "midnightblue")
  }

  ggplot(dag, aes(Time, log2FC, fill=Treatment)) +
    geom_hline(yintercept=0, linetype="dashed", colour="grey") +
    geom_line(aes(colour=Treatment), linewidth=1.1) + 
    geom_point(aes(colour=Treatment), size=2.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    facet_wrap(~feature, scales="free_y", ...) +
    labs(x="Time (min)", y="log2(fold-change)") + theme_bw() +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    scale_fill_manual(values=crscols, guide="none") +
    scale_color_manual(values=crscols) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  
}



ORA <- function(x, bg, gsets, minSize=5, excludeUnannotated=TRUE, reportOverlapping=TRUE){
  if(!is.list(gsets) && is.vector(gsets)){
    gsets <- list(set=gsets)
  }else if(excludeUnannotated){
    bg <- intersect(bg, unique(unlist(gsets)))
  }
  gsets <- lapply(gsets, unique)
  x <- unique(intersect(x,bg))
  gsets <- lapply(gsets, y=bg, FUN=intersect)
  gsets <- gsets[which(lengths(gsets)>=minSize)]
  res <- data.frame(Term=names(gsets), setSize=lengths(gsets))
  res$overlap <- sapply(gsets, y=x, FUN=function(x,y){ length(intersect(x,y)) })
  expected <- lengths(gsets)*length(x)/length(bg)
  res$log2.enrichment <- round(log2(res$overlap/expected),2)
  res$pValue <- fisher.test.p(res$overlap, lengths(gsets)-res$overlap,
                              length(x)-res$overlap, 
                              length(bg)-lengths(gsets)-length(x)+res$overlap,
                              alternative="greater")
  res <- res[which(res$overlap>0),]
  res$FDR <- p.adjust(res$pValue)
  if(reportOverlapping){
    res$overlap.items <- sapply(gsets[row.names(res)], y=x, FUN=function(x,y){
      paste(intersect(x,y),collapse=", ")
    })
  }
  row.names(res) <- NULL
  res[order(res$pValue, -res$log2.enrichment),]
}

#' fisher.test.p
#' 
#' Fast p-values from multiple Fisher's exact tests
#'
#' @param a,b,c,f Vectors containing the four positions of the 2x2 matrix
#' @param alternative greater, less, or 'two.sided' (default)
#'
#' @return A vector of p-values
#' @importFrom stats dhyper
#' @export
fisher.test.p <- function (a, b, c, d, 
                           alternative=c("two.sided", "greater", "less")){
  fn <- switch( match.arg(alternative), 
                less = function(x,m,n,k) phyper(x, m, n, k), 
                greater = function(x,m,n,k) phyper(x - 1, m, n, k, lower.tail=FALSE), 
                two.sided = function(x,m,n,k){
                  lo <- max(0, k - n)
                  support <- seq(from=lo, to=min(k, m))
                  d <- dhyper(support, m, n, k, log = TRUE)
                  d <- exp(d - max(d))
                  d <- d/sum(d)
                  sum(d[d <= d[x - lo + 1] * (1 + 10^(-7))])
                }
  )
  mapply(FUN=fn, a, a+c, b+d, a+b)
}

updateSEforPlots <- function(se, dark=FALSE){
  
  se$Treatment <- factor(se$Treatment)
  if(any(grepl("days",levels(se$Treatment)))){
    crscols <- c("Handling +ARS" = "lightgray", `10days CRS +ARS` = "slateblue", `20days CRS +ARS` = "midnightblue")
    levels(se$Treatment) <- gsub(" +ARS +ARS"," +ARS", paste0(levels(se$Treatment)," +ARS"), fixed=TRUE)
  }else{
    crscols <- c("Handling +ARS"="darkgrey", "CRS +ARS"=ifelse(dark,"midnightblue","#0000CD96"))
    levels(se$Treatment) <- names(crscols)
  }
  if(!is.null(metadata(se)$anno_colors$Treatment))
    metadata(se)$anno_colors$Treatment <- crscols
  se
}

