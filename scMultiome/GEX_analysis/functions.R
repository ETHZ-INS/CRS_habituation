#' Applies IHW with bulk significance as covariate to muscat pseudobulk 
#' analysis
#'
#' @param res 
#' @param bulkfdr A named vector of bulk significances
#' @param pb An optional pseudobulk object, in which the DEAs will be added
#'   as rowData columns.
#' @param alpha The aimed alpha for IHW
#'
#' @return A named list with the updated `res` and, if provided, `pb`
applyIHWonMuscatRes <- function(res, bulkfdr, pb=NULL, alpha=0.1){
  for(comp in names(res$table)){
    d <- dplyr::bind_rows(res$table[[comp]])
    d$coef <- comp
    d$bulkFDR <- bulkfdr[d$gene]
    d$bulkFDR[is.na(d$bulkFDR)] <- 2
    a <- IHW::ihw(d$p_val, d$bulkFDR, alpha=alpha)
    d$ihw <- IHW::adj_pvalues(a)
    res$table[[comp]] <- d <- split(d,d$cluster_id)
    if(!is.null(pb)){
      for(ct in names(d)){
        x <- d[[ct]]
        row.names(x) <- x$gene
        x <- x[,intersect(c("logFC","logCPM","p_val","p_adj.loc","p_adj.glb"),colnames(x))]
        rowData(pb)[[paste0("DEA.",ct,".",comp)]] <- x[row.names(pb),]
      }
    }
  }
  ret <- list(res=res)
  if(!is.null(pb)) ret$pb <- pb
  ret
}

#' computeResponseCoef
#'
#' @param d A data.frame of DEA results, including the columns gene, cluster_id,
#'   logFC and that specified by `sig`. The data.frame should be filtered for
#'   genes of the response.
#' @param pb An optional pseudobulk SummarizedExperiment with a log2FC assay, 
#'   and a celltype that corresponds to the `cluster_id` column of `d`. Only 
#'   needed for sample-wise coefs.
#' @param sampleWise Logical, whether to obtain (pseudobulk) sample-wise 
#'   coefficients. Enabled by default if `pb` is provided.
#' @param sig The column of `d` to use as significance.
#' @param tryRobust Logical, whether to try to use a robust fit
#'
#' @return A matrix of estimated coefficients and their significance for each
#'   sample or celltype.
computeResponseCoef <- function(d, pb=NULL, sampleWise=!is.null(pb), sig="ihw", tryRobust=FALSE,
                                linearScale=TRUE){
  d <- d[!is.na(d$logFC),]
  d[[sig]][is.na(d[[sig]])] <- 1
  lfcs <- reshape2::dcast(d, gene~cluster_id, value.var = "logFC", fun.aggregate = mean)
  geomean <- function(p) exp(mean(log(p)))
  fdrs <- reshape2::dcast(d, gene~cluster_id, value.var = "ihw", fun.aggregate = min)
  row.names(lfcs) <- row.names(fdrs) <- lfcs[,1]
  lfcs <- as.matrix(lfcs[,-1])
  fdrs <- as.matrix(fdrs[,-1])
  fdrs[is.na(fdrs)] <- 1
  weights <- (1-rowMeans(fdrs)) * (1-fdrs)
  meanlfcs <- rowSums(lfcs*(1-fdrs))/rowSums(1-fdrs)
  if(linearScale){
    lfcs <- exp(lfcs)
    meanlfcs <- exp(meanlfcs)
  }
  w <- which(!is.na(meanlfcs))
  if(!sampleWise){
    coef <- t(sapply(setNames(seq_len(ncol(lfcs)),colnames(lfcs)), FUN=function(i){
      mod <- NULL
      tryCatch({
        if(tryRobust)
          mod <- try(MASS::rlm(meanlfcs[w], lfcs[w,i], weights=weights[w,i]), silent=TRUE)
        if(!tryRobust || is(mod,"try-error") || !mod$converged){
          mod <- lm(lfcs[w,i]~0+meanlfcs[w], weights=weights[w,i])
        }
        coef(summary(mod))[1,1:ifelse(tryRobust,3,4)]
      }, error=function(e){
        rep(NA_real_,ifelse(tryRobust,3,4))
      })
    }))
    return(coef)
  }
  stopifnot(inherits(pb, "SummarizedExperiment"))
  lfcs <- assays(pb)$log2FC[names(which(!is.na(meanlfcs))),]
  if(linearScale) lfcs <- exp(lfcs)
  t(sapply(setNames(seq_len(ncol(pb)),colnames(pb)), FUN=function(i){
    iw <- weights[w,as.character(pb$celltype[i])]
    mod <- NULL
    tryCatch({
      if(tryRobust) 
        mod <- try(MASS::rlm(meanlfcs[w], lfcs[,i], weights=iw), silent=TRUE)
      if(!tryRobust || is(mod,"try-error") || !mod$converged){
        mod <- lm(lfcs[,i]~0+meanlfcs[w], weights=iw)
      }
      x <- coef(summary(mod))[1,1:ifelse(tryRobust,3,4)]
    }, error=function(e){
      rep(NA_real_,ifelse(tryRobust,3,4))
    })
  }))
}


#' computeCellResponseCoef
#'
#' @param d A data.frame of DEA results, including the columns gene, cluster_id,
#'   logFC and that specified by `sig`. The data.frame should be filtered for
#'   genes of the response.
#' @param sce A SingleCellExperiment with a `celltype` colData column 
#'   corresponding to the `cluster_id` of `d`.
#' @param ctrlPb A pseudobulk of controls SummarizedExperiment with a `celltype` 
#'   colData column corresponding to the `cluster_id` of `d`.
#' @param sig The column of `d` to use as significance.
#' @param ctfield The column of `ctrlPb` and `sce` with the celltype
#' @param indep.lfc The type of logFC to use as indenpendent variable. Either 
#'   'mean' (mean across cell types), 'celltype' (specific for that cell type),
#'   or 'weighted' (weighted average of the two, the default).
#' @param useintercept Whether to include an intercept in the model, recommended
#'   not.
#' @param robust Whether to estimate priors robustly.
#' @param weightedFit Whether to weigh observations, for fitting, by their 
#'   -log10(FDR) in the cell type of interest.
#'
#' @return A vector of estimated coefficients for each cell.
computeCellResponseCoef <- function(d, sce, ctrlPb, sig="ihw", ctfield="celltype",
                                    indep.lfc=c("weighted","mean","celltype"),
                                    useintercept=FALSE, robust=TRUE,
                                    weightedFit=FALSE, pseudocount=1){
  stopifnot(inherits(sce, "SingleCellExperiment"))
  stopifnot(inherits(ctrlPb, "SummarizedExperiment"))
  indep.lfc <- match.arg(indep.lfc)
  d <- d[!is.na(d$logFC),]
  d[[sig]][is.na(d[[sig]])] <- 1
  g <- unique(d$gene)
  d <- d[which(d$cluster_id %in% unique(sce[[ctfield]])),]
  
  lfcs <- reshape2::dcast(d, gene~cluster_id, value.var = "logFC", fun.aggregate = mean)
  geomean <- function(p) exp(mean(log(p)))
  fdrs <- reshape2::dcast(d, gene~cluster_id, value.var = sig, fun.aggregate = min)
  row.names(lfcs) <- row.names(fdrs) <- lfcs[,1]
  lfcs <- as.matrix(lfcs[,-1,drop=FALSE])
  fdrs <- as.matrix(fdrs[,-1,drop=FALSE])
  fdrs[is.na(fdrs) | is.infinite(fdrs)] <- 1
  lfcs[is.na(lfcs)] <- 0
  wg <- which(rowSums(fdrs==1)!=ncol(fdrs))
  if((dropped <- length(g)-length(wg))>0)
    message("Dropping ", dropped, " gene(s)")
  fdrs <- fdrs[wg,]
  lfcs <- lfcs[wg,]
  g <- g[wg]
  fdrw <- 1-sqrt(fdrs)
  # pbt <- ctrlPb[,ctrlPb[[ctfield]]==ct]
  # pbt <- pbt[colMeans(t(assay(pbt))/pbt$ncells)>0.5,]
  # 
  meanlfcs <- rowSums(lfcs*fdrw)/rowSums(fdrw)
  e <- 10000*t(as.matrix(counts(sce)[g,]))/sce$sum
  ac <- unlist(lapply(intersect(unique(ctrlPb[[ctfield]]), unique(d$cluster_id)),
                      FUN=function(ct){
    w <- ctrlPb[[ctfield]]==ct
    wCtl <- which(sce[[ctfield]]==ct & sce$TimePoint=="Control")
    ctrl <- colMeans(10000*t(pseudocount+as.matrix(counts(sce)[g,wCtl]))/sce$sum[wCtl])
    #ctrl <- colMeans(10000*t(assay(ctrlPb)[g,w])/colSums(assay(ctrlPb))[w])
    wsc <- colnames(sce)[sce[[ctfield]]==ct]
    es <- counts(sce)[g,wsc]
    #es <- 10000*t(1L+es)/colData(sce)[wsc,"sum"]
    #lfc <- log2(t(es))-log2(ctrl)
    x <- switch(indep.lfc,
                weighted=rowMeans(cbind(lfcs[,ct]*(1-fdrs[,ct]),
                                        meanlfcs*(fdrs[,ct])), na.rm=TRUE),
                mean=meanlfcs,
                celltype=lfcs[,ct]
    )
    if(useintercept){
      mm <- model.matrix(~1+x)
    }else{
      mm <- matrix(x, dimnames=list(NULL, "x"))
    }
    mm <- cbind(log(10000*ctrl/ctrl.sum), x2)
    res <- pbapply::pbapply(es,2,FUN=function(y){
      coef(glm.fit(mm, y, family=poisson(), intercept=FALSE, start=c(1,0)))
    })
    if(weightedFit){
      fit <- lmFit(lfc, mm, weights=1-sqrt(fdrs[,ct]))
      #pmin(2,pmax(0.1,-log10(fdrs[,ct]))))
    }else{
      fit <- lmFit(t(lfc), mm)
    }
    res <- topTable(eBayes(fit, robust=robust), coef = "x", number=Inf)
    setNames(res$logFC, row.names(res))
  }))
  ac[colnames(sce)]
}


#' computeCellResponseCoef
#'
#' @param d A data.frame of DEA results, including the columns gene, cluster_id,
#'   logFC and that specified by `sig`. The data.frame should be filtered for
#'   genes of the response.
#' @param sce A SingleCellExperiment with a `celltype` colData column 
#'   corresponding to the `cluster_id` of `d`.
#' @param ctrlPb A pseudobulk of controls SummarizedExperiment with a `celltype` 
#'   colData column corresponding to the `cluster_id` of `d`.
#' @param sig The column of `d` to use as significance.
#' @param ctfield The column of `ctrlPb` and `sce` with the celltype
#' @param indep.lfc The type of logFC to use as indenpendent variable. Either 
#'   'mean' (mean across cell types), 'celltype' (specific for that cell type),
#'   or 'weighted' (weighted average of the two, the default).
#' @param useintercept Whether to include an intercept in the model, recommended
#'   not.
#' @param robust Whether to estimate priors robustly.
#' @param weightedFit Whether to weigh observations, for fitting, by their 
#'   -log10(FDR) in the cell type of interest.
#'
#' @return A vector of estimated coefficients for each cell.
computeCellResponseCoefGlm <- function(d, sce, wCtrl, ctfield="celltype", sig="ihw",
                                    indep.lfc=c("mean","weighted","celltype"),
                                    lfc.conditionAgg=c("absmax","mean"),
                                    method=c("glmGamPoi","poisson"),
                                    BPPARAM=SerialParam(progress=TRUE)){
  stopifnot(inherits(sce, "SingleCellExperiment"))
  method <- match.arg(method)
  indep.lfc <- match.arg(indep.lfc)
  d <- d[!is.na(d$logFC),]
  d[[sig]][is.na(d[[sig]])] <- 1
  g <- unique(d$gene)
  d <- d[which(d$cluster_id %in% unique(sce[[ctfield]])),]
  
  if(match.arg(lfc.conditionAgg)=="absmax"){
    am <- function(x){ if(length(x)==0) return(0); x[which.max(abs(x))] }
  }else{
    am <- mean
  }
  lfcs <- reshape2::dcast(d, gene~cluster_id, value.var = "logFC", fun.aggregate = am)
  fdrs <- reshape2::dcast(d, gene~cluster_id, value.var = sig, fun.aggregate = min)
  row.names(lfcs) <- row.names(fdrs) <- lfcs[,1]
  lfcs <- as.matrix(lfcs[,-1,drop=FALSE])
  fdrs <- as.matrix(fdrs[,-1,drop=FALSE])
  fdrs[is.na(fdrs) | is.infinite(fdrs)] <- 1
  lfcs[is.na(lfcs)] <- 0
  wg <- which(rowSums(fdrs==1)!=ncol(fdrs))
  if((dropped <- length(g)-length(wg))>0)
    message("Dropping ", dropped, " gene(s)")
  fdrs <- fdrs[wg,,drop=FALSE]
  lfcs <- lfcs[wg,,drop=FALSE]
  g <- g[wg]
  fdrw <- 1-sqrt(fdrs)
  #meanlfcs <- rowMeans(lfcs)
  meanlfcs <- rowSums(lfcs*fdrw)/rowSums(fdrw)
  
  rs <- assay(scuttle::sumCountsAcrossCells(assay(sce)[,wCtrl], sce[[ctfield]][wCtrl]))
  rss <- sapply(split(sce$sum, sce[[ctfield]]), sum)[colnames(rs)]
  
  e <- 10000*t(as.matrix(counts(sce)[g,]))/sce$sum
  ac <- bplapply(intersect(unique(sce[[ctfield]]), unique(d$cluster_id)),
                 BPPARAM=BPPARAM, FUN=function(ct){
    topG <- row.names(rs)[head(order(-rs[,ct]),1000)]
    g2 <- union(topG, g)
    ctrl <- 10000*rs[g2,ct]/rss[ct]
    wsc <- which(sce[[ctfield]]==ct)
    e <- counts(sce)[g2,wsc]
    x <- switch(indep.lfc,
                weighted=rowMeans(cbind(lfcs[,ct]*(1-fdrs[,ct]),
                                        meanlfcs*(fdrs[,ct])), na.rm=TRUE),
                mean=meanlfcs,
                celltype=lfcs[,ct]
    )
    x2 <- setNames(rep(0,length(g2)),g2)
    x2[names(x)] <- x
    mm <- cbind(scale(log1p(ctrl)), as.numeric(x2))
    if(method=="poisson"){
      fit.ctl <- glm.control(epsilon=1e-4, maxit=50)
      o <- apply(e,2,FUN=function(y){
        tryCatch(
          coef(glm.fit(mm, y, family=poisson(), intercept=FALSE, start=c(1,0)))[2],
          error=function(e) NA_real_)
      })
    }else{
      o <- glmGamPoi::glm_gp(t(as.matrix(e)), mm, overdispersion="global", size_factors=1)$Beta[,2]
    }
    setNames(o, colnames(sce)[wsc])
  })
  unlist(ac)[colnames(sce)]
}


qtScale <- function(x, qt=c(0.01,0.99), maxMin=NULL, minMax=NULL, by=NULL){
  if(!is.null(by)){
    if(flag <- is.null(names(x))) names(x) <- seq_along(x)
    no <- names(x)
    x <- split(x, by)
    names(x) <- NULL
    x <- lapply(x, FUN=qtScale, qt=qt, maxMin=maxMin, minMax=minMax)
    x <- unlist(x)[no]
    if(flag) names(x) <- NULL
    return(x)
  }
  qt <- quantile(x,qt,na.rm=TRUE)
  if(!is.null(maxMin)) qt[1] <- min(maxMin, qt[1])
  if(!is.null(minMax)) qt[2] <- max(minMax, qt[2])
  setNames(pmin(1,pmax(0,x-qt[1])/(qt[2]-qt[1])), names(x))
}

celltypeVolcanos <- function(resFlattened, ihw.thres=0.1, alfc.thres=log2(2), maxG=10){
  res3 <- resFlattened
  res3$significant <- res3$ihw<ihw.thres & abs(res3$logFC)>alfc.thres
  res3$celltype <- res3$cluster_id
  levels(res3$celltype) <- gsub("\\.|/","\n",levels(res3$celltype))
  
  res4 <- res3[which(res3$significant),]
  res4 <- dplyr::bind_rows(lapply(split(res4,paste0(res4$Comparison,res4$celltype)), FUN=\(x){
    x[head(order(x$ihw),maxG),]
  }))
  
  
  ggplot(res3, aes(logFC, -log10(ihw), colour=significant, label=gene)) +
    ggrastr::geom_point_rast() +
    ggrepel::geom_text_repel(data=res4, colour="black", min.segment.length = 0) +
    facet_grid(Comparison~celltype) + theme_bw() + 
    coord_cartesian(xlim=round(range(res4$logFC),1))
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

ctGenProf <- function(se, g, y="logcpm", sem=FALSE, scaleTime=FALSE){
  d <- meltSE(se, g, rowDat.columns = NA)
  d$feature <- factor(d$feature, unique(d$feature))
  dag <- agg4Ribb(d, y, by=c("celltype2","feature","TimePoint","Treatment"), sem=sem)
  if(y=="logcpm") dag$lower <- pmax(dag$lower,0)
  dag <- dag[which(!(dag$TimePoint %in% c("24h"))),]
  dag$Time <- c(0,15,45,180)[as.integer(dag$TimePoint)]
  tv <- ifelse(scaleTime, "Time", "TimePoint")
  p <- ggplot(dag, aes(.data[[tv]], .data[[y]], fill=Treatment, group=Treatment))
  if(y %in% c("scaledLFC","log2FC")) p <- p + geom_hline(yintercept=0, linetype="dashed", colour="grey")
  p +
    geom_line(aes(colour=Treatment), linewidth=1.1) + 
    geom_point(aes(colour=Treatment), size=2.5) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4) +
    facet_grid(feature~celltype2, scales="free_y") +
    labs(x=ifelse(scaleTime,"Time (min)","Time")) + theme_bw() +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    scale_fill_discrete(guide="none") +
    theme(axis.text.x=element_text(angle=90, hjust=1),
          legend.position = "bottom")
  
}

