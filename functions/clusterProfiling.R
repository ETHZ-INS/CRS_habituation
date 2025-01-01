
#' goclassify
#'
#' @param cl A named vector indicating the cluster for each gene, or a named list
#' containing the genes for each cluster
#' @param go A named list of genes for each geneset
#' @param minSize The minimum size for a geneset to be considered
#' @param max_depth passed to xgboost
#' @param N the number of sets to select
#' @param nrounds passed to xgboost
#' @param by how to select the features
#' @param onlyAnnotated Logical; whether to restrict elements of `cl` to those that are
#' present in some set of `go`
#'
#' @return a list
#'
#' @examples
#' go <- msigdbr::msigdbr("Mus musculus", "C5")
#' go <- go[which(go$gs_subcat=="BP"),]
#' go <- split(go$gene_symbol, go$gs_name)
#' cl <- lapply(c(c1=50, c2=100), FUN=function(x) sample(unique(unlist(go)),x))
#' res <- goclassify(cl, go)
goclassify <- function( cl, go, minSize=5, max_depth=5, N=30, nrounds=50,
                        by=c("Frequency","Gain","Coverage"), onlyAnnotated=TRUE,
                        onlyInSets=TRUE){
  if(!is.list(cl)) cl <- split(names(cl),cl)
  if(onlyAnnotated) cl <- lapply(cl, y=unique(unlist(go)), FUN=intersect)
  if(onlyInSets) go <- lapply(go, y=unique(unlist(cl)), FUN=intersect)
  by <- match.arg(by)
  bm <- sapply(go, y=unlist(cl), FUN=function(x,y) as.numeric(y %in% x))
  row.names(bm) <- unlist(cl)
  bm <- bm[,colSums(bm)>=minSize]
  # library(glmnet)
  # fits = cv.glmnet( bm2, as.factor(labs), family = "multinomial",
  #                   type.multinomial = "grouped", standardize=FALSE )
  # co <- coef(fits, fits$lambda.min)
  # co <- row.names(co)[co[,1]!=0][-1]
  #
  suppressPackageStartupMessages(library(xgboost))
  labs <- rep(as.integer(as.factor(names(cl))), sapply(cl,length))-1
  fit <- xgboost(bm, labs, params=list( booster="gbtree", max_depth=max_depth,
                                         subsample=0.75, colsample_bytree=1,
                                         objective="multi:softprob",
                                         eval_metric="mlogloss",
                                         min_child_weight=3,
                                         num_class=length(cl) ),
                    nrounds = nrounds, verbose=0)
  vi <- xgb.importance(model=fit)
  co <- as.character(vi$Feature[order(vi$Gain, decreasing=TRUE)[seq_len(N)]])
  fns <- list( overlap=function(x,y) length(intersect(x,y)),
               proportion=function(x,y) length(intersect(x,y))/length(x),
               enrichment=function(x,y){
                 expected <- length(y)*length(x)/nrow(bm)
                 length(intersect(x,y))/expected
               },
               jaccard=function(x,y) length(intersect(x,y))/length(union(x,y)),
               ocoef=function(x,y) length(intersect(x,y))/min(length(x),length(y))
             )
  m <- lapply(fns, FUN=function(fn){
    t(sapply(go[co],FUN=function(y){
      sapply(cl, y=y, FUN=fn)
    }))
  })
  m$lengths <- sapply(cl,length)
  m
}

plot.goclassify <- function(m, n=min(15, nrow(m$overlap[!is.na(rownames(m$overlap)),])),
                            annotation_legend_param=list(), transpose=FALSE,
                            what=c("proportion","log2enr","enrichment","jaccard"),
                            show_size=TRUE, ...){
  m1 <- m[1:5]
  m1 <- lapply(m1, FUN=function(x) x[!is.na(rownames(x)),])
  m1[["lengths"]] <- m$lengths
  m <- m1
  n <- min(n, nrow(m$overlap))
  an <- data.frame(row.names=names(m$lengths), size=as.numeric(m$lengths))
  library(circlize)
  col_fun = colorRamp2(c(0, max(an$size)), c("white", "grey"))
  an <- HeatmapAnnotation(df=an,annotation_legend_param=annotation_legend_param,col = list(size = col_fun),
                          which=ifelse(transpose,"row","column"))
  what <- match.arg(what)
  if(what=="log2enr"){
    m$log2enr <- log2(m$enrichment+0.1)
  }
  row.names(m[[what]]) <- gsub("^GO_","",row.names(m[[what]]))
  if(na.omit(!is.null(breakStrings)))
    row.names(m[[what]]) <- breakStrings(gsub("_"," ",row.names(m[[what]])))
  suppressPackageStartupMessages(library(ComplexHeatmap))
  if(transpose){
    h <- Heatmap( t(m[[what]][seq_len(n),]), right_annotation=an,
                    name=what, cell_fun=function(j, i, x, y, width, height, fill){
                      grid.text(m$overlap[j,i], x, y, gp=gpar(fontsize=10))
                    }, ...)
    return(h)
  }
  h <- Heatmap( m[[what]][seq_len(n),], top_annotation=an, cluster_columns = F,
                name=what, cell_fun=function(j, i, x, y, width, height, fill){
                  grid.text(m$overlap[i,j], x, y, gp=gpar(fontsize=10))
                }, ...)
  #p <- draw(h, heatmap_legend_side="bottom", merge_legends = T)
  return(h)
}


plotClusterProfiler <- function(m, n = 15, color_by = c("FDR","enrichment")){
  # input is the list split into dataframes of the different clusters
  # the dataframes are already ordered by increasing pvalues so we take the top n rows of all lists
  # GeneRatio is observed proportion of genes, BgRatio expected proportion of genes to overlap
  # calculate enrichment score for all terms
  m <- lapply(m, FUN=function(y){
    f <- function(x){ x <- as.numeric(x); x[1]/x[2] }
    Enrichment <- sapply(strsplit(y$GeneRatio, split = "/"), FUN=f)/
      sapply(strsplit(y$BgRatio, split = "/"), FUN=f)
    cbind(y,Enrichment)
  })
  # isolate top n terms
  desc <- lapply(m, FUN=function(x) x[x$p.adjust < 0.05,])
  desc <- lapply(desc, FUN=function(x) x[1:n,])
  desc <- unique(unlist(lapply(desc, FUN=function(x) x$Description)))[1:n]
  # Enrichment scores
  mat <- matrix(0, nrow = n, ncol = length(names(m)))
  rownames(mat) <- desc
  colnames(mat) <- names(m)
  for (j in (1:ncol(mat))) { # for each cluster or list element
    for (i in (1:nrow(mat))) { # for each description in each cluster
      if(rownames(mat)[i] %in% m[[j]]$Description == TRUE){
        mat[i,j] <- m[[j]]$Enrichment[m[[j]]$Description == rownames(mat)[i]]
      }
    }
  }
  # colour by FDR (qvalues)
  fdr <- matrix(0, nrow = n, ncol = length(names(m)))
  rownames(fdr) <- desc
  colnames(fdr) <- names(m)
  for (j in (1:ncol(fdr))) { # for each cluster or list element
    for (i in (1:nrow(fdr))) { # for each description in each cluster
      if(rownames(fdr)[i] %in% m[[j]]$Description == TRUE){
        fdr[i,j] <- m[[j]]$qvalue[m[[j]]$Description == rownames(fdr)[i]]
      }
    }
  }
  # Counts
  c <- matrix(0, nrow = n, ncol = length(names(m)))
  rownames(c) <- desc
  colnames(c) <- names(m)
  for (j in (1:ncol(c))) { # for each cluster or list element
    for (i in (1:nrow(c))) { # for each description in each cluster
      if(rownames(c)[i] %in% m[[j]]$Description == TRUE){
        c[i,j] <- m[[j]]$Count[m[[j]]$Description == rownames(c)[i]]
      }
    }
  }
  rownames(c) <- NULL
  colnames(c) <- NULL
  # heatmap
  suppressPackageStartupMessages(library(ComplexHeatmap))
  if(color_by == "qvalue"){
    # symmetric colour scale
    suppressPackageStartupMessages(library(circlize))
    col_fun = colorRamp2(c(0, max(fdr)), c("blue", "red"))
    h <- Heatmap(fdr, cluster_columns = FALSE, cluster_rows = TRUE, name = "qvalue", col = col_fun, column_names_rot = 0, cell_fun=function(j, i, x, y, width, height, fill){
      grid.text(c[i,j], x, y, gp=gpar(fontsize=10))
    })
  }
  if(color_by == "enrichment"){
    # symmetric colour scale
    suppressPackageStartupMessages(library(circlize))
    a <- max(c(abs(max(mat)),abs(min(mat))))
    col_fun = colorRamp2(c(-a, 0, a), c("blue", "white", "red"))
    h <- Heatmap(mat, cluster_columns = FALSE, cluster_rows = TRUE, name = "log2enr", col = col_fun, column_names_rot = 0, cell_fun=function(j, i, x, y, width, height, fill){
      grid.text(c[i,j], x, y, gp=gpar(fontsize=10))
    })
  }
  draw(h, heatmap_legend_side="bottom")
}


#' farthestPoint
#'
#' Identifies the point farthest from a line passing through by the first and
#' last points. Used for automatization of the elbow method.
#'
#' @param y Monotonically inscreasing or decreasing values
#' @param x Optional x coordinates corresponding to `y` (defaults to seq)
#'
#' @return The value of `x` farthest from the diagonal.
#' @export
#'
#' @examples
#' y <- 2^(10:1)
#' plot(y)
#' x <- farthestPoint(y)
#' points(x,y[x],pch=16)
farthestPoint <- function(y, x=NULL){
  if(is.null(x)) x <- seq_len(length(y))
  d <- apply( cbind(x,y), 1,
              a=c(1,y[1]), b=c(length(y),rev(y)[1]),
              FUN=function(y, a, b){
                v1 <- a-b
                v2 <- y-a
                abs(det(cbind(v1,v2)))/sqrt(sum(v1*v1))
              })
  order(d,decreasing=TRUE)[1]
}


#' clusterEnrichment
#'
#' @param clusters A named list of sets of genes of interest
#' @param sets A named list of reference genesets (must use the same gene
#' identifiers as `clusters`)
#' @param universe The optional universe (by default, the union of clusters is
#' used)
#' @param minSetSize The minimum set size (in universe) to consider
#' @param threshold The p-value threshold
#' @param family The model family for deviance calculation.
#'
#' @return A data.frame
clusterEnrichment <- function(clusters, sets, universe=NULL, minSetSize=10, threshold=0.05, family=c("binomial","poisson")){
  family <- match.arg(family)
  if(is.null(universe)) universe <- unlist(clusters)
  sets <- lapply(sets, y=universe, FUN=intersect)
  clusters <- lapply(clusters, y=universe, FUN=intersect)
  sets <- sets[sapply(sets,length)>=minSetSize]
  sl <- sapply(sets, length)
  cs <- sapply(clusters, length)
  universe <- length(universe)
  b <- sapply(clusters, FUN=function(x){
    sapply(sets, FUN=function(y){
      length(intersect(x,y))
    })
  })
  dev <- sapply(1:nrow(b),sz=cs,mod=family,FUN=function(g, sz, mod){
    x <- b[g,]
    expected <- sz*sl[g]/universe
    enr <- log1p(x)-log1p(expected)
    p=sum(x)/sum(sz)
    # calculation of
    if(mod=="binomial"){
      term1<-sum(x*log(x/(sz*p)), na.rm=TRUE)
      nx<-sz-x
      term2<-sum(nx*log(nx/(sz*(1-p))), na.rm=TRUE)
      dev <- 2*(term1+term2)
    }else{
      dev <- 2*sum(x*log(x/(sz*p)),na.rm=TRUE)-2*sum(x-sz*p)
    }
    pval <- pchisq(dev, length(x)-1, lower.tail=FALSE)
    c(deviance=dev, p.value=pval, FDR=NA_real_, enr)
  })
  dev <- as.data.frame(t(dev))
  dev$FDR <- p.adjust(dev$p.value, method="fdr")
  colnames(dev)[4:ncol(dev)] <- paste("enrichment",colnames(dev)[4:ncol(dev)],sep=".")
  row.names(dev) <- row.names(b)
  dev <- dev[dev$p.value<threshold,]
  dev[order(dev$p.value),]
}


plotClusterEnrichment <- function(dev, k=5, top=NULL,
                                  color=colorRampPalette(c("blue","black", "yellow"))(50),
                                  ...){
  dev <- dev[order(dev$p.value),]
  if(is.null(top)){
    d <- 1-cor(t(dev[,grep("enrichment",colnames(dev))]))
    bg <- igraph::cluster_louvain(knn.graph(d,k))
    dev$cluster <- NA_integer_
    dev[bg$names,"cluster"] <- bg$membership
    s1 <- row.names(dev)[unique(apply(dev[,grep("enrichment",colnames(dev))],2,which.max))]
    s2 <- row.names(dev)[!duplicated(dev$cluster)]
    s <- union(s1,s2)
  }else{
    s <- row.names(dev)[seq_len(min(top,nrow(dev)))]
  }
  deve <- dev[s,grep("enrichment",colnames(dev))]
  colnames(deve) <- gsub("enrichment\\.","",colnames(deve))
  ComplexHeatmap::pheatmap(deve, color=color, border_color = NA, ...)
}



clusterGS <- function(gs, gsets, allgenes=NULL, k=2:min(20,nrow(gs)-1), thres=0.05){
  library(Matrix)
  gs <- gs[gs$padj<thres,]
  k <- k[k>1 & k<nrow(gs)]
  if(length(k)==0){
    warning("Not enough sets")
    return(gs)
  }
  gsets <- gsets[gs$pathway]
  ag <- unique(unlist(gsets))
  if(!is.null(allgenes)) ag <- intersect(ag,allgenes)
  gsets <- lapply(gsets, y=ag, intersect)
  gsets <- FactorList(CharacterList(gsets))
  bm <- sparseMatrix(j=rep(seq_along(gsets),lengths(gsets)),
                     i=unlist(IntegerList(gsets)), index1=TRUE,
                     dims=c(length(ag),length(gsets)))
  colnames(bm) <- names(gsets)
  d <- as.dist(1-jaccard(bm))
  cc <- lapply(k, FUN=function(k) kmeans(d,k, nstart=3))
  ve <- sapply(cc, FUN=function(x) x$betweenss/x$totss)
  k <- k[farthestPoint(ve)]
  cc <- kmeans(d, k, nstart = 3)$cluster
  gs$cluster <- factor(cc, levels=unique(cc))
  gs
}

jaccard <- function(m){
  nart <- ncol(m)
  jdist <- rep(0, nart*nart)
  dim(jdist) <- c(nart,nart)
  reg.col.sum <- colSums(m)
  reg.aggrement <- t(m) %*% m
  1- reg.aggrement / (reg.col.sum-t(t(reg.aggrement)-reg.col.sum))
}

getWordsFromString <- function(ss){
  for (i in c(" ", "\n", "\r", ";")) ss <- gsub(i, ",", ss,
                                                fixed = T)
  ss <- strsplit(ss, ",", fixed = T)[[1]]
  ss[which(ss != "")]
}

#' enrichmentTests
#'
#' Traditional overlap (Fisher test) between a logical signal (e.g.
#' differentially-expressed features) and the elements of sets.
#'
#' @param signal A named logical vector indicating which features are
#' differentially-expressed
#' @param sets A named list of genesets, or a sparse logical matrix, with sets as
#' columns and features as rows (dimensions must be named).
#' @param alternative Test alternative (defaults to 'greater' to test for
#' over-representation)
#' @param minSize minimum effective size of a set
#' @param maxSize maximum effective size of a set
#'
#' @return a data.frame.
#'
#' @export
enrichmentTests <- function(signal, sets, alternative=c("greater","less","two.sided"), minSize=5, maxSize=1000){
  library(Matrix)
  alternative <- match.arg(alternative)
  if(is.list(sets)){
    feats <- unique(unlist(sets))
    sets <- lapply(sets, FUN=function(x) as.integer(factor(x,feats)))
    m <- sparseMatrix(i=unlist(sets), j=rep(seq_along(sets),lengths(sets)), index1=TRUE)
    row.names(m) <- feats
    colnames(m) <- names(sets)
    sets <- m
  }
  sets <- sets[intersect(row.names(sets), names(signal)),]
  cs <- Matrix::colSums(sets)
  sets <- sets[,which(cs>=minSize & cs<=maxSize)]
  signal <- signal[row.names(sets)]
  sigAndSet <- colSums(sets[which(signal),,drop=FALSE])
  notSigAndSet <- colSums(sets[which(!signal),,drop=FALSE])
  sigAndNotSet <- sum(signal)-sigAndSet
  notNot <- nrow(sets)-sigAndSet-notSigAndSet-sigAndNotSet
  expected <- sum(signal)*(sigAndSet+notSigAndSet)/nrow(sets)
  res <- data.frame( overlap=sigAndSet,
                     enrichment=round(log2((sigAndSet+0.25)/(expected+0.25)),3),
                     pvalue=fisher.test.p(sigAndSet,notSigAndSet,
                                          sigAndNotSet,notNot,
                                          alternative=alternative) )
  res$FDR <- p.adjust(res$pvalue, method="fdr")
  res[order(res$FDR,res$pvalue),]
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


#' topGO wrapper
#'
#' @param genes The vector of genes of interest
#' @param universe The vector of all genes in the assay
#' @param collection The GO collection, either "BP", "MF", or "CC"
#' @param n The number of top terms to report
#' @param db The organism database
#' @param algorithm The algorithm to use.
#'
#' @return A data.frame of enrichment results
runTopGO <- function(genes, universe, collection=c("BP","MF","CC"), n=200,
                     db="org.Mm.eg.db", algorithm=c("classic","elim")){
  library(topGO)
  library(org.Mm.eg.db)
  collection <- match.arg(collection)
  algorithm <- match.arg(algorithm)
  universe <- unique(c(genes, universe))
  signature <- setNames(factor(as.integer(universe %in% genes)), universe)
  o <- new( "topGOdata", ontology=collection, allGenes=signature, nodeSize=5,
            annot=annFUN.org, mapping=db, ID="symbol" )
  res <- runTest(o, algorithm=algorithm, statistic = "Fisher" )
  d <- GenTable( o, Fisher=res, orderBy="Fisher", topNodes=n )
  d$q <- p.adjust(d$Fisher, n=1000)
  d
}
