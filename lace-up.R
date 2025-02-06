# LACE-UP
#
# Required packages: tidyverse, poLCA, umap, dbscan (if dbscan=TRUE)
# Input: a dataframe of binary data
# Ouput: a low-dimensional embedding with optional dbscan clustering
# Parameters:
#   Ks: vector of candidate values of K of interest
#   nnp: proportion of non-duplicate rows to use for UMAP's nearest neighbor parameter*
#   method: 'umap-learn' (the default) can implement Manhattan distance but requires the python package 'umap-learn'; 'naive' can be used with Euclidean distance
#   dist.metric: distance metric used in calculation of nearest neighbors graph**
#   dbscan: whether to automatically cluster the low-dimensional embedding with dbscan
#   eps.val: eps parameter to use in dbscan***
#   seed: set seed for UMAP reproducibility
#   plot: whether or not to output a plot of the low-dimensional embedding
#   verbose: whether or not to print progress updates
#
# *https://umap-learn.readthedocs.io/en/latest/parameters.html#n-neighbors
# **https://umap-learn.readthedocs.io/en/latest/parameters.html#metric
# ***https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html

library(tidyverse)
library(poLCA)
library(umap)

run.lca <- function(K,df, verbose){
  if(verbose){print(paste0("Running LCA with ", K, " classes"))}
  formula <- as.formula(paste0("cbind(",paste(colnames(df),collapse=","),")~1"))
  lca.model <- poLCA(formula, data=(df+1), nclass=K, na.rm=FALSE, verbose=FALSE)
  return(data.frame(lca.model$posterior))
}

get.pcs <- function(i,lca.models){
  pca <- prcomp(lca.models[[i]], center = TRUE, scale. = TRUE)
  pcs <- as.data.frame(pca$x)
  colnames(pcs) <- paste0(colnames(pcs),"-",i)
  return(pcs)
}

plot.laceup <- function(df, dbscan){
  if(!dbscan){ df$cluster <- "0" }
  gg <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
    geom_point() +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme_bw() +
    theme(legend.position = "none") +
    ggtitle("LACE-UP low-dimensional embedding")
  plot(gg)
  }

run.umap <- function(df, nnp, method, dist.metric, dbscan, eps.val, seed, plot, verbose){
  
  set.seed(seed)
  
  # calculate number of nearest neighbors based on unique rows
  nn <- nrow(df %>% distinct())*nnp
  
  # run UMAP
  if(verbose){print("Running UMAP")}
  umap_fit <- df %>%
    distinct() %>%
    scale() %>%
    umap(n_neighbors = nn,
         metric = dist.metric,
         method = method)
  
  # format output data
  umap_df <- umap_fit$layout %>%
    as.data.frame()
  
  # re-add duplicates
  tmp <- df %>% distinct()
  tmp$UMAP1 = umap_df$V1
  tmp$UMAP2 = umap_df$V2
  out <- df %>%
    left_join(tmp, by = colnames(df)) %>%
    dplyr::select(!starts_with("PC"))
  
  # optional: rub dbscan for automatic clustering
  if(dbscan){
    library(fpc)
    if(verbose){print("Running dbscan")}
    db <- dbscan(out,eps=eps.val)
    out$cluster <- as.character(db$cluster)
  }
  
  # optional: plot output embedding
  if(plot){ plot.laceup(out, dbscan) }
  
  return(out)
}

laceup <- function(df, Ks = 2:20, nnp = 0.1, method="umap-learn", dist.metric="manhattan", dbscan = TRUE, eps.val = 1, seed = 1, plot = TRUE, verbose = TRUE){

  # run LCA on input dataframe for all candidate values of K
  lca.models <- lapply(Ks, run.lca, df, verbose)
  
  # extract and concatenate PCs of posterior probability matrix
  if(verbose){print("Extracting and concatenating principal components")}
  m <- bind_cols(lapply(1:length(lca.models), get.pcs, lca.models))
  
  # run UMAP on PC matrix with optional dbscan clustering
  out <- run.umap(m, nnp, method, dist.metric, dbscan, eps.val, seed, plot, verbose)
  
  return(out)
}






