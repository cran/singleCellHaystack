## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 5,
  fig.align = 'center'
)

dat.tsne <- singleCellHaystack::dat.tsne
dat.expression <- singleCellHaystack::dat.expression

## ----example1-----------------------------------------------------------------
library(singleCellHaystack)
set.seed(1234)

# Turn the expression data into detection (gene detected = TRUE, not detected = FALSE)
# let's define detection as having more than 1 read
dat.detection <- dat.expression > 1

# run the main 'haystack' analysis
# inputs are:
# 1) the coordinates of the cells in the input space (here: dat.tsne)
# 2) the detection data (dat.detection)
# 3) the method used ("2D" since we are using a 2D input space here)
res <- haystack(dat.tsne, detection = dat.detection, method = "2D")

# the returned results 'res' is of class 'haystack'
class(res)

## ----example2-----------------------------------------------------------------
# show top 10 DEGs
show_result_haystack(res.haystack = res, n=10)

# alternatively: use a p-value threshold
#show_result_haystack(res.haystack = res, p.value.threshold = 1e-10)

## ----example3-----------------------------------------------------------------
# visualize one of the surprizing genes
plot_gene_haystack(
  dat.tsne,
  expression = dat.expression,
  gene = "gene_497",
  detection = dat.detection,
  high.resolution = TRUE,
  point.size = 2
)

## ----example4-----------------------------------------------------------------
# get the top most significant genes, and cluster them by their distribution pattern in the 2D plot
sorted.table <- show_result_haystack(res.haystack = res, p.value.threshold = 1e-10)
gene.subset <- row.names(sorted.table)

# k-means clustering
km <- kmeans_haystack(dat.tsne, detection=dat.detection, genes=gene.subset, k=5)
km.clusters <- km$cluster

# alternatively: hierarchical clustering
#hc <- hclust_haystack(dat.tsne, detection=dat.detection, genes=gene.subset)
#hc.clusters <- cutree(hc,k = 5)

## ----example5-----------------------------------------------------------------
# visualize cluster distributions
plot_gene_set_haystack(dat.tsne, detection = dat.detection, 
                        genes=names(km.clusters[km.clusters==1]), point.size=2)

