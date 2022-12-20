## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = 'center'
)

library(ggplot2)
theme_set(cowplot::theme_cowplot())

## ----example1-----------------------------------------------------------------
library(singleCellHaystack)
set.seed(1234)

# run the main 'haystack' analysis
# inputs are:
# 1) the coordinates of the cells in the input space (here: dat.tsne)
# 2) the expression data (dat.expression)
res <- haystack(dat.tsne, dat.expression)

# the returned results 'res' is of class 'haystack'
class(res)

## ----example2-----------------------------------------------------------------
# show top 10 DEGs
show_result_haystack(res.haystack = res, n=10)

# alternatively: use a p-value threshold
#show_result_haystack(res.haystack = res, p.value.threshold = 1e-10)

## -----------------------------------------------------------------------------
d <- cbind(dat.tsne, t(dat.expression))
d[1:4, 1:4]

## ----fig.width=6, fig.height=4------------------------------------------------
library(ggplot2)

ggplot(d, aes(tSNE1, tSNE2, color=gene_497)) +
  geom_point() +
  scale_color_distiller(palette="Spectral")

## ----example4-----------------------------------------------------------------
# get the top most significant genes, and cluster them by their distribution pattern in the 2D plot
sorted.table <- show_result_haystack(res.haystack = res, p.value.threshold = 1e-10)
gene.subset <- row.names(sorted.table)

# k-means clustering
#km <- kmeans_haystack(dat.tsne, dat.expression[gene.subset, ], grid.coordinates=res$info$grid.coordinates, k=5)
#km.clusters <- km$cluster

# alternatively: hierarchical clustering
hm <- hclust_haystack(dat.tsne, dat.expression[gene.subset, ], grid.coordinates=res$info$grid.coordinates)

## ----fig.width=6, fig.height=8------------------------------------------------
ComplexHeatmap::Heatmap(dat.expression[gene.subset, ], show_column_names=FALSE, cluster_rows=hm, name="expression")

## -----------------------------------------------------------------------------
hm.clusters <- cutree(hm, k=4)
table(hm.clusters)

## -----------------------------------------------------------------------------
for (cluster in unique(hm.clusters)) {
  d[[paste0("cluster_", cluster)]] <- colMeans(dat.expression[names(which(hm.clusters == cluster)), ])
}

## ----fig.width=8, fig.height=6------------------------------------------------
lapply(c("cluster_1", "cluster_2", "cluster_3", "cluster_4"), function(cluster) {
  ggplot(d, aes(tSNE1, tSNE2, color=.data[[cluster]])) +
  geom_point() +
  scale_color_distiller(palette="Spectral")
}) |> patchwork::wrap_plots()

