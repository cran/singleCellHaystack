#' The main Haystack function
#'
#' @param x a matrix or other object from which coordinates of cells can be extracted.
#' @param assay name of assay data for Seurat method.
#' @param slot name of slot for assay data for Seurat method.
#' @param coord name of coordinates slot for specific methods.
#' @param dims dimensions from coord to use. By default, all.
#' @param cutoff cutoff for detection.
#' @param method choose between highD (default) and 2D haystack.
#' @param expression a matrix with expression data of genes (rows) in cells (columns)
#' @param weights.advanced.Q If NULL naive sampling is used. If a vector is given (of length = no. of cells) sampling is done according to the values in the vector.
#' @param dir.randomization If NULL, no output is made about the random sampling step. If not NULL, files related to the randomizations are printed to this directory.
#' @param scale Logical (default=TRUE) indicating whether input coordinates in x should be scaled to mean 0 and standard deviation 1.
#' @param grid.points An integer specifying the number of centers (gridpoints) to be used for estimating the density distributions of cells. Default is set to 100.
#' @param grid.method The method to decide grid points for estimating the density in the high-dimensional space. Should be "centroid" (default) or "seeding".
#' @param ... further parameters passed down to methods.
#'
#' @return An object of class "haystack"
#' @export
#'
haystack <- function(x, ...) {
  UseMethod("haystack")
}

#' @rdname haystack
#' @export
haystack.matrix <- function(x, expression, weights.advanced.Q = NULL, dir.randomization = NULL, scale = TRUE, grid.points = 100, grid.method = "centroid", ...) {

  haystack_continuous_highD(x,
                            expression = expression,
                            weights.advanced.Q = weights.advanced.Q,
                            dir.randomization = dir.randomization,
                            scale = scale,
                            grid.points = grid.points,
                            grid.method = grid.method, ...)
}

#' @rdname haystack
#' @export
haystack.data.frame <- function(x, expression, weights.advanced.Q = NULL, dir.randomization = NULL, scale = TRUE, grid.points = 100, grid.method = "centroid", ...) {
  haystack(as.matrix(x), expression = expression, weights.advanced.Q = weights.advanced.Q, dir.randomization = dir.randomization, scale = scale, grid.points = grid.points, grid.method = grid.method, ...)
}

#' @rdname haystack
#' @export
haystack.Seurat <- function(x, coord, assay = "RNA", slot = "data", dims = NULL, cutoff = 1, method = NULL, weights.advanced.Q = NULL, ...) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package \"SeuratObject\" needed for this function to work. Please install it.", call. = FALSE)
  }

  if (missing(coord)) stop("Please specify an embedding. One of: ", paste0(SeuratObject::Reductions(x), collapse=", "))

  if (packageVersion("Seurat") >= "5.0.0")
    expression <- SeuratObject::GetAssayData(x, layer = slot, assay = assay)
  else
    expression <- SeuratObject::GetAssayData(x, slot = slot, assay = assay)
  coord <- SeuratObject::Embeddings(x, coord)

  if (! is.null(dims)) {
    coord <- coord[, dims, drop = FALSE]
  }

  haystack(coord, expression = expression, weights.advanced.Q = weights.advanced.Q, ...)
}

#' @rdname haystack
#' @export
haystack.SingleCellExperiment <- function(x, assay = "counts", coord = "TSNE", dims = NULL, cutoff = 1, method = NULL, weights.advanced.Q = NULL, ...) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package \"SummarizedExperiment\" needed for this function to work. Please install it.", call. = FALSE)
  }

  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package \"SingleCellExperiment\" needed for this function to work. Please install it.", call. = FALSE)
  }

  if (missing(coord)) stop("Please specify an embedding. One of: ", paste0(SingleCellExperiment::reducedDimNames(x), collapse=", "))

  expression <- SummarizedExperiment::assay(x, assay)
  coord <- SingleCellExperiment::reducedDim(x, coord)

  if (! is.null(dims)) {
    coord <- coord[, dims, drop = FALSE]
  }

  haystack(coord, expression = expression, weights.advanced.Q = weights.advanced.Q, ...)
}

#' Visualizing the detection/expression of a gene in a 2D plot
#'
#' @param x a matrix or other object from which coordinates of cells can be extracted.
#' @param dim1 column index or name of matrix for x-axis coordinates.
#' @param dim2 column index or name of matrix for y-axis coordinates.
#' @param assay name of assay data for Seurat method.
#' @param slot name of slot for assay data for Seurat method.
#' @param coord name of coordinates slot for specific methods.
#' @param ... further parameters passed to plot_gene_haystack_raw().
#'
#' @export
#'
plot_gene_haystack <- function(x, ...) {
  UseMethod("plot_gene_haystack")
}

#' @rdname plot_gene_haystack
#' @export
plot_gene_haystack.matrix <- function(x, dim1 = 1, dim2 = 2, ...) {
  plot_gene_haystack_raw(x[, dim1], x[, dim2], ...)
}

#' @rdname plot_gene_haystack
#' @export
plot_gene_haystack.data.frame <- function(x, dim1 = 1, dim2 = 2, ...) {
  plot_gene_haystack_raw(x[, dim1], x[, dim2], ...)
}

#' @rdname plot_gene_haystack
#' @export
plot_gene_haystack.SingleCellExperiment <- function(x, dim1 = 1, dim2 = 2, assay = "counts", coord = "TSNE", ...) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package \"SummarizedExperiment\" needed for this function to work. Please install it.", call. = FALSE)
  }

  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package \"SingleCellExperiment\" needed for this function to work. Please install it.", call. = FALSE)
  }

  y <- SummarizedExperiment::assay(x, assay)
  z <- SingleCellExperiment::reducedDim(x, coord)
  plot_gene_haystack_raw(z[, dim1], z[, dim2], expression = y, ...)
}

#' @rdname plot_gene_haystack
#' @export
plot_gene_haystack.Seurat <- function(x, dim1 = 1, dim2 = 2, assay = "RNA", slot = "data", coord = "tsne", ...) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package \"SeuratObject\" needed for this function to work. Please install it.", call. = FALSE)
  }

  y <- SeuratObject::GetAssayData(x, slot = slot, assay = assay)
  z <- SeuratObject::Embeddings(x, coord)
  plot_gene_haystack_raw(z[, dim1], z[, dim2], expression = y, ...)
}

#' Visualizing the detection/expression of a set of genes in a 2D plot
#'
#' @param x a matrix or other object from which coordinates of cells can be extracted.
#' @param dim1 column index or name of matrix for x-axis coordinates.
#' @param dim2 column index or name of matrix for y-axis coordinates.
#' @param assay name of assay data for Seurat method.
#' @param slot name of slot for assay data for Seurat method.
#' @param coord name of coordinates slot for specific methods.
#' @param ... further parameters passed to plot_gene_haystack_raw().
#'
#' @export
#'
plot_gene_set_haystack <- function(x, ...) {
  UseMethod("plot_gene_set_haystack")
}

#' @rdname plot_gene_set_haystack
#' @export
plot_gene_set_haystack.matrix <- function(x, dim1 = 1, dim2 = 2, ...) {
  plot_gene_set_haystack_raw(x[, dim1], x[, dim2], ...)
}

#' @rdname plot_gene_set_haystack
#' @export
plot_gene_set_haystack.data.frame <- function(x, dim1 = 1, dim2 = 2, ...) {
  plot_gene_set_haystack_raw(x[, dim1], x[, dim2], ...)
}

#' @rdname plot_gene_set_haystack
#' @export
plot_gene_set_haystack.SingleCellExperiment <- function(x, dim1 = 1, dim2 = 2, assay = "counts", coord = "TSNE", ...) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package \"SummarizedExperiment\" needed for this function to work. Please install it.", call. = FALSE)
  }

  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package \"SingleCellExperiment\" needed for this function to work. Please install it.", call. = FALSE)
  }

  y <- SummarizedExperiment::assay(x, assay)
  z <- SingleCellExperiment::reducedDim(x, coord)
  plot_gene_set_haystack_raw(z[, dim1], z[, dim2], detection = y > 1, ...)
}

#' @rdname plot_gene_set_haystack
#' @export
plot_gene_set_haystack.Seurat <- function(x, dim1 = 1, dim2 = 2, assay = "RNA", slot = "data", coord = "tsne", ...) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package \"SeuratObject\" needed for this function to work. Please install it.", call. = FALSE)
  }

  y <- SeuratObject::GetAssayData(x, slot = slot, assay = assay)
  z <- SeuratObject::Embeddings(x, coord)
  plot_gene_set_haystack_raw(z[, dim1], z[, dim2], detection = y > 1, ...)
}


#' Function for hierarchical clustering of genes according to their expression distribution in 2D or multi-dimensional space
#'
#' @param x a matrix or other object from which coordinates of cells can be extracted.
#' @param expression expression matrix.
#' @param grid.coordinates coordinates of the grid points.
#' @param hclust.method method used with hclust.
#' @param cor.method method used with cor.
#' @param ... further parameters passed down to methods.
#'
#' @export
#'
hclust_haystack <- function(x, expression, grid.coordinates, hclust.method="ward.D", cor.method="spearman", ...) {
  UseMethod("hclust_haystack")
}

#' @rdname hclust_haystack
#' @export
hclust_haystack.matrix <- function(x, expression, grid.coordinates, hclust.method="ward.D", cor.method="spearman", ...) {
  hclust_haystack_continuous(x, expression, grid.coordinates=grid.coordinates, hclust.method=hclust.method, cor.method=cor.method, ...)
}

#' @rdname hclust_haystack
#' @export
hclust_haystack.data.frame <- function(x, expression, grid.coordinates, hclust.method="ward.D", cor.method="spearman", ...) {
  hclust_haystack_continuous(x, expression, grid.coordinates=grid.coordinates, hclust.method=hclust.method, cor.method=cor.method, ...)
}

#' Function for k-means clustering of genes according to their expression distribution in 2D or multi-dimensional space
#'
#' @param x a matrix or other object from which coordinates of cells can be extracted.
#' @param expression expression matrix.
#' @param grid.coordinates coordinates of the grid points.
#' @param k number of clusters.
#' @param ... further parameters passed down to methods.
#'
#' @export
#'
kmeans_haystack <- function(x, expression, grid.coordinates, k, ...) {
  UseMethod("kmeans_haystack")
}

#' @rdname kmeans_haystack
#' @export
kmeans_haystack.matrix <- function(x, expression, grid.coordinates, k, ...) {
  kmeans_haystack_continuous(x, expression, grid.coordinates=grid.coordinates, k=k, ...)
}

#' @rdname kmeans_haystack
#' @export
kmeans_haystack.data.frame <- function(x, expression, grid.coordinates, k, ...) {
  kmeans_haystack_continuous(x, expression, grid.coordinates=grid.coordinates, k=k, ...)
}
