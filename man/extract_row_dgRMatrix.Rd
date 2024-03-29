% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparse.R
\name{extract_row_dgRMatrix}
\alias{extract_row_dgRMatrix}
\title{Returns a row of a sparse matrix of class dgRMatrix.
Function made by Ben Bolker and Ott Toomet (see https://stackoverflow.com/questions/47997184/)}
\usage{
extract_row_dgRMatrix(m, i = 1)
}
\arguments{
\item{m}{a sparse matrix of class dgRMatrix}

\item{i}{the index of the row to return}
}
\value{
A row (numerical vector) of the sparse matrix
}
\description{
Returns a row of a sparse matrix of class dgRMatrix.
Function made by Ben Bolker and Ott Toomet (see https://stackoverflow.com/questions/47997184/)
}
