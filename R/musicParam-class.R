#!/usr/bin/env R

# Author: Sean Maden

#' musicParam-class
#'
#' Applies the MuSiC::music.basic() implementation of the MuSiC deconvolution 
#' algorithm.
#' 
#' @details Main constructor for class \linkS4class{musicParam}.
#' @rdname musicParam-class
#' @seealso 
#' \linkS4class{deconParam}
#' 
#' @examples
#' exampleDataList <- getDeconvolutionExampleData()
#' newParam <- musicParam(exampleDataList[["bulkExpression"]],exampleDataList[["referenceExpression"]],exampleDataList[["cellScaleFactors"]])
#' deconvolution(newParam)
#' 
#' @references 
#' 
#' Wang, Xuran and Jiaxin Fan. MuSiC: Multi-subject single cell deconvolution. 
#' (2022) GitHub, R package version 1.0.0. URL: https://github.com/xuranw/MuSiC.
#' 
#' Wang, X., Park, J., Susztak, K. et al. Bulk tissue cell type deconvolution 
#' with multi-subject single-cell expression reference. Nat Commun 10, 380 
#' (2019). https://doi.org/10.1038/s41467-018-08023-x
#' 
#' @aliases 
#' MuSiCParam-class
#'
setClass("musicParam", contains="referencebasedParam", 
         slots=c(sigma = "matrix", nu = "numeric", 
                 eps = "numeric", iter.max = "numeric"))

#' Make new object of class musicParam
#'
#' Main constructor for class \linkS4class{musicParam}.
#'
#' @param bulkExpression Bulk mixed signals matrix of samples, which can be 
#' matched to single-cell samples.
#' @param referenceExpression Signature matrix of cell type-specific signals. 
#' If not provided, can be computed from a
#' provided ExpressionSet containing single-cell data.
#' @param cellScaleFactors Cell size factor transformations of length equal to 
#' the K cell types to deconvolve.
#' @param sigma Additional argument for algorithm.
#' @param nu Additional argument for algorithm.
#' @param iter.max Additional argument for algorithm.
#' @param eps Additional argument for algorithm.
#' @param returnInfo Whether to return metadata and original method outputs 
#' with predicted proportions.
#'
#' @examples
#' exampleDataList <- getDeconvolutionExampleData()
#' newParam <- musicParam(
#' exampleDataList[["bulkExpression"]],exampleDataList[["referenceExpression"]],
#' exampleDataList[["cellScaleFactors"]]
#' )
#' deconvolution(newParam)
#'
#' @returns New object of class \linkS4class{musicParam}.
#'
#' @details Takes standard inputs for the MuSiC algorithm
#' 
#' @export
musicParam <- function(bulkExpression, referenceExpression, 
                       cellScaleFactors = NULL, sigma = NULL, nu = NULL, 
                       iter.max = NULL, eps = NULL, returnInfo = FALSE) {
  if(is(cellScaleFactors, "NULL")){
    cellScaleFactors <- rep(1, ncol(referenceExpression))}
  if(is(sigma, "NULL")){
    sigma <- matrix(0, ncol = 1, nrow = nrow(referenceExpression))}
  if(is(nu, "NULL")){nu <- 1e-10}
  if(is(iter.max, "NULL")){iter.max <- 1000}
  if(is(eps, "NULL")){eps <- 0}
  new("musicParam", bulkExpression = bulkExpression, 
      referenceExpression = referenceExpression, 
      cellScaleFactors = cellScaleFactors, 
      sigma = sigma, nu = nu, iter.max = iter.max, eps = eps, 
      returnInfo = returnInfo)
}

#' Deconvolution method for class \linkS4class{musicParam}
#' 
#' Main deconvolution method for the \linkS4class{musicParam} to run the 
#' \code{music.basic()} implementation of the MuSiC algorithm.
#' 
#' @param object An object of class \linkS4class{musicParam}.
#' 
#' @examples
#' exampleDataList <- getDeconvolutionExampleData()
#' newParam <- musicParam(exampleDataList[["bulkExpression"]],exampleDataList[["referenceExpression"]],exampleDataList[["cellScaleFactors"]])
#' deconvolution(newParam)
#' 
#' @returns Either a vector of predicted proportions, or a list containing 
#' predictions, metadata, and original outputs.
#'
#' @references 
#' 
#' Wang, Xuran and Jiaxin Fan. MuSiC: Multi-subject single cell deconvolution. 
#' (2022) GitHub, R package version 1.0.0. URL: https://github.com/xuranw/MuSiC.
#' 
#' Wang, X., Park, J., Susztak, K. et al. Bulk tissue cell type deconvolution 
#' with multi-subject single-cell expression reference. Nat Commun 10, 380 
#' (2019). https://doi.org/10.1038/s41467-018-08023-x
#'
#' @export
setMethod("deconvolution", signature(object = "musicParam"), function(object){
  lparam <- callNextMethod()
  # instantiate objects
  nu <- object[["nu"]]
  iter.max <- object[["iter.max"]]
  eps <- object[["eps"]]
  sigma <- object[["sigma"]]
  bulkExpression <- lparam[["bulkExpression"]]
  referenceExpression <- lparam[["referenceExpression"]]
  cellScaleFactors <- lparam[["cellScaleFactors"]]
  sigma <- as.matrix(sigma)
  bulkSamplesIndexVector <- seq(ncol(bulkExpression))
  # run deconvolution algorithm
  result <- lapply(bulkSamplesIndexVector, function(index){
    MuSiC::music.basic(X=referenceExpression, 
                       Y=bulkExpression[,index,drop=F], 
                       S=cellScaleFactors, 
                       Sigma=sigma, nu=nu, iter.max=iter.max, eps=eps)
  })
  # get standardized return list
  names(result) <- colnames(bulkExpression)
  predictions <- lapply(result, function(iter){iter$p.weight})
  returnList <- parseDeconvolutionPredictionsResults(
    predictions, colnames(referenceExpression), colnames(bulkExpression))
  if(object[["returnInfo"]]){
    returnList <- list(
      predictions=predictions,
      result.info = result,
      metadata = parametersList[["metadata"]])}
  return(returnList)
})
