#' @title Give Row Means of a Matrix-like Object, Based on a Grouping Variable
#' @description Compute Row (weighted) means across columns of a numeric matrix-like object for each level of a grouping variable.
#' @details This function is a wrapper for \pkg{analytics} function \code{rowmean} which allows one to compute the (weighted) mean instead of the sum,
#' while handling integer overflow.
#'
#' Note: although data frames ara allowed, keep in mind that data frames do not allow duplicate row names. Hence if you have a dataframe with more than 1 group, you may want to use the function as.matrix() to convert it to an object of class matrix
#'
#' To compute the mean over all the rows of a matrix (i.e. a single group) use colMeans, which should be even faster.
#' @param M a matrix, data frame or vector of numeric data. Missing values are allowed. A numeric vector will be treated as a column vector.
#' @param group a vector or factor giving the grouping, with one element per row of M. Default: rownames of M.
#' @param w a vector giving the weights that must be applied to each of the stacked blocks of an original object
#' @param reord if TRUE, then the result will be in order of sort(unique(group)), if FALSE (the default), it will be in the order that groups were encountered.
#' @param na_rm logical (TRUE or FALSE). Should NA (including NaN) values be discarded?
#' @param big is your object big and integer overflow is likely? If TRUE, then M is multiplied by 1.0 to ensure values are of type double (perhaps taking more RAM).
#' @param ... other arguments to be passed to or from methods.
#' @return A matrix-like object containing the means by group. There will be one row per unique value of group.
#' If object supplied in fact (explicitly) had just one group, base function
#' \code{colMeans} is called for maximum efficiency and a numeric vector containing
#' the mean of each column is returned.
#'
#' @author Albert Dorador
#' @export
#' @seealso
#' \code{\link[analytics]{rowmean}}
#' \code{\link[base]{rowsum}}
#' @examples
#'
#' A <- matrix(1:8, ncol = 4)
#' colnames(A) <- c("A", "B", "A", "B")
#' colmean(A)
#' colmean(A, w = c(0.2,0.8))
#'

colmean <- function(M, group=colnames(M), w=FALSE, reord=FALSE, na_rm=FALSE, big=TRUE, ...){

return(t(rowmean(t(M), group = group, w = w, reord = reord, na_rm = na_rm, big = big, ...)))

}
