#' @title Give Row sums of a Matrix-like Object, Based on a Grouping Variable
#' @description Compute Row sums across columns of a numeric matrix-like object for each level of a grouping variable.
#' @details This function is a wrapper for base function \code{rowsum} and is its "column" version.
#'
#' @param M a matrix, data frame or vector of numeric data. Missing values are allowed. A numeric vector will be treated as a column vector.
#' @param group a vector or factor giving the grouping, with one element per row of M. Default: rownames of M.
#' @param reord if TRUE, then the result will be in order of sort(unique(group)), if FALSE (the default), it will be in the order that groups were encountered.
#' @param na_rm logical (TRUE or FALSE). Should NA (including NaN) values be discarded?
#' @param big is your object big and integer overflow is likely? If TRUE, then M is multiplied by 1.0 to ensure values are of type double (perhaps taking more RAM).
#' @param ... other arguments to be passed to or from methods.
#' @return A matrix-like object containing the sums by group. There will be one row per unique value of group.
#'
#' @author Albert Dorador
#' @export
#' @seealso
#' \code{\link[base]{rowsum}}
#' @examples
#'
#' A <- matrix(1:8, ncol = 4)
#' colnames(A) <- c("A", "B", "A", "B")
#' colsum(A)
#'

colsum <- function(M, group=colnames(M), reord=FALSE, na_rm=FALSE, big=TRUE, ...){
  if (big) {
    M <- M*1.0 #to prevent integer overflow
  }
  return(t(rowsum(t(M), group = group, reorder = reord, na_rm = na_rm, ...)))

}
