#' @title Weighted Biplot
#' @description \code{Wbiplot} produces a biplot with any weight distribution between Row and Column markers.
#' This way the full spectrum from perfect row resolution (Row-metric preserving biplot)
#' to perfect column resolution (Column-metric preserving biplot) is available.
#' @details This function makes use of function \code{Matpow} from package \pkg{powerplus} to
#' be able to raise any valid matrix (see \code{Matpow} documentation) to any real power between 0 and 1 included.
#' @param df a dataframe with numeric values only
#' @param numer1 numerator of first exponent (can be a decimal)
#' @param denom1 denominator of first exponent (default: 1)
#' @param numer2 numerator of second exponent (can be a decimal)
#' @param denom2 denominator of second exponent (default: 1)
#' @param cx graphical magnification factor (default: 0.5)
#' @return A biplot of a dataframe with the specified weights.
#' Weights can either be supplied as two fractions, or as two decimal numbers.
#' @author Albert Dorador
#' @export
#' @import powerplus
#' @seealso \code{\link[powerplus]{Matpow}}
#' @examples
#'
#' require(graphics)
#'
#' # Exemple 1: Row metric preserving
#' Wbiplot(USArrests, numer1 = 1, numer2 = 0, cx = 0.6)
#'
#' # Exemple 2: Column metric preserving
#' Wbiplot(USArrests, numer1 = 0, numer2 = 1, cx = 0.6)
#'
#' # Comparison with function \code{biplot} from package \pkg{stats}
#' biplot(princomp(USArrests), cex = 0.6)
#'
#' # Example 3: Custom, 50-50
#' Wbiplot(USArrests, numer1 = 0.5, numer2 = 0.5)
#'
#' # Example 4: Custom, 20-80
#' Wbiplot(USArrests, numer1 = 0.2, numer2 = 0.8)
#'

Wbiplot <- function(df, numer1, denom1 = 1, numer2, denom2 = 1, cx = 0.5){
  stopifnot((numer1/denom1+numer2/denom2)==1 && (numer1/denom1)<=1 && (numer2/denom2)<=1 && (numer1/denom1)>=0 && (numer2/denom2)>=0)
  for(i in 1:ncol(df)){
    df[, i] <- as.numeric(df[, i])
  }
  X <- scale(df, scale = T)
  Dec <- svd(X)
  U <- Dec$u
  D <- diag(Dec$d)
  V <- Dec$v
  diag(svd(D)$d)
  D1 <- Matpow(D, numer1, denom1)
  D2 <- Matpow(D, numer2, denom2)
  Fps <- U %*% D1
  Gsp <- V %*% D2
  biplot(Fps, Gsp, xlab = "Comp. 1", ylab = "Comp. 2", xlabs = row.names(df), ylabs = colnames(df),
         main = paste0("Custom biplot: ", numer1/denom1, "-", numer2/denom2),
         cex = cx)
}
