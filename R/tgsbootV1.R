#' @title Bootstrap for Stationary Data
#' @description Generate bootstrap samples for stationary data, using a
#' truncated geometric distribution to more accurately determine the value of the
#' \code{p} parameter involved in the algorithm.
#' @details The value of the \code{b} parameter involved in the stationary bootstrap
#' algorithm is determined using the heuristic laid out in Politis & White (2004).
#' Then, the value of the \code{p} parameter is found by numerically solving a
#' polynomial of order N+1 in variable \code{q}, where \eqn{q = 1-p} and N is the
#' length of the data supplied.
#' The previous polynomial is derived using the expectation of a truncated geometric
#' distribution (for the stochastic block length), shown in Olatayo (2014).
#'
#' The general structure of the algorithm is similar to the one laid out in James & Yang (2010).
#'
#' @references
#' Politis, D.N. and Romano, J.P. (1994), 'The stationary bootstrap', Journal of the
#' American Statistical Association 89(428), 1303-1313.
#'
#' Politis, D.N. and White, H. (2004), 'Automatic block-length selection for the
#' dependent bootstrap', Econometric Reviews 23(1), 53-70.
#'
#' Olatayo, T.O. (2014), 'Truncated geometric bootstrap method for timeseries
#' stationary process', Applied Mathematics 5, 2057-2061.
#'
#' James, J. and Yang, L. (2010), 'Stop-losses, maximum drawdown-at-risk and
#' replicating financial time series with the stationary bootstrap',
#' Quantitative Finance 10(1), 1-12.
#'
#' @param tseries a numeric vector or time series giving the original data.
#' @param nb the number of bootstrap series to compute.
#' @param b.info if TRUE, the value of the \code{b} parameter found is returned as well.
#' The default is FALSE.
#' @return If b.info is FALSE, a matrix or time series with nb columns and length(tseries) rows
#' containing the bootstrap data. Each column contains one bootstrap sample.
#' If b.info is TRUE, a list with two fields: one containing the bootstrap data, and
#' another containing the \code{b} value found.
#'
#' @author Albert Dorador
#' @export
#' @importFrom np b.star
#' @examples
#' set.seed(123)
#' x = rnorm(1e4)
#' boot = tgsboot(x)
#' boot = tgsboot(x, b.info = TRUE)
#' boot = tgsboot(x, nb = 2)
#' boot = tgsboot(x, nb = 2, b.info = TRUE)
#'

tgsboot <- function(tseries, nb = 1, b.info = FALSE){
  if (!(nb > 0 && is.numeric(nb) && (nb == floor(nb))))
    stop("nb must be a positive integer.")
  N <- length(tseries)
  b <- b.star(tseries, round = TRUE)[1]
  f <- function(q){
    (N-b)*q^(N+1) + (b-1-N)*q^N + (1+b)*q + 1 - b
  }
  q <- uniroot(f, c(0,1))$root
  p <- 1-q
  names(tseries) <- seq_along(tseries) #create named vector to recover index

  if (nb == 1) {
    x_boot <- numeric(N)

    x_boot.i <- sample(x = tseries, size = 1, replace = TRUE)
    n <- as.integer(names(x_boot.i)) #index in original data
    x_boot[1] <- x_boot.i
    for (i in 2:N) {
      if (n == N) {  # a)
        step4 <- sample(c(0,1), size = 1, replace = TRUE, prob = c(p, (1-p)))
        if (step4 == 0) {
          x_boot.i <- sample(x = tseries, size = 1, replace = TRUE)
          n <- as.integer(names(x_boot.i))
        } else {
          x_boot.i <- tseries[1]
          n <- 1
        }
        x_boot[i] <- x_boot.i
      } else {  # b)
        step4 <- sample(c(0,1), size = 1, replace = TRUE, prob = c(p, (1-p)))
        if (step4 == 0) {
          x_boot.i <- sample(x = tseries, size = 1, replace = TRUE)
          n <- as.integer(names(x_boot.i))
        } else {
          x_boot.i <- tseries[n+1]
          n <- n + 1
        }
        x_boot[i] <- x_boot.i
      }
    }

    if (b.info) {
      return(list(xBoot = x_boot, b_star = b))
    } else {
      return(x_boot)
    }
  } else {
    x_boot <- matrix(nrow = N, ncol = nb)

    for (j in 1:nb) {
      x_boot.i <- sample(x = tseries, size = 1, replace = TRUE)
      n <- as.integer(names(x_boot.i)) #index in original data
      x_boot[1,j] <- x_boot.i
      for (i in 2:N) {
        if (n == N) {  # a)
          step4 <- sample(c(0,1), size = 1, replace = TRUE, prob = c(p, (1-p)))
          if (step4 == 0) {
            x_boot.i <- sample(x = tseries, size = 1, replace = TRUE)
            n <- as.integer(names(x_boot.i))
          } else {
            x_boot.i <- tseries[1]
            n <- 1
          }
          x_boot[i,j] <- x_boot.i
        } else {  # b)
          step4 <- sample(c(0,1), size = 1, replace = TRUE, prob = c(p, (1-p)))
          if (step4 == 0) {
            x_boot.i <- sample(x = tseries, size = 1, replace = TRUE)
            n <- as.integer(names(x_boot.i))
          } else {
            x_boot.i <- tseries[n+1]
            n <- n + 1
          }
          x_boot[i,j] <- x_boot.i
        }
      }
    }

    if (b.info) {
      return(list(xBoot = x_boot, b_star = b))
    } else {
      return(x_boot)
    }
  }
}
