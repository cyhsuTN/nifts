#' Calculating Shape and Scale Parameters for Gamma Distribution
#' @description A function to calculate the shape (k) and scale parameters of a gamma distribution
#' given a median survival and a survival probability at a specified time.
#'
#' @import stats
#'
#' @param m Median survival time.
#' @param t1 Time t1.
#' @param surv.prob Survival probability at time t1.

#' @examples # paramGamma(m=1, t1=2.5, surv.prob=0.1)

#' @export
paramGamma <- function(m, t1, surv.prob) {

  shapefn <- function(b) {
    shape1 <- exp(b)

    scalefn <- function(s) {
      scale1 <- exp(s)
      pgamma(m, scale=scale1, shape=shape1, lower.tail=F)
    }
    fs.surv <- function(s) scalefn(s) - 0.5
    solve.s <- uniroot(fs.surv, interval = c(0, 5), extendInt="yes")
    scale1 <- exp(solve.s$root)

    pgamma(t1, shape=shape1, scale=scale1, lower.tail=F)
  }

  fb.surv <- function(b) shapefn(b) - surv.prob

  solve.b <- uniroot(fb.surv, interval = c(-5, 5), extendInt="yes")
  b <- solve.b$root
  shape1 <- exp(b)

  scalefn <- function(s) {
    scale1 <- exp(s)
    pgamma(m, scale=scale1, shape=shape1, lower.tail=F)
  }
  fs.surv <- function(s) scalefn(s) - 0.5
  solve.s <- uniroot(fs.surv, interval = c(0.01, 5), extendInt="yes")
  scale1 <- exp(solve.s$root)

  return(c(shape_k=shape1, scale=scale1))
}

