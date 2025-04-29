#' Calculating Shape and Scale Parameters for Weibull Distribution
#' @description A function to calculate the shape and scale parameters of a Weibull distribution
#' given a median survival and a survival probability at a specified time.
#'
#' @import stats
#'
#' @param m Median survival time.
#' @param t1 Time t1.
#' @param surv.prob Survival probability at time t1.

#' @examples # paramWeibull(m=1, t1=2.5, surv.prob=0.1)

#' @export
paramWeibull <- function(m, t1, surv.prob) {

  shapefn <- function(b) {
    shape1 <- exp(b)
    scale1 <- m/(log(2)^(1/shape1))
    pweibull(t1, shape=shape1, scale=scale1, lower.tail=F)
  }

  fb.surv <- function(b) shapefn(b) - surv.prob

  solve.b <- uniroot(fb.surv, interval = c(-5, 5), extendInt="yes")
  b <- solve.b$root
  shape1 <- exp(b)
  scale1 <- m/(log(2)^(1/shape1))

  return(c(shape=shape1, scale=scale1))
}


