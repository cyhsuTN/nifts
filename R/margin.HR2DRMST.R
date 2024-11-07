#' Converting HR Margin
#' @description A function to convert HR margin to DRMST margin
#'
#' @import stats
#'
#' @param m1 Median survival time in the active control group.
#' @param shape The shape parameter of Weibull distributions for the event times
#' of two groups. The shape parameters are assumed to be the same.
#' @param tau A value to specify the truncation time point for the RMST calculation.
#' @param theta H1: HR21 < 1/theta (a margin for the HR of new treatment to active control).
#' Superiority if theta = 1; Non-inferiority if 0 < theta < 1.
#' When theta is given, the HR margin is converted to the DRMST margin.
#' i.e., H1: RMST2 - RMST1 > - the DRMST margin.

#' @examples # margin.HR2DRMST(m1=1, shape=1, tau=2.5, theta=0.833)

#' @export
margin.HR2DRMST <- function(m1=1, shape=1, tau=3, theta=0.833) {

  shape1 <- shape

  rate1 <- (log(2)^(1/shape1))/m1

  scale1 <- 1/rate1

  integratef <- function(t, scale, shape) {

    Surv.weibull0 <- function(c) {
      pweibull(c, shape=shape, scale=scale, lower.tail=FALSE)
    }
    integrate(Surv.weibull0, 0, t)$value

  }
  u1 <- integratef(tau, scale1, shape1)

  if(is.null(theta)) stop("Input theta. theta must be between 0 < theta <= 1")

  if(theta<=0 | theta>1) stop("theta must be between 0 < theta <= 1")

  HRm1 <- 1/theta
  scalem <- scale1*theta^(1/shape1)
  um <- integratef(tau, scalem, shape1)
  margin <- - (um - u1)

  return(list(margin=margin, tau=tau))

}
