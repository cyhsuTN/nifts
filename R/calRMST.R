#' Calculate RMST
#' @description A function to calculate RMST
#'
#' @import stats
#' @importFrom flexsurv pgengamma.orig rgengamma.orig
#'
#' @param m Median survival time.
#' @param shape The shape parameter in generalized gamma distributions
#' (flexsurv::pgengamma.orig) for event time.
#' @param k The k parameter in generalized gamma distributions
#' (flexsurv::pgengamma.orig) for event time.
#' @param tau A value to specify the truncation time point for the RMST calculation.

#' @examples # calRMST(m=1, shape=1, k=1, tau=3)

#' @export
calRMST <- function(m=1, shape=1, k=1, tau=NULL) {

  scale1 <- calScale(m, shape=shape, k=k)

  integratef <- function(t, scale, shape, k) {

    Surv.gengamma.orig0 <- function(c) {
      pgengamma.orig(c, shape=shape, scale=scale, k=k, lower.tail=FALSE)
    }
    integrate(Surv.gengamma.orig0, 0, t)$value

  }
  u1 <- integratef(tau, scale1, shape, k)

  return(u1)

}
