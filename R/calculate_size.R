#' Calculating Sample Sizes in Non-Inferiority Trials with Treatment Switching
#' @description A function to calculate required sample sizes throughout a fitted
#' monotone increasing spline on specified bounds of sample sizes.
#'
#' @import scam
#'
#' @param nL Lower bound of testing interval.
#' @param nU Upper bound of testing interval.
#' @param B Number of intervals splitting upper and lower bound. The default is 10
#' @param epwr Expected power for test. The default is 0.8.
#' @param r The ratio of participants in the new treatment group to
#' those in the active control group. The default ratio is 1.
#' @param m1 Median survival time in the active control group.
#' @param m2 Median survival time in the new treatment group.
#' @param shape The shape parameter of Weibull distributions for the event times
#' of two groups. The shape parameters are assumed to be the same.
#' @param p.s Proportion of patients who may switch from active control to
#' new treatment after evaluation. The default is 0.3
#' @param r.s The ratio of the mean switching time to the mean survival time in the
#' active control group. The default is 0.3.
#' @param rho.s Correlation between switching time and event time in active control.
#' @param s.dist An option for the distribution of switching times.
#' s.dist="gamma", "beta", "unif", "indepExp", or a constant value.
#' @param f1 Preserved fraction (0 < f1 <= 1) for the RMST of active control.
#' If f1 = 0.8, RMST of new treatment is not lower than 80% RMST of active control.
#' That is, H1: RMST of new treatment > 0.8 RMST of active control, equivalent to
#' H1: (RMST of new treatment - RMST of active control) >
#' - (1 - f1)(RMST of active control). margin = (1 - f1)(RMST of active control).
#' @param m0 Median survival time in placebo group.
#' @param f2 Preserved fraction (0 < f2 <= 1) of the efficacy of active control compared to placebo.
#' H1: (RMST of new treatment - RMST of RMST of placebo) >
#' f2 (RMST of active control - RMST of placebo), equivalent to
#' H1: (RMST of new treatment - RMST of active control) >
#' - (1 - f2)(RMST of active control - RMST of placebo).
#' margin = (1 - f2)(RMST of active control - RMST of placebo).
#' @param margin A margin (>0) for DRMST.
#' @param censoring.rate Censoring rate in active control under no treatment switching.
#' The default is c("AC.only").
#' @param Ta Accrual duration.
#' @param Te Trial duration.
#' @param tau A value to specify the truncation time point for the RMST calculation.
#' The default is NULL, denoting tau = Te.
#' @param one.sided.alpha One-sided significance level. The default is 0.025.
#' @param TXswitch direction of treatment switching (TXswitch="1to2" or "2to1").
#' The default TXswitch = "1to2".
#' @param n_simulations Number of simulations. The default is 1000.
#' @param seed set.seed(seed).

#' @examples # calculate_size(nL=20, nU=200, r=1, m1=1, m2=1.1, f1=0.8, tau=2.5)
#' @examples # calculate_size(nL=20, nU=200, r=1, m1=1, m2=1.1, m0=0.5, f2=0.5, tau=2.5)
#'
#' @examples # ab <- paramWeibull(m1=1, t1=2.5, surv.rate=0.1)
#' @examples # calculate_size(nL=20, nU=200, r=1, m1=1, m2=1.5, shape=ab[1],
#' @examples #                m0=0.5, f2=0.5, p.s=0.3, tau=2.5)

#' @export
calculate_size <- function(nL=10, nU=100, B=10, epwr=0.8,
                           r=1,
                           m1=1,
                           m2=1.5,
                           shape=1,
                           f1=NULL,
                           m0=0.5,
                           f2=0.75,
                           margin=NULL,
                           p.s=0.3,
                           r.s=0.3,
                           rho.s=0.3,
                           s.dist="gamma",
                           censoring.rate=c("AC.only"),
                           Ta=1.5,
                           Te=3,
                           tau=NULL,
                           one.sided.alpha=0.025,
                           TXswitch=c("1to2", "2to1")[1],
                           n_simulations=1000,
                           seed=2024) {
  powers <- numeric(B)
  E1s <- numeric(B)
  E2s <- numeric(B)

  w <- round((nU - nL) / B)
  for (k in 0:B) {
    n <- nL + k * w
    outs <- calculate_power(n, r, m1, m2, shape, f1, m0, f2, margin, p.s, r.s, rho.s, s.dist, censoring.rate,
                            Ta, Te, tau, one.sided.alpha, TXswitch, seed, n_simulations)
    powers[k+1] <- outs$Power
    E1s[k+1] <- outs$E1
    E2s[k+1] <- outs$E2
  }

  ss <- data.frame(cbind(n=nL + (0:B) * w, power=powers, E1=E1s, E2=E2s))

  xx <- ss[,1]
  yy <- ss[,2]

  nxx <- seq(min(xx), max(xx), by = 1)

  fit.scam <- scam(yy ~ s(xx, bs="mpi")) #mpi: monotone increasing; mpd: decreasing
  fit.pred <- predict(fit.scam, data.frame(xx = nxx))

  fit.scam1 <- scam(ss$E1 ~ s(xx, bs="mpi"))
  fit.pred1 <- predict(fit.scam1, data.frame(xx = nxx))

  fit.scam2 <- scam(ss$E2 ~ s(xx, bs="mpi"))
  fit.pred2 <- predict(fit.scam2, data.frame(xx = nxx))

  nss <- data.frame(cbind(n = nxx,
                          power = round(fit.pred,4),
                          E1 = round(fit.pred1,2),
                          E2 = round(fit.pred2,2)))


  # To find what sample size in nxx closest to epwr
  idex_close_epwr <- order(abs(fit.pred - epwr))[1]
  # check if the power at req.size is larger than epwr
  idex_close_epwr <- ifelse(nss$power[idex_close_epwr]<epwr, idex_close_epwr+1, idex_close_epwr)
  req.size <- nss$n[idex_close_epwr]
  Power <- nss$power[idex_close_epwr]
  E1 <- nss$E1[idex_close_epwr]
  E2 <- nss$E2[idex_close_epwr]


  if(is.na(req.size)){
    req.size = 'Increase the Upper Bound (nU) to get the required size'
  }

  list(Result = data.frame(req.size = req.size,
                           Power = round(Power,3),
                           E1 = round(E1,1),
                           E2 = round(E2,1)),
       nss = nss)
}



