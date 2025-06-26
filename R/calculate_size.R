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
#' @param r A ratio of participants in the experimental group to
#' those in the active control group. The default ratio is 1.
#' @param m1 Median survival time in the active control group.
#' @param m2 Median survival time in the experimental group.
#' @param shape The shape parameter in generalized gamma distributions (flexsurv::pgengamma.orig)
#' for the event times of two treatment groups.
#' The same shape parameter is used for the two treatment groups.
#' @param k The k parameter in generalized gamma distributions (flexsurv::pgengamma.orig)
#' for the event times of two treatment groups.
#' The same k parameter is used for the two treatment groups.
#' @param f1 Preserved fraction (0 < f1 < 1) for the RMST of the active control group.
#' If f1 = 0.8, the RMST of the experimental group is not lower than
#' 80% RMST of the active control group.
#' H1: the RMST of the experimental group > 80% RMST of the active control group.
#' It is equivalent to
#' H1: (the RMST of the experimental group - the RMST of the active control group) >
#' - (1 - f1)(the RMST of the active control group).
#' margin = (1 - f1)(the RMST of the active control group).
#' @param m0 Median survival time in the hypothetical placebo group.
#' @param f2 Preserved fraction (0 < f2 < 1) of the efficacy of the active control group
#' compared to the hypothetical placebo group.
#' H1: (the RMST of the experimental group - the RMST of hypothetical placebo group) >
#' f2 (the RMST of the active control group - the RMST of hypothetical placebo group).
#' It is equivalent to
#' H1: (the RMST of the experimental group - the RMST of the active control group) >
#' - (1 - f2)(the RMST of the active control group - the RMST of hypothetical placebo group).
#' margin = (1 - f2)(the RMST of the active control group - the RMST of hypothetical placebo group).
#' @param margin A margin (>0) for DRMST.
#' @param p.s Switching probability of participants who may switch treatment
#' from the control group to the experimental group after evaluation.
#' The default is 0.2
#' @param r.s A ratio of the mean switching time to the mean survival time in the
#' active control group. The default is 0.5.
#' @param rho.s Correlation between switching time and event time in the active control group.
#' @param s.dist An option for the distribution of switching times.
#' s.dist = "gamma", "beta", "unif", "indepExp", or a constant value.
#' @param entry Entry patterns. entry = "increasing", "decreasing", or "unif".
#' @param censoring.prob Censoring probability in the active control group
#' under no treatment switching. The default is c("AC.only").
#' @param lossfu.dist The distribution for dropout censoring: uniform ("unif") or
#' exponential ("exp").
#' @param Ta Accrual duration.
#' @param Te Trial duration.
#' @param tau A value to specify the truncation time point for RMST calculation.
#' The default is NULL, denoting tau = Te.
#' @param one.sided.alpha One-sided significance level. The default is 0.025.
#' @param TXswitch Direction of treatment switching (TXswitch="1to2" or "2to1").
#' The default TXswitch = "1to2".
#' @param af A numerical multiplier for the accelerated factor.
#' @param n_simulations Number of simulations. The default is 5000.
#' @param seed Simulation seed for reproducibility. set.seed(seed).

#' @examples # calculate_size(nL=20, nU=200, r=1, m1=1, m2=1.1, f1=0.8, tau=3)
#' @examples # calculate_size(nL=20, nU=200, r=1, m1=1, m2=1.1, m0=0.5, f2=0.5, tau=3)
#'
#' @examples # ab <- paramWeibull(m=1, t1=2.5, surv.prob=0.1)
#' @examples # calculate_size(nL=20, nU=200, r=1, m1=1, m2=1.1, shape=ab[1],
#' @examples #                m0=0.5, f2=0.5, p.s=0.3, tau=3)

#' @export
calculate_size <- function(nL=10, nU=100, B=10, epwr=0.8,
                           r=1,
                           m1=1,
                           m2=1.1,
                           shape=1,
                           k=1,
                           f1=NULL,
                           m0=0.5,
                           f2=0.5,
                           margin=NULL,
                           p.s=0.2,
                           r.s=0.5,
                           rho.s=0.775,
                           s.dist=c("unif", "beta", "gamma", "indepExp")[3],
                           entry=c("unif", "decreasing", "increasing")[1],
                           censoring.prob=c("AC.only"),
                           lossfu.dist=c("unif", "exp")[1],
                           Ta=3,
                           Te=5,
                           tau=NULL,
                           one.sided.alpha=0.025,
                           TXswitch=c("1to2", "2to1")[1],
                           af=1,
                           n_simulations=1000,
                           seed=2024) {
  powers <- numeric(B)
  E1s <- numeric(B)
  E2s <- numeric(B)

  w <- round((nU - nL) / B)
  for (b in 0:B) {
    n <- nL + b * w
    outs <- calculate_power(n=n, r=r, m1=m1, m2=m2, shape=shape, k=k,
                            f1=f1, m0=m0, f2=f2, margin=margin,
                            p.s=p.s, r.s=r.s, rho.s=rho.s, s.dist=s.dist,
                            entry=entry, censoring.prob=censoring.prob, lossfu.dist=lossfu.dist,
                            Ta=Ta, Te=Te, tau=tau, one.sided.alpha=one.sided.alpha,
                            TXswitch=TXswitch, af=af, n_simulations=n_simulations, seed=seed)
    powers[b+1] <- outs$Power
    E1s[b+1] <- outs$Expected.number.events[1]
    E2s[b+1] <- outs$Expected.number.events[2]
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

  Result = data.frame(req.size = req.size,
                      Power = round(Power,3),
                      E1 = round(E1,1),
                      E2 = round(E2,1))
  nss <- nss[nss$power>=0 & nss$power<=1,]

  list(Result = Result, nss = nss)
}



