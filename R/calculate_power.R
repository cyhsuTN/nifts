#' Calculating Powers in Non-Inferiority Trials with Treatment Switching
#' @description A function to calculate powers in non-inferiority trials with treatment switching
#' by comparing restricted mean survival time between two groups
#'
#' @import survival
#' @import stats
#' @importFrom flexsurv pgengamma.orig rgengamma.orig
#'
#' @param n Number of participants in the active control group.
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

#' @return \item{Power}{The proportion of rejecting null hypothesis using the adjusted margin}
#' @return \item{Power.with.unadj.margin}{The proportion of rejecting null hypothesis
#' using the unadjusted margin}
#' @return \item{E1}{The mean number of events in the active control group}
#' @return \item{E2}{The mean number of events in the experimental group}
#' @return \item{mRMST1}{The mean DRMST value across simulations in the active control group}
#' @return \item{mRMST2}{The mean DRMST value across simulations in the experimental group}
#' @return \item{mDRMST}{The mean DRMST value across simulations}
#' @return \item{margin.unadj}{unadjusted margin}
#' @return \item{margin.adj}{adjusted margin}

#' @examples # calculate_power(141, r=1, m1=1, m2=1.1, f1=0.8, p.s=0.3, tau=3)
#' @examples # calculate_power(141, r=1, m1=1, m2=1.1, m0=0.5, f2=0.5, p.s=0.3, tau=3)
#'
#' @examples # mar <- margin.HR2DRMST(m1=1, shape=1, tau=3, theta=0.833)
#' @examples # calculate_power(141, r=1, m1=1, m2=1.1, margin=mar$margin, p.s=0.3, tau=3)
#'
#' @examples # ab <- paramWeibull(m=1, t1=2.5, surv.prob=0.1)
#' @examples # calculate_power(141, r=1, m1=1, m2=1.1, shape=ab[1],
#' @examples #                 m0=0.5, f2=0.5, p.s=0.3, tau=3)

#' @export
calculate_power <- function(n,
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
                            n_simulations=5000,
                            seed=2024) {

  shape1 <- shape2 <- shape0 <- shape

  scale1 <- calScale(m1, shape=shape1, k=k)
  scale2 <- calScale(m2, shape=shape2, k=k)

  if(!is.null(m0)) {
    scale0 <- calScale(m0, shape=shape0, k=k)
  }

  if(is.null(tau)) tau <- Te


  integratef <- function(t, scale, shape, k) {

    Surv.gengamma.orig0 <- function(c) {
      pgengamma.orig(c, shape=shape, scale=scale, k=k, lower.tail=FALSE)
    }
    integrate(Surv.gengamma.orig0, 0, t)$value

  }
  u1 <- integratef(tau, scale1, shape1, k)


  if(is.null(margin)) {
    if(is.numeric(f1)) {
      M2 <- (1-f1)*u1
    } else {
      u0 <- integratef(tau, scale0, shape0, k)
      M2 <- (1-f2)*(u1-u0)
    }
    margin <- M2
  }


  if(p.s>0) {

    if(is.numeric(s.dist)) s <- s.dist else s <- NULL

    if(is.null(s)) {

      if(TXswitch=="1to2") {
        scale_s <- scale1; shape_s <- shape1
      } else if(TXswitch=="2to1") {
        scale_s <- scale2; shape_s <- shape2
      }

      if(k==1) {
        ## Weibull variance and mean
        V <- (scale_s)^2*(gamma(1+2/shape_s)-gamma(1+1/shape_s)^2)
        mu <- (scale_s)*gamma(1+1/shape_s)
      } else {
        V <- (scale_s)^2*(gamma(k+2/shape_s)/gamma(k)-(gamma(k+1/shape_s)/gamma(k))^2)
        mu <- scale_s * gamma(k+1/shape_s)/gamma(k)
      }

      if(s.dist=="gamma") {
        ab <- sw.gamma(r.s, rho.s, V=V, mu=mu)
        a <- ab[1]; b <- ab[2]
        X.sfun <- function(n) rgamma(n, shape = a, rate = b)
      } else if(s.dist=="beta") {
        ab <- sw.beta(r.s, rho.s, V=V, mu=mu)
        a <- ab[1]; b <- ab[2]
        X.sfun <- function(n) rbeta(n, shape1 = a, shape2 = b)
      } else if(s.dist=="unif") {
        X.sfun <- function(n) runif(n)
      } else if(s.dist=="indepExp") {
        b <- 1/(r.s*mu)
        X.sfun <- function(n) rexp(n, rate = b)
      }
    }

  }


  ## Distribution for entry
  if(Ta>0) {
    a.en <- 2/Ta; b.en <- 2/Ta^2
    if(entry=="increasing") {
      a.star <- 0; b.star <- b.en
    } else if(entry=="decreasing") {
      a.star <- a.en; b.star <- -b.en
    } else if(entry=="unif") {
      a.star <- a.en/2; b.star <- 0
    } else stop("What patern for entry: increasing, decreasing, or unif?")
    fv <- function(v) a.star + b.star*v

    if(b.star==0) {
      inv.Fv <- function(u) u/a.star
    } else {
      inv.Fv <- function(u) (-a.star+sqrt(a.star^2+2*b.star*u))/b.star
    }

  } else inv.Fv <- function(u) return(rep(0, length(u)))


  ## Distribution for dropout censoring
  if(lossfu.dist=="unif") {
    lambda1 <- NULL
    h <- unif.censoring(censoring.prob, scale1, shape1, k, fv, Ta, Te)
    lfu.dist <- function(n) runif(n, min=0, max=h)
  } else if(lossfu.dist=="exp") {
    h <- NULL
    lambda1 <- exp.censoring(censoring.prob, scale1, shape1, k, fv, Ta, Te)
    lfu.dist <- function(n) rexp(n, rate=lambda1)
  } else stop("What distribution for loss follow-up: unif or exp?")



  z <- qnorm(1-one.sided.alpha)

  n1 <- n
  n2 <- n*r

  lbd <- NULL
  mean1 <- NULL
  mean2 <- NULL
  rmst1 <- NULL
  rmst2 <- NULL
  drmst <- NULL

  if (!is.null(seed)){
    set.seed(seed)
  }

  for (i in 1:n_simulations) {

    ### Generate event times from gengamma distributions
    event_times_group1 <- rgengamma.orig(n1, shape=shape1, scale=scale1, k=k)
    event_times_group2 <- rgengamma.orig(n2, shape=shape2, scale=scale2, k=k)

    entry_time_group1 <- inv.Fv(runif(n1, 0, 1))
    entry_time_group2 <- inv.Fv(runif(n2, 0, 1))

    event_times <- c(event_times_group1, event_times_group2)

    entry_times <- c(entry_time_group1, entry_time_group2)

    data <- data.frame(
      Group = factor(rep(0:1, c(n1, n2))),
      EventTime = event_times,
      EntryTime = entry_times
    )

    if (!is.null(h) | !is.null(lambda1)) {
      drop_censoring_time_group1 <- lfu.dist(n1)
      drop_censoring_time_group2 <- lfu.dist(n2)
      drop_censoring_time <- c(drop_censoring_time_group1, drop_censoring_time_group2)
      DropCensoring <- drop_censoring_time
      AdministrativeCensoring <- Te - data$EntryTime
      data$CensoringTime <- pmin(DropCensoring, AdministrativeCensoring)
    } else {
      data$CensoringTime <- Te - data$EntryTime
    }

    data$Y <- pmin(data$EventTime, data$CensoringTime)
    data$Status <- 1*I(data$EventTime < data$CensoringTime)

    if(p.s>0) {

      if(TXswitch=="1to2") {
        n_to_move <- ceiling(p.s * n1)
        indices_to_move <- 1:n_to_move
        A <- af*m2/m1
      } else if(TXswitch=="2to1") {
        n_to_move <- ceiling(p.s * n2)
        indices_to_move <- n1+(1:n_to_move)
        A <- af*m1/m2
      }

      if(n_to_move>0) {
        if (is.null(s)) {
          if(s.dist=="indepExp") {
            sv <- X.sfun(n_to_move)
          } else {
            sv <- data$EventTime[indices_to_move] * X.sfun(n_to_move)
          }
          idx.cross <- which(data$Y[indices_to_move] > sv)
          T.cross <- indices_to_move[idx.cross]
        } else {
          sv <- rep(s, n_to_move)
          idx.cross <- which(data$Y[indices_to_move] > sv)
          T.cross <- indices_to_move[idx.cross]
        }
      } else {
        T.cross <- NULL
      }

      if(length(T.cross)>0) {
        data$EventTime[T.cross] <- sv[idx.cross] + ((data$EventTime[T.cross] - sv[idx.cross]) * A)
        data$Y[T.cross] <- pmin(data$EventTime[T.cross], data$CensoringTime[T.cross])
        data$Status[T.cross] <- 1*I(data$EventTime[T.cross] < data$CensoringTime[T.cross])
      }

    }


    rmean <- tau
    res.km <- summary(survfit(Surv(Y, Status)~Group, data=data), rmean=rmean)
    rmst1est <- res.km$table[1,5]
    rmst1_SE <- res.km$table[1,6]
    rmst2est <- res.km$table[2,5]
    rmst2_SE <- res.km$table[2,6]
    lower_bound <- rmst2est - rmst1est - z*sqrt(rmst1_SE^2+rmst2_SE^2)

    lbd[i] <- lower_bound
    mean1[i] <- sum(data$Status[data$Group == 0])
    mean2[i] <- sum(data$Status[data$Group == 1])
    rmst1[i] <- rmst1est
    rmst2[i] <- rmst2est
    drmst[i] <- rmst2est - rmst1est

  }
  df_status_mean <- data.frame(lbd, mean1, mean2, rmst1, rmst2, drmst)


  if(p.s>0) {
    margin.adj <- mean(df_status_mean$rmst1) - u1 + margin ## adjusted margin
  } else {
    margin.adj <- margin ## unadjusted margin
  }
  power <- round(mean(df_status_mean$lbd > - margin.adj),4)
  power0 <- round(mean(df_status_mean$lbd > - margin),4)

  E1 <- round(mean(df_status_mean$mean1),1)
  E2 <- round(mean(df_status_mean$mean2),1)
  RMST1 <- round(mean(df_status_mean$rmst1),3)
  RMST2 <- round(mean(df_status_mean$rmst2),3)
  DRMST <- round(mean(df_status_mean$drmst),3)

  return(list(Power=power,
              Power.with.unadj.margin=power0,
              Expected.number.events=c(E1=E1, E2=E2),
              RMST.info=c(mRMST1=RMST1, mRMST2=RMST2, mDRMST=DRMST),
              Margin=c(margin.unadj=margin, margin.adj=margin.adj)))

}



## Not export: a function to transform median to scale
calScale <- function(m, shape, k) {

  if(k==1) {
    b <- m/(log(2)^(1/shape))
  } else {
    scalefn <- function(b) {
      scale <- exp(b)
      pgengamma.orig(m, scale=scale, shape=shape, k=k, lower.tail=F)
    }

    fb.surv <- function(b) scalefn(b) - 0.5

    solve.b <- uniroot(fb.surv, interval = c(0, 5), extendInt="yes")
    b <- exp(solve.b$root)
  }

  return(c(scale=b))
}



## Not export: exponential dropout censoring
exp.censoring <- function(censoring.prob, scale1, shape1, k, fv, Ta, Te) {

  if(is.numeric(censoring.prob) & censoring.prob < 1 & censoring.prob > 0) {

    censoring2.fun <- function(scale1, shape1, k, Ta, Te) {
      if(Ta==0) {
        pgengamma.orig(Te, shape=shape1, scale=scale1, k=k, lower.tail=FALSE)
      } else {

        Surv.gengamma.orig2 <- function(v) {
          pgengamma.orig(Te-v, shape=shape1, scale=scale1, k=k, lower.tail=FALSE) * fv(v)
        }
        integrate(Surv.gengamma.orig2, 0, Ta)$value
      }
    }
    AC.rate <- censoring2.fun(scale1, shape1, k, Ta, Te)

    if(censoring.prob <= AC.rate) {
      stop(paste0("The censoring.prob value must be not smaller than ", round(AC.rate,3), " (AC only)." ))
    } else {
      censoring.fun <- function(lambda1, scale1, shape1, k, Ta, Te) {
        if(Ta==0) {

          Surv.gengamma.orig0 <- function(c) {
            pgengamma.orig(c, shape=shape1, scale=scale1, k=k, lower.tail=FALSE) *
              dexp(c, rate=lambda1)
          }

          integrate(Surv.gengamma.orig0, 0, Te)$value +
            pgengamma.orig(Te, shape=shape1, scale=scale1, k=k, lower.tail=FALSE) *
              pexp(Te, rate=lambda1, lower.tail=FALSE)

        } else {

          Surv.gengamma.orig1 <- function(v) {
            sapply(v, function(v1) {
              survfun <- function(c) {
                pgengamma.orig(c, shape=shape1, scale=scale1, k=k, lower.tail=FALSE) *
                  dexp(c, rate=lambda1)
              }
              integrate(survfun, 0, Te-v1)$value
            }) * fv(v)
          }

          Surv.gengamma.orig2 <- function(v) {
            pgengamma.orig(Te-v, shape=shape1, scale=scale1, k=k, lower.tail=FALSE) *
              pexp(Te-v, rate=lambda1, lower.tail=FALSE) * fv(v)
          }

          integrate(Surv.gengamma.orig1, 0, Ta)$value + integrate(Surv.gengamma.orig2, 0, Ta)$value

        }

      }
      fc <- function(lambda1) censoring.fun(lambda1, scale1, shape1, k, Ta, Te) - censoring.prob
      lambda1 <- round(uniroot(fc, interval = c(0.001, 10*Te), extendInt = c("yes"))$root, 3)
    }
  } else {
    lambda1 <- NULL ## administrative censoring only
  }
  return(lambda1)

}


## Not export: uniform dropout censoring
unif.censoring <- function(censoring.prob, scale1, shape1, k, fv, Ta, Te) {

  if(is.numeric(censoring.prob) & censoring.prob < 1 & censoring.prob > 0) {

    censoring2.fun <- function(scale1, shape1, k, Ta, Te) {
      if(Ta==0) {
        pgengamma.orig(Te, shape=shape1, scale=scale1, k=k, lower.tail=FALSE)
      } else {

        Surv.gengamma.orig2 <- function(v) {
          pgengamma.orig(Te-v, shape=shape1, scale=scale1, k=k, lower.tail=FALSE) * fv(v)
        }
        integrate(Surv.gengamma.orig2, 0, Ta)$value
      }
    }
    AC.rate <- censoring2.fun(scale1, shape1, k, Ta, Te)

    if(censoring.prob <= AC.rate) {
      stop(paste0("The censoring.prob value must be not smaller than ", round(AC.rate,3), " (AC only)." ))
    } else {
      censoring.fun <- function(h, scale1, shape1, k, Ta, Te) {
        if(Ta==0) {

          Surv.gengamma.orig0 <- function(c) {
            pgengamma.orig(c, shape=shape1, scale=scale1, k=k, lower.tail=FALSE)
          }

          if(0 <= Te-h) {
            integrate(Surv.gengamma.orig0, 0, h)$value/h
          } else {
            integrate(Surv.gengamma.orig0, 0, Te)$value/h + (1 - Te/h)*Surv.gengamma.orig0(Te)
          }
        } else {

          Surv.gengamma.orig0 <- function(c) {
            pgengamma.orig(c, shape=shape1, scale=scale1, k=k, lower.tail=FALSE)
          }

          Surv.gengamma.orig1 <- function(v) {
            sapply(v, function(v1) {
              survfun <- function(c) {
                pgengamma.orig(c, shape=shape1, scale=scale1, k=k, lower.tail=FALSE)
              }
              integrate(survfun, 0, Te-v1)$value/h
            }) * fv(v)
          }

          Surv.gengamma.orig2 <- function(v) {
            pgengamma.orig(Te-v, shape=shape1, scale=scale1, k=k, lower.tail=FALSE) *
              (1 - (Te-v)/h) * fv(v)
          }

          if(0 < Te-h & Te-h < Ta) {
            integrate(Surv.gengamma.orig0, 0, h)$value/h * integrate(fv, 0, Te-h)$value +
              integrate(Surv.gengamma.orig1, Te-h, Ta)$value +
              integrate(Surv.gengamma.orig2, Te-h, Ta)$value

          } else if(Ta <= Te-h) {
            integrate(Surv.gengamma.orig0, 0, h)$value/h
          } else {
            integrate(Surv.gengamma.orig1, 0, Ta)$value + integrate(Surv.gengamma.orig2, 0, Ta)$value

          }
        }

      }
      fc <- function(h) censoring.fun(h, scale1, shape1, k, Ta, Te) - censoring.prob
      h <- round(uniroot(fc, interval = c(0.001, 10*Te), extendInt = c("yes"))$root, 3)
    }
  } else {
    h <- NULL ## administrative censoring only
  }
  return(h)

}
