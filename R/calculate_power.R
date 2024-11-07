#' Calculating Powers in Non-Inferiority Trials with Treatment Switching
#' @description A function to calculate powers in non-inferiority trials with treatment switching
#' by comparing restricted mean survival time between two groups
#'
#' @import survival
#' @import stats
#'
#' @param n Number of participants in active control.
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
#' @param n_simulations Number of simulations. The default is 5000.
#' @param seed set.seed(seed).

#' @return \item{power}{The proportion of rejecting null hypothesis}
#' @return \item{E1}{The mean number of events in active control group}
#' @return \item{E2}{The mean number of events in new treatment group}
#' @return \item{mRMST1}{The mean DRMST value across simulations in the active control group}
#' @return \item{mRMST2}{The mean DRMST value across simulations in the new treatment group}
#' @return \item{mDRMST}{The mean DRMST value across simulations}
#' @return \item{margin}{margin}

#' @examples # calculate_power(141, r=1, m1=1, m2=1.1, f1=0.8, p.s=0.3, tau=2.5)
#' @examples # calculate_power(141, r=1, m1=1, m2=1.1, m0=0.5, f2=0.5, p.s=0.3, tau=2.5)
#'
#' @examples # mar <- margin.HR2DRMST(m1=1, shape=1, tau=2.5, theta=0.833)
#' @examples # calculate_power(141, r=1, m1=1, m2=1.1, margin=mar$margin, p.s=0.3, tau=2.5)
#'
#' @examples # ab <- paramWeibull(m1=1, t1=2.5, surv.rate=0.1)
#' @examples # calculate_power(141, r=1, m1=1, m2=1.1, shape=ab[1],
#' @examples #                 m0=0.5, f2=0.5, p.s=0.3, tau=2.5)

#' @export
calculate_power <- function(n,
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
                            n_simulations=5000,
                            seed=2024) {

  shape1 <- shape2 <- shape0 <- shape

  rate1 <- (log(2)^(1/shape1))/m1
  rate2 <- (log(2)^(1/shape2))/m2

  scale1 <- 1/rate1
  scale2 <- 1/rate2

  if(!is.null(m0)) {
    rate0 <- (log(2)^(1/shape0))/m0
    scale0 <- 1/rate0
  }

  if(is.null(tau)) tau <- Te

  if(is.null(margin)) {
    integratef <- function(t, scale, shape) {

      Surv.weibull0 <- function(c) {
        pweibull(c, shape=shape, scale=scale, lower.tail=FALSE)
      }
      integrate(Surv.weibull0, 0, t)$value

    }
    u1 <- integratef(tau, scale1, shape1)

    if(is.numeric(f1)) {
      M2 <- (1-f1)*u1
    } else {
      u0 <- integratef(tau, scale0, shape0)
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

      ## Weibull variance and variance
      V <- (scale_s)^2*(gamma(1+2/shape_s)-gamma(1+1/shape_s)^2)
      mu <- (scale_s)*gamma(1+1/shape_s)

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



  if(is.numeric(censoring.rate)) {

    censoring2.fun <- function(scale1, Ta, Te) {
      if(Ta==0) {
        pweibull(Te, shape=shape1, scale=scale1, lower.tail=FALSE)
      } else {
        Surv.weibull2 <- function(v) {
          pweibull(Te-v, shape=shape1, scale=scale1, lower.tail=FALSE)
        }
        integrate(Surv.weibull2, 0, Ta)$value/Ta
      }
    }
    AC.rate <- censoring2.fun(scale1, Ta, Te)

    if(censoring.rate <= AC.rate) {
      stop(paste0("The censoring.rate value must be not smaller than ", round(AC.rate,3), " (AC only)." ))
    } else {
      censoring.fun <- function(h, scale1, Ta, Te) {
        if(Ta==0) {

          Surv.weibull0 <- function(c) {
            pweibull(c, shape=shape1, scale=scale1, lower.tail=FALSE)
          }

          if(0 <= Te-h) {
            integrate(Surv.weibull0, 0, h)$value/h
          } else {
            integrate(Surv.weibull0, 0, Te)$value/h + (1 - Te/h)*Surv.weibull0(Te)
          }
        } else {

          Surv.weibull0 <- function(c) {
            pweibull(c, shape=shape1, scale=scale1, lower.tail=FALSE)
          }

          Surv.weibull1 <- function(v) {
            sapply(v, function(v1) {
              survfun <- function(c) {
                pweibull(c, shape=shape1, scale=scale1, lower.tail=FALSE)
              }
              integrate(survfun, 0, Te-v1)$value/h
            })
          }

          Surv.weibull2 <- function(v) {
            pweibull(Te-v, shape=shape1, scale=scale1, lower.tail=FALSE) *
              (1 - (Te-v)/h)
          }

          if(0 < Te-h & Te-h < Ta) {
            integrate(Surv.weibull0, 0, h)$value/h * (Te-h)/Ta +
              integrate(Surv.weibull1, Te-h, Ta)$value/Ta +
              integrate(Surv.weibull2, Te-h, Ta)$value/Ta

          } else if(Ta <= Te-h) {
            integrate(Surv.weibull0, 0, h)$value/h
          } else {
            integrate(Surv.weibull1, 0, Ta)$value/Ta + integrate(Surv.weibull2, 0, Ta)$value/Ta

          }
        }

      }
      fc <- function(h) censoring.fun(h, scale1, Ta, Te) - censoring.rate
      h <- round(uniroot(fc, interval = c(0.001, 10*Te), extendInt = c("yes"))$root, 3)
    }
  } else {
    h <- NULL
  }




  z <- qnorm(1-one.sided.alpha)

  n1 <- n
  n2 <- n*r

  reject <- NULL
  mean1 <- NULL
  mean2 <- NULL
  rmst1 <- NULL
  rmst2 <- NULL
  drmst <- NULL

  if (!is.null(seed)){
    set.seed(seed)
  }

  for (i in 1:n_simulations) {

    ### Generate event times from Weibull distributions
    event_times_group1 <- rweibull(n1, shape=shape1, scale=scale1)
    event_times_group2 <- rweibull(n2, shape=shape2, scale=scale2)

    entry_time_group1 <- runif(n1, 0, Ta)
    entry_time_group2 <- runif(n2, 0, Ta)

    event_times <- c(event_times_group1, event_times_group2)
    entry_times <- c(entry_time_group1, entry_time_group2)

    data <- data.frame(
      Group = factor(rep(0:1, c(n1, n2))),
      EventTime = event_times,
      EntryTime = entry_times
    )

    if (!is.null(h)) {
      drop_censoring_time_group1 <- runif(n1, 0, h)
      drop_censoring_time_group2 <- runif(n2, 0, h)
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
        A <- m2/m1
      } else if(TXswitch=="2to1") {
        n_to_move <- ceiling(p.s * n2)
        indices_to_move <- n1+(1:n_to_move)
        A <- m1/m2
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


    rmean <- tau #min(max(data$Y[data$Group==0]), max(data$Y[data$Group==1]), tau)
    res.km <- summary(survfit(Surv(Y, Status)~Group, data=data), rmean=rmean)
    rmst1est <- res.km$table[1,5]
    rmst1_SE <- res.km$table[1,6]
    rmst2est <- res.km$table[2,5]
    rmst2_SE <- res.km$table[2,6]
    lower_bound <- rmst2est - rmst1est - z*sqrt(rmst1_SE^2+rmst2_SE^2)

    reject[i] <- 1*(lower_bound > - margin)
    mean1[i] <- sum(data$Status[data$Group == 0])
    mean2[i] <- sum(data$Status[data$Group == 1])
    rmst1[i] <- rmst1est
    rmst2[i] <- rmst2est
    drmst[i] <- rmst2est - rmst1est

  }
  df_status_mean <- data.frame(reject, mean1, mean2, rmst1, rmst2, drmst)
  E1 <- round(mean(df_status_mean$mean1),1)
  E2 <- round(mean(df_status_mean$mean2),1)
  power <- round(mean(df_status_mean$reject),4)
  RMST1 <- round(mean(df_status_mean$rmst1),3)
  RMST2 <- round(mean(df_status_mean$rmst2),3)
  DRMST <- round(mean(df_status_mean$drmst),3)

  return(list(Power=power, E1=E1, E2=E2,
              mRMST1=RMST1, mRMST2=RMST2, mDRMST=DRMST,
              margin=margin))

}


