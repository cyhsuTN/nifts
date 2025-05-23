% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_power.R
\name{calculate_power}
\alias{calculate_power}
\title{Calculating Powers in Non-Inferiority Trials with Treatment Switching}
\usage{
calculate_power(
  n,
  r = 1,
  m1 = 1,
  m2 = 1.1,
  shape = 1,
  k = 1,
  f1 = NULL,
  m0 = 0.5,
  f2 = 0.5,
  margin = NULL,
  p.s = 0.2,
  r.s = 0.5,
  rho.s = 0.775,
  s.dist = c("unif", "beta", "gamma", "indepExp")[3],
  entry = c("unif", "decreasing", "increasing")[1],
  censoring.prob = c("AC.only"),
  lossfu.dist = c("unif", "exp")[1],
  Ta = 3,
  Te = 5,
  tau = NULL,
  one.sided.alpha = 0.025,
  TXswitch = c("1to2", "2to1")[1],
  af = 1,
  n_simulations = 5000,
  seed = 2024
)
}
\arguments{
\item{n}{Number of participants in the active control group.}

\item{r}{A ratio of participants in the experimental group to
those in the active control group. The default ratio is 1.}

\item{m1}{Median survival time in the active control group.}

\item{m2}{Median survival time in the experimental group.}

\item{shape}{The shape parameter in generalized gamma distributions (flexsurv::pgengamma.orig)
for the event times of two treatment groups.
The same shape parameter is used for the two treatment groups.}

\item{k}{The k parameter in generalized gamma distributions (flexsurv::pgengamma.orig)
for the event times of two treatment groups.
The same k parameter is used for the two treatment groups.}

\item{f1}{Preserved fraction (0 < f1 < 1) for the RMST of the active control group.
If f1 = 0.8, the RMST of the experimental group is not lower than
80% RMST of the active control group.
H1: the RMST of the experimental group > 80% RMST of the active control group.
It is equivalent to
H1: (the RMST of the experimental group - the RMST of the active control group) >
- (1 - f1)(the RMST of the active control group).
margin = (1 - f1)(the RMST of the active control group).}

\item{m0}{Median survival time in the hypothetical placebo group.}

\item{f2}{Preserved fraction (0 < f2 < 1) of the efficacy of the active control group
compared to the hypothetical placebo group.
H1: (the RMST of the experimental group - the RMST of hypothetical placebo group) >
f2 (the RMST of the active control group - the RMST of hypothetical placebo group).
It is equivalent to
H1: (the RMST of the experimental group - the RMST of the active control group) >
- (1 - f2)(the RMST of the active control group - the RMST of hypothetical placebo group).
margin = (1 - f2)(the RMST of the active control group - the RMST of hypothetical placebo group).}

\item{margin}{A margin (>0) for DRMST.}

\item{p.s}{Switching probability of participants who may switch treatment
from the control group to the experimental group after evaluation.
The default is 0.2}

\item{r.s}{A ratio of the mean switching time to the mean survival time in the
active control group. The default is 0.5.}

\item{rho.s}{Correlation between switching time and event time in the active control group.}

\item{s.dist}{An option for the distribution of switching times.
s.dist = "gamma", "beta", "unif", "indepExp", or a constant value.}

\item{entry}{Entry patterns. entry = "increasing", "decreasing", or "unif".}

\item{censoring.prob}{Censoring probability in the active control group
under no treatment switching. The default is c("AC.only").}

\item{lossfu.dist}{The distribution for dropout censoring: uniform ("unif") or
exponential ("exp").}

\item{Ta}{Accrual duration.}

\item{Te}{Trial duration.}

\item{tau}{A value to specify the truncation time point for RMST calculation.
The default is NULL, denoting tau = Te.}

\item{one.sided.alpha}{One-sided significance level. The default is 0.025.}

\item{TXswitch}{Direction of treatment switching (TXswitch="1to2" or "2to1").
The default TXswitch = "1to2".}

\item{af}{A numerical multiplier for the accelerated factor.}

\item{n_simulations}{Number of simulations. The default is 5000.}

\item{seed}{Simulation seed for reproducibility. set.seed(seed).}
}
\value{
\item{Power}{The proportion of rejecting null hypothesis using the adjusted margin}

\item{Power.with.unadj.margin}{The proportion of rejecting null hypothesis
using the unadjusted margin}

\item{E1}{The mean number of events in the active control group}

\item{E2}{The mean number of events in the experimental group}

\item{mRMST1}{The mean DRMST value across simulations in the active control group}

\item{mRMST2}{The mean DRMST value across simulations in the experimental group}

\item{mDRMST}{The mean DRMST value across simulations}

\item{margin.unadj}{unadjusted margin}

\item{margin.adj}{adjusted margin}
}
\description{
A function to calculate powers in non-inferiority trials with treatment switching
by comparing restricted mean survival time between two groups
}
\examples{
# calculate_power(141, r=1, m1=1, m2=1.1, f1=0.8, p.s=0.3, tau=3)
# calculate_power(141, r=1, m1=1, m2=1.1, m0=0.5, f2=0.5, p.s=0.3, tau=3)

# mar <- margin.HR2DRMST(m1=1, shape=1, tau=3, theta=0.833)
# calculate_power(141, r=1, m1=1, m2=1.1, margin=mar$margin, p.s=0.3, tau=3)

# ab <- paramWeibull(m=1, t1=2.5, surv.prob=0.1)
# calculate_power(141, r=1, m1=1, m2=1.1, shape=ab[1],
#                 m0=0.5, f2=0.5, p.s=0.3, tau=3)
}
