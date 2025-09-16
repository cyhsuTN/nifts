nifts: a simulation-based approach to calculate powers and sample sizes required for non-inferiority trials that allow treatment switching in ITT analysis using DRMST
================
Austin Shih & Chih-Yuan Hsu

Austin Shih, Chih‑Yuan Hsu, Yu Shyr (2025). Power and sample size calculation for non-inferiority trials with treatment switching in intention-to-treat analysis comparing RMSTs. *BMC Medical Research Methodology*. 25:157

## Installation

Download nifts_0.4.3.tar.gz and locally install it, or execute the following code:
``` r
library(devtools)
install_github("cyhsuTN/nifts")
```

## Usage
``` r
library(nifts)
```
### Power Calculation
##### -- Option 1 for RMST margin
** Preserved fraction of the RMST of the active control group **
``` r
calculate_power(141, r=1, m1=1, m2=1.1, f1=0.8, p.s=0.3, tau=2.5, Ta=1.5, Te=3)
```
##### -- Option 2 for RMST margin
``` r
calculate_power(141, r=1, m1=1, m2=1.1, m0=0.5, f2=0.5, p.s=0.3, tau=2.5, Ta=1.5, Te=3)
```
##### -- Option 3 for RMST margin
``` r
mar <- margin.HR2DRMST(m1=1, shape=1, tau=2.5, theta=0.833)
calculate_power(141, r=1, m1=1, m2=1.1, margin=mar$margin, p.s=0.3, tau=2.5, Ta=1.5, Te=3)
```
When p.s = 0, there is NO treatment switching, and the study design simplifies to a standard non-inferiority trial.

### Sample Size Calculation
##### Corresponding to Option 1
``` r
pow <- calculate_size(nL=100, nU=200, r=1, m1=1, m2=1.1, f1=0.8, p.s=0.3, tau=2.5, Ta=1.5, Te=3)
pow$Result
plotPower(pow, x.by=10, y.by=0.05)
```
##### Corresponding to Option 2
``` r
pow <- calculate_size(nL=100, nU=200, r=1, m1=1, m2=1.1, m0=0.5, f2=0.5, p.s=0.3, tau=2.5, Ta=1.5, Te=3)
pow$Result
plotPower(pow, x.by=10, y.by=0.05)
```
##### Corresponding to Option 3
``` r
mar <- margin.HR2DRMST(m1=1, shape=1, tau=2.5, theta=0.833)
pow <- calculate_size(nL=200, nU=400, r=1, m1=1, m2=1.1, margin=mar$margin, p.s=0.3,
                      tau=2.5, Ta=1.5, Te=3, n_simulations = 2000)
pow$Result
plotPower(pow, x.by=20, y.by=0.05)
```

### When event times NOT follow exponential distributions
Weibull distributions are assumed for event times of two treatment groups.
The shape and scale parameters are calculated by a given m and survival probability at t1.
``` r
ab <- paramWeibull(m=1, t1=2.5, surv.prob=0.1)
calculate_power(141, r=1, m1=1, m2=1.1, shape=ab[1], k=1, m0=0.5, f2=0.5, p.s=0.3,
                tau=2.5, Ta=1.5, Te=3)
```

``` r
ab <- paramWeibull(m=1, t1=2.5, surv.prob=0.1)
pow <- calculate_size(nL=20, nU=200, r=1, m1=1, m2=1.1, shape=ab[1], k=1, m0=0.5, f2=0.5,
                      p.s=0.3, tau=2.5, Ta=1.5, Te=3)
plotPower(pow, x.by=20, y.by=0.1)
```

Gamma distributions are assumed for event times of two treatment groups.
The shape and scale parameters are calculated by a given m and survival probability at t1.
``` r
ab <- paramGamma(m=1, t1=2.5, surv.prob=0.1)
calculate_power(141, r=1, m1=1, m2=1.1, shape=1, k=ab[1], m0=0.5, f2=0.5, p.s=0.3,
                tau=2.5, Ta=1.5, Te=3)
```

``` r
ab <- paramGamma(m=1, t1=2.5, surv.prob=0.1)
pow <- calculate_size(nL=20, nU=200, r=1, m1=1, m2=1.1, shape=1, k=ab[1], m0=0.5, f2=0.5,
                      p.s=0.3, tau=2.5, Ta=1.5, Te=3)
plotPower(pow, x.by=20, y.by=0.1)
```

