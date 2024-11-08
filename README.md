nifts: a simulation-based approach to calculate powers and sample sizes required for non-inferiority trials that allow treatment switching in ITT analysis using DRMST
================
Austin Shih & Chih-Yuan Hsu

## Installation

Download nifts_0.3.6.tar.gz and locally install it, or execute the following code:
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
``` r
calculate_power(141, r=1, m1=1, m2=1.1, f1=0.8, p.s=0.3, tau=2.5, Ta=1.5, Te=3)
```
##### -- Option 2 for RMST margin
``` r
calculate_power(141, r=1, m1=1, m2=1.1, m0=0.5, f2=0.5, p.s=0.3, tau=2.5, Ta=1.5, Te=3)
```
#### -- Option 3 for RMST margin
``` r
mar <- margin.HR2DRMST(m1=1, shape=1, tau=2.5, theta=0.833)
calculate_power(141, r=1, m1=1, m2=1.1, margin=mar$margin, p.s=0.3, tau=2.5, Ta=1.5, Te=3)
```

### Sample Size Calculation
``` r
mar <- margin.HR2DRMST(m1=1, shape=1, tau=2.5, theta=0.833)
pow <- calculate_size(nL=200, nU=400, r=1, m1=1, m2=1.1, margin=mar$margin, p.s=0.3,
                      tau=2.5, Ta=1.5, Te=3)
plotPower(pow, x.by=40, y.by=0.1)
```

### When event times NOT follow exponential distributions
Weibull distributions are assumed for event times of two treatment groups.
The shape parameter is calculated by a given m1 and survival rate at t1 in control group.
``` r
ab <- paramWeibull(m1=1, t1=2.5, surv.rate=0.1)
calculate_power(141, r=1, m1=1, m2=1.1, shape=ab[1], m0=0.5, f2=0.5, p.s=0.3,
                tau=2.5, Ta=1.5, Te=3)
```

``` r
ab <- paramWeibull(m1=1, t1=2.5, surv.rate=0.1)
pow <- calculate_size(nL=20, nU=200, r=1, m1=1, m2=1.1, shape=ab[1], m0=0.5, f2=0.5,
                      p.s=0.3, tau=2.5, Ta=1.5, Te=3)
plotPower(pow, x.by=20, y.by=0.1)
```
