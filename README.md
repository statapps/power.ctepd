# power.ctepd

## Power and sample size for clinical trials and epidemiology

"power.ctepd" is a R package for power and sample size calculations.
Please use the following steps to install 'power.ctepd' package:

1. First, you need to install the 'devtools' package. You can skip this step if you have '
devtools' inst
alled in your R. Invoke R and then type

  install.packages("devtools")

2. Load the devtools package.

  library(devtools)

3. Install "power.ctepd" package with R command

  install_github("statapps/ctepd")


Version history,

1. version 0.05 is the initial version, with support for one sample survival distribution (2016).
2. version 0.06: add sample size for Simon's pick winner design(2016).
3. version 0.07: add sample size for case control SNP analysis (2017).
4. version 0.08: add sample size for false discovery rate (FDR) with t-tests (2017).
