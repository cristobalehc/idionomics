## *Idionomics*

#### Making idionomic analyses accessible and easy to use. 

This is the first version of the "idionomics" package. Currently, it has the following wrapper functions with robot names -development version only-: 

* imputatron_2000: Impute data via copy-means.
* i_standarbot_300: Person-mean center your time-series.
* IARIMAXoid_Pro: Run the I-ARIMAX algorithm.
  * 21-05-2024: Added argument to perform a comparison with an HLM growth curve model with linear time. 

#### Usage:

You can install the development version of idionomics.

``` r
# install.packages("devtools")
devtools::install_github("cristobalehc/idionomics")
```

#### A word of caution:

This is a developer version, so it can and will contain errors. The first commit was on 02-05-2024, so it is a new and roughly untested project.  
