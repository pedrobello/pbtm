---
title: "The PBTM R-Package"
output:
  html_document: 
    keep_md: yes
---



## Population-based Threshold Models (PBTM) Calculator - beta version

This R package combines existing PBTMs calculation, parameters output and respective plots. 

## Installation 

### Install development version from Github
devtools::install_github("pedrobello/PBTM")

Initial step:
Set working folder - Add the folder location of your data files


```r
devtools::install_github("pedrobello/PBTM")
```

```
## Downloading GitHub repo pedrobello/PBTM@master
```

```
## Rcpp       (1.0.1 -> 1.0.2) [CRAN]
## tibble     (2.1.1 -> 2.1.3) [CRAN]
## gtable     (0.2.0 -> 0.3.0) [CRAN]
## lazyeval   (0.2.1 -> 0.2.2) [CRAN]
## scales     (0.5.0 -> 1.0.0) [CRAN]
## pillar     (1.3.1 -> 1.4.2) [CRAN]
## colorspace (1.3-2 -> 1.4-1) [CRAN]
```

```
## Installing 7 packages: Rcpp, tibble, gtable, lazyeval, scales, pillar, colorspace
```

```
## 
## The downloaded binary packages are in
## 	/var/folders/cb/v7wys_0d343_wcz4958r9mzm0000gp/T//Rtmpc2AUAi/downloaded_packages
##      checking for file ‘/private/var/folders/cb/v7wys_0d343_wcz4958r9mzm0000gp/T/Rtmpc2AUAi/remotes50e2ca54fcd/pedrobello-PBTM-a58d71c/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/cb/v7wys_0d343_wcz4958r9mzm0000gp/T/Rtmpc2AUAi/remotes50e2ca54fcd/pedrobello-PBTM-a58d71c/DESCRIPTION’
##   ─  preparing ‘PBTM’:
##      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
##   ─  checking for LF line-endings in source and make files and shell scripts
##   ─  checking for empty or unneeded directories
##   ─  building ‘PBTM_0.1.0.tar.gz’
##      
## 
```

## Including Plots

You can also embed plots, for example:

![](PBTM_files/figure-html/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
