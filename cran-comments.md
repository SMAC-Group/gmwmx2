## R-CMD-check on local Ubuntu 22.04

==> Rcpp::compileAttributes()

* Updated R/RcppExports.R

==> roxygen2::roxygenize('.', roclets = c('rd', 'collate', 'namespace'))

ℹ Loading gmwmx2
Documentation completed

==> R CMD build gmwmx2

* checking for file ‘gmwmx2/DESCRIPTION’ ... OK
* preparing ‘gmwmx2’:
* checking DESCRIPTION meta-information ... OK
* cleaning src
* installing the package (it is needed to build vignettes)
* creating vignettes ... OK
* cleaning src
* checking for LF line-endings in source and make files and shell scripts
* checking for empty or unneeded directories
Removed empty directory ‘gmwmx2/.agents’
Removed empty directory ‘gmwmx2/.codex’
* building ‘gmwmx2_0.0.6.tar.gz’

==> R CMD check gmwmx2_0.0.6.tar.gz

* using log directory ‘/home/lionel/github_repo/gmwmx2.Rcheck’
* using R version 4.6.1 (2026-06-24)
* using platform: x86_64-pc-linux-gnu
* R was compiled by
    cc (Ubuntu 11.4.0-1ubuntu1~22.04.3) 11.4.0
    GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04.3) 11.4.0
* running under: Ubuntu 22.04.5 LTS
* using session charset: UTF-8
* current time: 2026-10-04 15:06:43 UTC
* checking for file ‘gmwmx2/DESCRIPTION’ ... OK
* this is package ‘gmwmx2’ version ‘0.0.6’
* package encoding: UTF-8
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for executable files ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking for sufficient/correct file permissions ... OK
* checking whether package ‘gmwmx2’ can be installed ... OK
* used C++ compiler: ‘g++ (Ubuntu 11.4.0-1ubuntu1~22.04.3) 11.4.0’
* checking installed package size ... INFO
  installed size is 11.9Mb
  sub-directories of 1Mb or more:
    doc      2.1Mb
    libs     5.9Mb
    papers   3.5Mb
* checking package directory ... OK
* checking ‘build’ directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking code files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* checking whether the package can be loaded ... OK
* checking whether the package can be loaded with stated dependencies ... OK
* checking whether the package can be unloaded cleanly ... OK
* checking whether the namespace can be loaded with stated dependencies ... OK
* checking whether the namespace can be unloaded cleanly ... OK
* checking loading without being on the library search path ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... OK
* checking Rd files ... OK
* checking Rd metadata ... OK
* checking Rd cross-references ... OK
* checking for missing documentation entries ... OK
* checking for code/documentation mismatches ... OK
* checking Rd \usage sections ... OK
* checking Rd contents ... OK
* checking for unstated dependencies in examples ... OK
* checking contents of ‘data’ directory ... OK
* checking data for non-ASCII characters ... OK
* checking LazyData ... OK
* checking data for ASCII and uncompressed saves ... OK
* checking line endings in C/C++/Fortran sources/headers ... OK
* checking line endings in Makefiles ... OK
* checking compilation flags in Makevars ... OK
* checking for GNU extensions in Makefiles ... OK
* checking for portable use of $(BLAS_LIBS) and $(LAPACK_LIBS) ... OK
* checking use of PKG_*FLAGS in Makefiles ... OK
* checking compiled code ... OK
* checking installed files from ‘inst/doc’ ... OK
* checking files in ‘vignettes’ ... OK
* checking examples ... OK
* checking for unstated dependencies in vignettes ... OK
* checking package vignettes ... OK
* checking re-building of vignette outputs ... OK
* checking PDF version of manual ... OK
* DONE
Status: OK


R CMD check succeeded

## R-CMD-check on GitHub actions

All jobs pass on 

- macOS-latest (release)
- windows-latest (release)
- ubuntu-latest (devel)
- ubuntu-latest (release)
- ubuntu-latest (oldrel-1)

see https://github.com/SMAC-Group/gmwmx2/actions/workflows/R-CMD-check.yaml


## Downstream dependencies

There are currently no downstream dependencies for this package.
