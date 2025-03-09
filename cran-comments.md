# Local check on Ubuntu 22.04

==> Rcpp::compileAttributes()

* Updated R/RcppExports.R

==> devtools::check()

══ Documenting ════════════════════════════════════════════════════════════════
ℹ Updating gmwmx2 documentation
ℹ Loading gmwmx2

══ Building ═══════════════════════════════════════════════════════════════════
Setting env vars:
• CFLAGS    : -Wall -pedantic -fdiagnostics-color=always
• CXXFLAGS  : -Wall -pedantic -fdiagnostics-color=always
• CXX11FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX14FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX17FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX20FLAGS: -Wall -pedantic -fdiagnostics-color=always
── R CMD build ────────────────────────────────────────────────────────────────
✔  checking for file ‘/home/lionel/github_repo/gmwmx2/DESCRIPTION’ (679ms)
─  preparing ‘gmwmx2’:
✔  checking DESCRIPTION meta-information ...
─  cleaning src
─  installing the package to build vignettes (418ms)
✔  creating vignettes (5m 48.3s)
─  cleaning src
─  checking for LF line-endings in source and make files and shell scripts (1.9s)
─  checking for empty or unneeded directories
─  building ‘gmwmx2_0.0.1.tar.gz’
   
══ Checking ═══════════════════════════════════════════════════════════════════
Setting env vars:
• _R_CHECK_CRAN_INCOMING_USE_ASPELL_           : TRUE
• _R_CHECK_CRAN_INCOMING_REMOTE_               : FALSE
• _R_CHECK_CRAN_INCOMING_                      : FALSE
• _R_CHECK_FORCE_SUGGESTS_                     : FALSE
• _R_CHECK_PACKAGES_USED_IGNORE_UNUSED_IMPORTS_: FALSE
• NOT_CRAN                                     : true
── R CMD check ────────────────────────────────────────────────────────────────
─  using log directory ‘/home/lionel/github_repo/gmwmx2.Rcheck’ (822ms)
─  using R version 4.4.2 (2024-10-31)
─  using platform: x86_64-pc-linux-gnu
─  R was compiled by
       gcc (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
       GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
─  running under: Ubuntu 22.04.5 LTS
─  using session charset: UTF-8
─  using options ‘--no-manual --as-cran’ (410ms)
✔  checking for file ‘gmwmx2/DESCRIPTION’
─  this is package ‘gmwmx2’ version ‘0.0.1’
─  package encoding: UTF-8
✔  checking package namespace information ...
✔  checking package dependencies (3.9s)
✔  checking if this is a source package ...
✔  checking if there is a namespace
✔  checking for executable files (690ms)
✔  checking for hidden files and directories
✔  checking for portable file names
✔  checking for sufficient/correct file permissions
─  checking whether package ‘gmwmx2’ can be installed ... [201s/202s] OK (3m 21.6s)
─  used C++ compiler: ‘g++ (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0’
N  checking installed package size
     installed size is  7.4Mb
     sub-directories of 1Mb or more:
       doc    1.7Mb
       libs   5.4Mb
✔  checking package directory ...
✔  checking for future file timestamps (423ms)
✔  checking ‘build’ directory
✔  checking DESCRIPTION meta-information (656ms)
✔  checking top-level files
✔  checking for left-over files
✔  checking index information (1.1s)
✔  checking package subdirectories (2.4s)
✔  checking code files for non-ASCII characters ...
✔  checking R files for syntax errors ...
✔  checking whether the package can be loaded (7.6s)
✔  checking whether the package can be loaded with stated dependencies (6.4s)
✔  checking whether the package can be unloaded cleanly (5.9s)
✔  checking whether the namespace can be loaded with stated dependencies (5.7s)
✔  checking whether the namespace can be unloaded cleanly (6.5s)
✔  checking loading without being on the library search path (6.6s)
✔  checking dependencies in R code (12.7s)
✔  checking S3 generic/method consistency (6.6s)
✔  checking replacement functions (6.1s)
✔  checking foreign function calls (6.3s)
─  checking R code for possible problems ... [36s/36s] NOTE (36s)
   download_station_ngl: no visible binding for global variable
     ‘step_type_code’
   download_station_ngl: no visible binding for global variable
     ‘date_YYMMDD’
   download_station_ngl: no visible binding for global variable
     ‘type_equipment_change’
   Undefined global functions or variables:
     date_YYMMDD step_type_code type_equipment_change
✔  checking Rd files (825ms)
✔  checking Rd metadata ...
✔  checking Rd line widths ...
✔  checking Rd cross-references (342ms)
✔  checking for missing documentation entries (5.9s)
✔  checking for code/documentation mismatches (17s)
✔  checking Rd \usage sections (6.9s)
✔  checking Rd contents ...
✔  checking for unstated dependencies in examples (389ms)
✔  checking contents of ‘data’ directory ...
✔  checking data for non-ASCII characters (646ms)
✔  checking LazyData
✔  checking data for ASCII and uncompressed saves (594ms)
✔  checking line endings in C/C++/Fortran sources/headers
✔  checking line endings in Makefiles
✔  checking compilation flags in Makevars ...
✔  checking for GNU extensions in Makefiles
✔  checking for portable use of $(BLAS_LIBS) and $(LAPACK_LIBS)
✔  checking use of PKG_*FLAGS in Makefiles
✔  checking use of SHLIB_OPENMP_*FLAGS in Makefiles ...
✔  checking pragmas in C/C++ headers and code
✔  checking compilation flags used
✔  checking compiled code ...
✔  checking installed files from ‘inst/doc’ (881ms)
✔  checking files in ‘vignettes’ (939ms)
─  checking examples ... [38s/39s] OK (40.3s)
   Examples with CPU (user + system) or elapsed time > 5s
                             user system elapsed
   plot.fit_gnss_ts_ngl    11.613  2.190  14.192
   gmwmx2                   6.020  0.606   6.424
   summary.fit_gnss_ts_ngl  4.338  1.047   5.500
✔  checking for unstated dependencies in vignettes (1s)
✔  checking package vignettes ...
─  checking re-building of vignette outputs ... [132s/140s] OK (2m 20.2s)
✔  checking for non-standard things in the check directory
✔  checking for detritus in the temp directory
   
   See
     ‘/home/lionel/github_repo/gmwmx2.Rcheck/00check.log’
   for details.


── R CMD check results ────────────────────────────────────── gmwmx2 0.0.1 ────
Duration: 8m 58.6s

❯ checking installed package size ... NOTE
    installed size is  7.4Mb
    sub-directories of 1Mb or more:
      doc    1.7Mb
      libs   5.4Mb

❯ checking R code for possible problems ... [36s/36s] NOTE
  download_station_ngl: no visible binding for global variable
    ‘step_type_code’
  download_station_ngl: no visible binding for global variable
    ‘date_YYMMDD’
  download_station_ngl: no visible binding for global variable
    ‘type_equipment_change’
  Undefined global functions or variables:
    date_YYMMDD step_type_code type_equipment_change

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

R CMD check succeeded

## R-CMD-check on GitHub actions 

All jobs pass on 

- macOS-latest (release)
- ubuntu-latest (devel)
- ubuntu-latest (oldrel-1)
- ubuntu-latest (release)
- windows-latest (release)

see https://github.com/SMAC-Group/gmwmx2/actions/workflows/R-CMD-check.yaml


# Downstream dependencies
There are currently no downstream dependencies for this package.
