## R-CMD-check on local Ubuntu 22.04

==> Rcpp::compileAttributes()

* Updated R/RcppExports.R

==> devtools::check()

══ Documenting ══════════════════════════════════════════════════════════════════════════
ℹ Updating gmwmx2 documentation
ℹ Loading gmwmx2

══ Building ═════════════════════════════════════════════════════════════════════════════
Setting env vars:
• CFLAGS    : -Wall -pedantic -fdiagnostics-color=always
• CXXFLAGS  : -Wall -pedantic -fdiagnostics-color=always
• CXX11FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX14FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX17FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX20FLAGS: -Wall -pedantic -fdiagnostics-color=always
── R CMD build ──────────────────────────────────────────────────────────────────────────
✔  checking for file ‘/home/lionel/github_repo/gmwmx2/DESCRIPTION’ ...
─  preparing ‘gmwmx2’:
✔  checking DESCRIPTION meta-information ...
─  cleaning src
─  installing the package to build vignettes
✔  creating vignettes (59.5s)
─  cleaning src
─  checking for LF line-endings in source and make files and shell scripts (586ms)
─  checking for empty or unneeded directories
   Removed empty directory ‘gmwmx2/.agents’
   Removed empty directory ‘gmwmx2/.codex’
─  building ‘gmwmx2_0.0.5.tar.gz’
   
══ Checking ═════════════════════════════════════════════════════════════════════════════
Setting env vars:
• _R_CHECK_CRAN_INCOMING_USE_ASPELL_           : TRUE
• _R_CHECK_CRAN_INCOMING_REMOTE_               : FALSE
• _R_CHECK_CRAN_INCOMING_                      : FALSE
• _R_CHECK_FORCE_SUGGESTS_                     : FALSE
• _R_CHECK_PACKAGES_USED_IGNORE_UNUSED_IMPORTS_: FALSE
• NOT_CRAN                                     : true
── R CMD check ──────────────────────────────────────────────────────────────────────────
─  using log directory ‘/home/lionel/github_repo/gmwmx2.Rcheck’
─  using R version 4.5.2 (2025-10-31)
─  using platform: x86_64-pc-linux-gnu
─  R was compiled by
       gcc (Ubuntu 11.4.0-1ubuntu1~22.04.2) 11.4.0
       GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04.2) 11.4.0
─  running under: Ubuntu 22.04.5 LTS
─  using session charset: UTF-8
─  using options ‘--no-manual --as-cran’
✔  checking for file ‘gmwmx2/DESCRIPTION’
─  this is package ‘gmwmx2’ version ‘0.0.5’
─  package encoding: UTF-8
✔  checking package namespace information ...
✔  checking package dependencies (1s)
✔  checking if this is a source package ...
✔  checking if there is a namespace
✔  checking for executable files ...
✔  checking for hidden files and directories ...
✔  checking for portable file names ...
✔  checking for sufficient/correct file permissions
─  checking whether package ‘gmwmx2’ can be installed ... [50s/50s] OK (50.1s)
─  used C++ compiler: ‘g++ (Ubuntu 11.4.0-1ubuntu1~22.04.3) 11.4.0’
─  checking installed package size ... INFO
     installed size is 11.4Mb
     sub-directories of 1Mb or more:
       doc      1.7Mb
       libs     5.8Mb
       papers   3.5Mb
✔  checking package directory
N  checking for future file timestamps (1m 0.4s)
   unable to verify current time
✔  checking ‘build’ directory
✔  checking DESCRIPTION meta-information ...
✔  checking top-level files
✔  checking for left-over files
✔  checking index information ...
✔  checking package subdirectories (651ms)
✔  checking code files for non-ASCII characters ...
✔  checking R files for syntax errors ...
✔  checking whether the package can be loaded (954ms)
✔  checking whether the package can be loaded with stated dependencies (798ms)
✔  checking whether the package can be unloaded cleanly (809ms)
✔  checking whether the namespace can be loaded with stated dependencies (759ms)
✔  checking whether the namespace can be unloaded cleanly (951ms)
✔  checking loading without being on the library search path (992ms)
✔  checking dependencies in R code (2.6s)
✔  checking S3 generic/method consistency (925ms)
✔  checking replacement functions (733ms)
✔  checking foreign function calls (892ms)
✔  checking R code for possible problems (8s)
✔  checking Rd files ...
✔  checking Rd metadata ...
✔  checking Rd line widths ...
✔  checking Rd cross-references ...
✔  checking for missing documentation entries (816ms)
✔  checking for code/documentation mismatches (2.3s)
✔  checking Rd \usage sections (1.1s)
✔  checking Rd contents ...
✔  checking for unstated dependencies in examples ...
✔  checking contents of ‘data’ directory ...
✔  checking data for non-ASCII characters ...
✔  checking LazyData
✔  checking data for ASCII and uncompressed saves ...
✔  checking line endings in C/C++/Fortran sources/headers
✔  checking line endings in Makefiles
✔  checking compilation flags in Makevars ...
✔  checking for GNU extensions in Makefiles ...
✔  checking for portable use of $(BLAS_LIBS) and $(LAPACK_LIBS)
✔  checking use of PKG_*FLAGS in Makefiles
✔  checking use of SHLIB_OPENMP_*FLAGS in Makefiles ...
✔  checking pragmas in C/C++ headers and code
✔  checking compilation flags used
✔  checking compiled code ...
✔  checking installed files from ‘inst/doc’ ...
✔  checking files in ‘vignettes’ ...
─  checking examples ... [27s/47s] OK (46.9s)
   Examples with CPU (user + system) or elapsed time > 5s
                                 user system elapsed
   plot.gmwmx2_fit_gnss_ts_ngl 14.488  7.254  24.228
   plot.gnss_ts_ngl             0.793  0.050   7.547
✔  checking for unstated dependencies in vignettes (364ms)
✔  checking package vignettes ...
─  checking re-building of vignette outputs ... [14s/20s] OK (20.5s)
✔  checking for non-standard things in the check directory ...
✔  checking for detritus in the temp directory
   
   See
     ‘/home/lionel/github_repo/gmwmx2.Rcheck/00check.log’
   for details.
   
── R CMD check results ──────────────────────────────────────────────── gmwmx2 0.0.5 ────
Duration: 3m 26.9s

❯ checking for future file timestamps ... NOTE
  unable to verify current time

0 errors ✔ | 0 warnings ✔ | 1 note ✖

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
