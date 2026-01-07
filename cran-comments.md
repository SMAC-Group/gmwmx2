# Local check on Ubuntu 22.04
==> Rcpp::compileAttributes()

* Updated R/RcppExports.R

==> devtools::check()

══ Documenting ═══════════════════════════════════════════════════
ℹ Updating gmwmx2 documentation
ℹ Loading gmwmx2

══ Building ══════════════════════════════════════════════════════
Setting env vars:
• CFLAGS    : -Wall -pedantic -fdiagnostics-color=always
• CXXFLAGS  : -Wall -pedantic -fdiagnostics-color=always
• CXX11FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX14FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX17FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX20FLAGS: -Wall -pedantic -fdiagnostics-color=always
── R CMD build ───────────────────────────────────────────────────
✔  checking for file ‘/home/lionel/github_repo/gmwmx2/DESCRIPTION’ (478ms)
─  preparing ‘gmwmx2’: (347ms)
✔  checking DESCRIPTION meta-information ...
─  cleaning src
─  installing the package to build vignettes
✔  creating vignettes (4m 19.1s)
─  cleaning src
─  checking for LF line-endings in source and make files and shell scripts (1.2s)
─  checking for empty or unneeded directories
─  building ‘gmwmx2_0.0.4.tar.gz’
   
══ Checking ══════════════════════════════════════════════════════
Setting env vars:
• _R_CHECK_CRAN_INCOMING_USE_ASPELL_           : TRUE
• _R_CHECK_CRAN_INCOMING_REMOTE_               : FALSE
• _R_CHECK_CRAN_INCOMING_                      : FALSE
• _R_CHECK_FORCE_SUGGESTS_                     : FALSE
• _R_CHECK_PACKAGES_USED_IGNORE_UNUSED_IMPORTS_: FALSE
• NOT_CRAN                                     : true
── R CMD check ───────────────────────────────────────────────────
─  using log directory ‘/home/lionel/github_repo/gmwmx2.Rcheck’ (386ms)
─  using R version 4.5.2 (2025-10-31)
─  using platform: x86_64-pc-linux-gnu
─  R was compiled by
       gcc (Ubuntu 11.4.0-1ubuntu1~22.04.2) 11.4.0
       GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04.2) 11.4.0
─  running under: Ubuntu 22.04.5 LTS
─  using session charset: UTF-8
─  using options ‘--no-manual --as-cran’
✔  checking for file ‘gmwmx2/DESCRIPTION’
─  this is package ‘gmwmx2’ version ‘0.0.4’
─  package encoding: UTF-8
✔  checking package namespace information ...
✔  checking package dependencies (1.7s)
✔  checking if this is a source package ...
✔  checking if there is a namespace
✔  checking for executable files ...
✔  checking for hidden files and directories ...
✔  checking for portable file names
✔  checking for sufficient/correct file permissions
─  checking whether package ‘gmwmx2’ can be installed ... [91s/91s] OK (1m 30.8s)
─  used C++ compiler: ‘g++ (Ubuntu 11.4.0-1ubuntu1~22.04.2) 11.4.0’
─  checking installed package size ... INFO
     installed size is  7.5Mb
     sub-directories of 1Mb or more:
       doc    1.7Mb
       libs   5.5Mb
✔  checking package directory
N  checking for future file timestamps (32.4s)
   unable to verify current time
✔  checking ‘build’ directory
✔  checking DESCRIPTION meta-information (482ms)
✔  checking top-level files
✔  checking for left-over files
✔  checking index information (534ms)
✔  checking package subdirectories (1.5s)
✔  checking code files for non-ASCII characters ...
✔  checking R files for syntax errors ...
✔  checking whether the package can be loaded (3.8s)
✔  checking whether the package can be loaded with stated dependencies (3s)
✔  checking whether the package can be unloaded cleanly (3.3s)
✔  checking whether the namespace can be loaded with stated dependencies (3.2s)
✔  checking whether the namespace can be unloaded cleanly (3.4s)
✔  checking loading without being on the library search path (3.5s)
✔  checking dependencies in R code (6.9s)
✔  checking S3 generic/method consistency (3.3s)
✔  checking replacement functions (3.4s)
✔  checking foreign function calls (4s)
─  checking R code for possible problems ... [23s/23s] OK (22.8s)
✔  checking Rd files (510ms)
✔  checking Rd metadata ...
✔  checking Rd line widths ...
✔  checking Rd cross-references ...
✔  checking for missing documentation entries (5.5s)
✔  checking for code/documentation mismatches (13.4s)
✔  checking Rd \usage sections (4.5s)
✔  checking Rd contents ...
✔  checking for unstated dependencies in examples ...
✔  checking contents of ‘data’ directory ...
✔  checking data for non-ASCII characters (471ms)
✔  checking LazyData
✔  checking data for ASCII and uncompressed saves (465ms)
✔  checking line endings in C/C++/Fortran sources/headers
✔  checking line endings in Makefiles
✔  checking compilation flags in Makevars ...
✔  checking for GNU extensions in Makefiles
✔  checking for portable use of $(BLAS_LIBS) and $(LAPACK_LIBS)
✔  checking use of PKG_*FLAGS in Makefiles
✔  checking use of SHLIB_OPENMP_*FLAGS in Makefiles ...
✔  checking pragmas in C/C++ headers and code ...
✔  checking compilation flags used ...
✔  checking compiled code ...
✔  checking installed files from ‘inst/doc’ (628ms)
✔  checking files in ‘vignettes’ (629ms)
─  checking examples ... [27s/35s] OK (35.7s)
   Examples with CPU (user + system) or elapsed time > 5s
                            user system elapsed
   plot.fit_gnss_ts_ngl    6.821  0.993   9.151
   summary.fit_gnss_ts_ngl 4.558  1.136   6.730
   gmwmx2                  3.713  0.336   5.077
✔  checking for unstated dependencies in vignettes (845ms)
✔  checking package vignettes ...
─  checking re-building of vignette outputs ... [160s/181s] OK (3m 0.7s)
✔  checking for non-standard things in the check directory
✔  checking for detritus in the temp directory
   
   See
     ‘/home/lionel/github_repo/gmwmx2.Rcheck/00check.log’
   for details.
   
── R CMD check results ───────────────────────── gmwmx2 0.0.4 ────
Duration: 7m 15.8s

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


# Downstream dependencies
There are currently no downstream dependencies for this package.
