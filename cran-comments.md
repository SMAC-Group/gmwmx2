# Local check on Ubuntu 22.04

==> Rcpp::compileAttributes()

* Updated R/RcppExports.R

==> devtools::check()

══ Documenting ════════════════════════════════════════════
ℹ Updating gmwmx2 documentation
ℹ Loading gmwmx2

══ Building ═══════════════════════════════════════════════
Setting env vars:
• CFLAGS    : -Wall -pedantic -fdiagnostics-color=always
• CXXFLAGS  : -Wall -pedantic -fdiagnostics-color=always
• CXX11FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX14FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX17FLAGS: -Wall -pedantic -fdiagnostics-color=always
• CXX20FLAGS: -Wall -pedantic -fdiagnostics-color=always
── R CMD build ────────────────────────────────────────────
✔  checking for file ‘/home/lionel/github_repo/gmwmx2/DESCRIPTION’ (368ms)
─  preparing ‘gmwmx2’: (388ms)
✔  checking DESCRIPTION meta-information ...
─  cleaning src
─  installing the package to build vignettes
✔  creating vignettes (3m 43.6s)
─  cleaning src
─  checking for LF line-endings in source and make files and shell scripts (1.1s)
─  checking for empty or unneeded directories
─  building ‘gmwmx2_0.0.2.tar.gz’
   
══ Checking ═══════════════════════════════════════════════
Setting env vars:
• _R_CHECK_CRAN_INCOMING_USE_ASPELL_           : TRUE
• _R_CHECK_CRAN_INCOMING_REMOTE_               : FALSE
• _R_CHECK_CRAN_INCOMING_                      : FALSE
• _R_CHECK_FORCE_SUGGESTS_                     : FALSE
• _R_CHECK_PACKAGES_USED_IGNORE_UNUSED_IMPORTS_: FALSE
• NOT_CRAN                                     : true
── R CMD check ────────────────────────────────────────────
─  using log directory ‘/home/lionel/github_repo/gmwmx2.Rcheck’ (531ms)
─  using R version 4.4.3 (2025-02-28)
─  using platform: x86_64-pc-linux-gnu
─  R was compiled by
       gcc (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
       GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0
─  running under: Ubuntu 22.04.5 LTS
─  using session charset: UTF-8
─  using options ‘--no-manual --as-cran’
✔  checking for file ‘gmwmx2/DESCRIPTION’
─  this is package ‘gmwmx2’ version ‘0.0.2’
─  package encoding: UTF-8
✔  checking package namespace information ...
✔  checking package dependencies (2.8s)
✔  checking if this is a source package
✔  checking if there is a namespace
✔  checking for executable files (392ms)
✔  checking for hidden files and directories
✔  checking for portable file names
✔  checking for sufficient/correct file permissions
─  checking whether package ‘gmwmx2’ can be installed ... [135s/136s] OK (2m 15.7s)
─  used C++ compiler: ‘g++ (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0’
N  checking installed package size
     installed size is  7.4Mb
     sub-directories of 1Mb or more:
       doc    1.7Mb
       libs   5.4Mb
✔  checking package directory ...
✔  checking for future file timestamps (433ms)
✔  checking ‘build’ directory
✔  checking DESCRIPTION meta-information (509ms)
✔  checking top-level files
✔  checking for left-over files
✔  checking index information (759ms)
✔  checking package subdirectories (2.1s)
✔  checking code files for non-ASCII characters ...
✔  checking R files for syntax errors ...
✔  checking whether the package can be loaded (5s)
✔  checking whether the package can be loaded with stated dependencies (3.3s)
✔  checking whether the package can be unloaded cleanly (3.8s)
✔  checking whether the namespace can be loaded with stated dependencies (3.8s)
✔  checking whether the namespace can be unloaded cleanly (4.6s)
✔  checking loading without being on the library search path (5s)
✔  checking dependencies in R code (8.3s)
✔  checking S3 generic/method consistency (4.2s)
✔  checking replacement functions (3.9s)
✔  checking foreign function calls (3.7s)
─  checking R code for possible problems ... [25s/25s] OK (25.2s)
✔  checking Rd files (851ms)
✔  checking Rd metadata ...
✔  checking Rd line widths ...
✔  checking Rd cross-references (485ms)
✔  checking for missing documentation entries (5.5s)
✔  checking for code/documentation mismatches (10.6s)
✔  checking Rd \usage sections (3.5s)
✔  checking Rd contents ...
✔  checking for unstated dependencies in examples ...
✔  checking contents of ‘data’ directory ...
✔  checking data for non-ASCII characters (544ms)
✔  checking LazyData
✔  checking data for ASCII and uncompressed saves (344ms)
✔  checking line endings in C/C++/Fortran sources/headers
✔  checking line endings in Makefiles
✔  checking compilation flags in Makevars ...
✔  checking for GNU extensions in Makefiles ...
✔  checking for portable use of $(BLAS_LIBS) and $(LAPACK_LIBS)
✔  checking use of PKG_*FLAGS in Makefiles
✔  checking use of SHLIB_OPENMP_*FLAGS in Makefiles
✔  checking pragmas in C/C++ headers and code ...
✔  checking compilation flags used
✔  checking compiled code ...
✔  checking installed files from ‘inst/doc’ (440ms)
✔  checking files in ‘vignettes’ (430ms)
─  checking examples ... [35s/38s] OK (38.1s)
   Examples with CPU (user + system) or elapsed time > 5s
                             user system elapsed
   plot.fit_gnss_ts_ngl    12.120  2.274  14.801
   gmwmx2                   4.866  0.468   5.354
   summary.fit_gnss_ts_ngl  3.819  1.037   5.142
✔  checking for unstated dependencies in vignettes (808ms)
✔  checking package vignettes ...
─  checking re-building of vignette outputs ... [158s/170s] OK (2m 50.5s)
✔  checking for non-standard things in the check directory
✔  checking for detritus in the temp directory
   
   See
     ‘/home/lionel/github_repo/gmwmx2.Rcheck/00check.log’
   for details.
   
── R CMD check results ────────────────── gmwmx2 0.0.2 ────
Duration: 7m 29.6s

❯ checking installed package size ... NOTE
    installed size is  7.4Mb
    sub-directories of 1Mb or more:
      doc    1.7Mb
      libs   5.4Mb

0 errors ✔ | 0 warnings ✔ | 1 note ✖

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
