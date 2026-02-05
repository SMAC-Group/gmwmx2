# Download estimated velocities provided by the Nevada Geodetic Laboratory for all stations.

Download estimated velocities provided by the Nevada Geodetic Laboratory
for all stations.

## Usage

``` r
download_estimated_velocities_ngl(verbose = FALSE)
```

## Arguments

- verbose:

  A `boolean` that controls the level of detail in the output of the
  `wget` command used to load data. Default is `FALSE`.

## Value

Return a `data.frame` with all stations name, information about the time
series for each station, estimated velocities and estimated standard
deviation of the estimated velocities.

## Examples

``` r
df_estimated_velocities <- download_estimated_velocities_ngl()
head(df_estimated_velocities)
#>    station_name midas_version_label time_series_duration_year
#>          <char>              <char>                     <num>
#> 1:         00NA              MIDAS5                   10.4969
#> 2:         01NA              MIDAS5                   11.4743
#> 3:         02NA              MIDAS5                    8.2738
#> 4:         0ABI              MIDAS5                   15.2992
#> 5:         0ABN              MIDAS5                    1.4840
#> 6:         0ABY              MIDAS5                    9.4703
#>    east_velocity_m_yr north_velocity_m_yr up_velocity_m_yr
#>                 <num>               <num>            <num>
#> 1:           0.036213            0.058799        -0.001111
#> 2:           0.035826            0.059595        -0.001008
#> 3:           0.036171            0.059940        -0.000560
#> 4:           0.015129            0.014934         0.005377
#> 5:           0.015129            0.017643         0.005281
#> 6:           0.018483            0.013913         0.005170
#>    east_velocity_unc_m_yr north_velocity_unc_m_yr up_velocity_unc_m_yr
#>                     <num>                   <num>                <num>
#> 1:               0.000233                0.000262             0.000919
#> 2:               0.000250                0.000256             0.000929
#> 3:               0.000254                0.000279             0.001063
#> 4:               0.000164                0.000165             0.000624
#> 5:               0.000772                0.000883             0.003350
#> 6:               0.000206                0.000186             0.000738
#>     latitude longitude    height
#>        <num>     <num>     <num>
#> 1: -12.46664 -229.1560 104.84337
#> 2: -12.47822 -229.0180 105.37638
#> 3: -12.35592 -229.1183 117.65070
#> 4:  68.35434 -341.1836 431.36226
#> 5:  65.03368 -338.6671  52.76520
#> 6:  58.65891 -343.8204  60.52854
```
