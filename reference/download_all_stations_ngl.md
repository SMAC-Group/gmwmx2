# Download all stations name and location from the Nevada Geodetic Laboratory

Download all stations name and location from the Nevada Geodetic
Laboratory

## Usage

``` r
download_all_stations_ngl(verbose = FALSE)
```

## Arguments

- verbose:

  A `boolean` that controls the level of detail in the output of the
  `wget` command used to load data. Default is `FALSE`.

## Value

Return a `data.frame` with all stations name, latitude, longitude and
heights.

## Examples

``` r
df_all_stations <- download_all_stations_ngl()
head(df_all_stations)
#>    station_name  latitude longitude    height
#>          <char>     <num>     <num>     <num>
#> 1:         00NA -12.46664 -229.1560 104.85105
#> 2:         01NA -12.47822 -229.0180 105.40857
#> 3:         02NA -12.35592 -229.1183 117.65247
#> 4:         0ABI  68.35434 -341.1836 431.38847
#> 5:         0ABN  65.03368 -338.6671  52.76211
#> 6:         0ABY  58.65891 -343.8204  60.54753
```
