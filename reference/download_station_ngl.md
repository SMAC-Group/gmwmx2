# Download GNSS position time series and steps reference from the Nevada Geodetic Laboratory with IGS14 or IGS20 reference frame.

Download GNSS position time series and steps reference from the Nevada
Geodetic Laboratory with IGS14 or IGS20 reference frame.

## Usage

``` r
download_station_ngl(station_name, verbose = FALSE, reference_frame = "IGS20")
```

## Arguments

- station_name:

  A `string` specifying the station name.

- verbose:

  A `boolean` that controls the level of detail in the output of the
  `wget` command used to load data. Default is `FALSE`.

- reference_frame:

  A `string` with value either "IGS14" or "IGS20" that specify which
  reference frame to use. Default is "IGS20".

## Value

A `list` of class `gnss_ts_ngl` that contains three `data.frame`: The
`data.frame` `df_position` which contains the position time series
extracted from the .tenv3 file available from the Nevada Geodetic
Laboratory, the `data.frame` `df_equipment_software_changes` which
specify the equipment or software changes for that stations and the
`data.frame` `df_earthquakes` that specify the earthquakes associated
with that station.

## Examples

``` r
station_1LSU <- download_station_ngl("1LSU")
attributes(station_1LSU)
#> $names
#> [1] "df_position"                   "df_equipment_software_changes"
#> [3] "df_earthquakes"               
#> 
#> $class
#> [1] "gnss_ts_ngl"
#> 
```
