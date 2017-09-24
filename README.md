# InterpolateGRIB

### Description:

Script that interpolates data in a GRIB file. 

The script requires two input arguments:

- The input GRIB file.
- The message record.

For example, the script can be run with:

```
python interpolate.py data/gfs.t00z.pgrb2.0p25.f000.grd 7
```
The number seven indicates the record number used to do the interpolations and it refers here to the temperature (isobaric 100 Pa). The program will display the original data, the reprocessed data to a coarser grid and the interpolated data using three interpolation methods.
