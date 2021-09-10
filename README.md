# FitThinDisk
`FitThinDisk` is a Python program to fit the spectral line cubes with a geometrically thin disk model.   

## Download
`git clone https://github.com/yunan77/FitThinDisk.git`

## Requirements
* NumPy
* Astropy
* SciPy

## Description
`FitThinDisk` uses `scipy.optimize.least_squares` to minimizes the value of &chi;<sup>2</sup> between model and data by giving the spectral line cubes. Trust Region Reflective algorithm is selected to perform minimization in order to handle the fitting with bounds.

The input cubes must be three dimensional in FITS format, and the spectral axes should be in velocity instead of frequency. 

There are seven parameters can be fit by this program. In the initial file `params.ini`, users can choose which parameters they want to fit by setting variation and can determine the bounds of the parameters.

## Model Parameters
### Fixed Parameters
* `qr`, power-law index of the rotation velocity as a function of radius, `vr0(r/r0)^qr`, r0 is an arbitrary reference radius and vr0 is the rotation velocity at the reference radius 
* `qd`, power-law index of the disk intensity as a function of radius, `(r/r0)^qd`
* `r0`, reference radius (arcsec)
* `x0`, horizontal offset of the disk center position from the center position of the image (pix)
* `y0`, vertical offset of the disk center position from the center position of the image (pix)
* `pixsize`, size of a pixel in the image (arcsec)
* `dist`, source distance (pc) 
* `bmaj`, major axes width of the beam (arcsec)
* `bmin`, minor axes width of the beam (arcsec)
* `bpa`, position angle of the beam major axis (deg)
### Free Parameters
* `incl`, inclination angel of the disk (deg)(0 for face-on, 90 for edge-on) 
* `pa`, position angel of the major axes of the disk (deg) 
* `vr0sin`, roataion velocity at the reference radius projected onto the line of sight (km/s) 
* `linewidth`, FWHM of Gaussian line profile (km/s)
* `vlsr`, source system velocity (km/s) 
* `ri`, inner radii of the disk (arcsec) 
* `ro`, outer radii of the disk (arcsec)

## Input File
Parameters are read from an input file `params.ini`.

In `[free]`, the parameters are designed by four parts, `initial/fixed value`, `variation`, `min` and `max`.

The variation is used to decide whether the parameter should be fit or not. For `True`, the parameter will be fit and the first value will be taken as initial guess in fitting process; for `False`, the parameter will be fixed at the first value.

The `min` and `max` are designed as bounds of the parameter during fitting process. 

```ini
[file]
imgpath = /home/user/
imglist = ["cube1.fits",
           "cube2.fits",
           "cube3.fits"
          ]
[fixed]
pixsize = 0.02
dist = 3080.
qr = -0.5
qd = -1
r0 = 0.1
bmaj = 0.15
bmin = 0.10
bpa = -50
x0 = 5.72
y0 = -9.14

[free]
vlsr = -15.5     ,True,-30,-10
incl = 25        ,True,0,180
pa = -30         ,True,-90,0
vr0sin = 3       ,False,-inf,inf
linewidth = 5    ,False,0,inf
ri = 0           ,False,0,inf
ro = 5.0         ,False,0,inf
```
## Usage
`python fitdisk.py`

### Output example 
```
Fitting success!
[variables]
vlsr        =  -20  +/- 0.01 (init = -15.5)  (bounds = [-30,-10])
incl        =  50   +/- 0.01 (init = 25)     (bounds = [0,180])
pa          =  -60  +/- 0.01 (init = -30)    (bounds = [-90,0])
vr0sin      =  3             (fixed)
linewidth   =  5             (fixed)
ri          =  0             (fixed)
ro          =  5.0           (fixed)
```
