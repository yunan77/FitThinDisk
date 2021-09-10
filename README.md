# FitThinDisk
`FitThinDisk` is a Python program to fit the spectral line cubes with a geometrically thin disk model. By giving the spectral line cubes in FITS format, it minimizes the value of &chi;<sup>2</sup> between model and data, and generate the model cube in FITS format in the directory `modelfits`. The input cubes must be three dimensional and the spectral axes should be in velocity instead of frequency.   

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

