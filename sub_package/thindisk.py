#!/usr/local/bin/python3
import numpy as np
import configparser
import math
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from scipy.signal import convolve as scipy_convolve
from .getparams import ReadCube,Parameters

#################################################
## please setup input parameters in params.ini ##
#################################################

## read fixed parameters ##
params = Parameters('params.ini')

sampsize = params.fixed['samplesize']
pixsize = params.fixed['pixsize']
dist = params.fixed['dist']
qr = params.fixed['qr']
qd = params.fixed['qd']
r0 = params.fixed['r0']
bmaj = params.fixed['bmaj']
bmin = params.fixed['bmin']
bpa = params.fixed['bpa']
x0 = params.fixed['x0']
y0 = params.fixed['y0']
samfrac = int(pixsize/sampsize)
xleft = int(params.fixed['imgrange'][0])
xright = int(params.fixed['imgrange'][1])+1
ydown = int(params.fixed['imgrange'][2])
yup = int(params.fixed['imgrange'][3])+1
xlre = xleft*samfrac
xrre = xright*samfrac
ydre = ydown*samfrac
yure = yup*samfrac

####################
## thindisk model ##
####################
def thindisk(vlsr,incl,pa,vr0,linewidth,ri,ro,imgpar):
    # input image parameters
    data,npix,nchan,velwidth,velstart = imgpar
    
    # ra,dec,vel grids on the plane of the sky 
    samppix = int(npix*pixsize/sampsize)
    ra_offset = ((np.arange(samppix)-(samppix-1)/2)*sampsize - x0*pixsize)*-1
    dec_offset = (np.arange(samppix)-(samppix-1)/2)*sampsize - y0*pixsize
    ra_grid, dec_grid = np.meshgrid(ra_offset,dec_offset)
    vel = (np.arange(nchan)*velwidth + velstart)/1000 # km/s

    # x,y grids on the plane of the disk 
    pa = pa*math.pi/180.
    x_grid = +ra_grid*np.cos(pa) - dec_grid*np.sin(pa) 
    y_grid = +ra_grid*np.sin(pa) + dec_grid*np.cos(pa)
    costheta = np.zeros((samppix, samppix))
    r = np.zeros((samppix, samppix))
    mask = y_grid != 0
    incl = incl*math.pi/180.
    r = np.sqrt(x_grid**2 + (y_grid/np.cos(incl))**2) # acrsec
    costheta[mask] = x_grid[mask]/r[mask]
    r = r*dist # AU

    # calculate line peak intensity 
    peakint = np.zeros((samppix, samppix))
    mask = r != 0
    peakint[mask] = 10*pow(r[mask]/(r0*dist), qd)

    # calculate disk velocity 
    vtheta = np.zeros((samppix,samppix))
    vtheta[mask] = vr0*pow(r[mask]/(r0*dist), qr)

    # calculate projected velocity along the line of sight 
    vproj = np.zeros((samppix,samppix))
    vproj[mask] = vtheta[mask]*costheta[mask]*np.sin(incl)
    vproj = vproj + vlsr

    # Compute the synthetic datacube 
    mask = r != 0
    sigma = linewidth/(2*np.sqrt(2*np.log(2)))
    intensity = np.zeros((nchan, samppix, samppix)) 
    intensity[:,ydre:yure,xlre:xrre] = peakint[np.newaxis,ydre:yure,xlre:xrre] \
                *np.exp(-(vel[:,np.newaxis,np.newaxis] - vproj[ydre:yure,xlre:xrre])**2
                  /(2*sigma**2))

    # inner and outer radius
    for i in range(np.size(intensity,0)):
        mask = r < ri*dist
        intensity[i][mask] = 0
        mask = r > ro*dist
        intensity[i][mask] = 0

    # convolution 
    x_sigma = bmaj/sampsize / (2*np.sqrt(2*np.log(2)))
    y_sigma = bmin/sampsize / (2*np.sqrt(2*np.log(2)))
    bpa_r = bpa*math.pi/180.
    beam = Gaussian2DKernel(x_stddev=x_sigma, y_stddev=y_sigma, theta=bpa_r)
    
    for i in range(np.size(intensity,0)):
        if npix % 2 == 0:
            intensity[i,ydre:yure,xlre:xrre] = scipy_convolve(intensity[i,ydre:yure,xlre:xrre],beam,mode='same')
        else:
            #deal with singularity
            intensity[i,ydre:yure,xlre:xrre] = convolve(intensity[i,ydre:yure,xlre:xrre],beam)
                
    # resample pixel
    resinten = np.zeros((nchan, npix, npix))
    for k in range(ydown,yup):
        for n in range(xleft,xright):
            new = 0
            for i in range(samfrac):
                for j in range(samfrac):
                    new += intensity[:,k*samfrac+i,n*samfrac+j]
            resinten[:,k,n] = new    

    # scaling intensity 
    alpha = np.nansum(resinten*data)/np.nansum(pow(resinten,2))
    resinten = alpha*resinten
    
    return resinten

