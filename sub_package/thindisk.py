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

####################
## thindisk model ##
####################
def alg(x):
    return x/np.sqrt(1+x**2)

def thindisk(vlsr,incl,pa,vr0sin,linewidth,ri,ro,imgpar):
    # input image parameters
    data,npix,nchan,velwidth,velstart = imgpar
    
    # ra,dec,vel grids on the plane of the sky 
    ra_offset = (np.arange(npix)-(npix-1)/2 - x0)*pixsize
    dec_offset = (np.arange(npix)-(npix-1)/2 - y0)*pixsize
    ra_grid, dec_grid = np.meshgrid(ra_offset,dec_offset)
    vel = (np.arange(nchan)*velwidth + velstart)/1000 # km/s

    # x,y grids on the plane of the disk 
    pa = pa*math.pi/180.
    x_grid = +ra_grid*np.cos(pa) + dec_grid*np.sin(pa) 
    y_grid = -ra_grid*np.sin(pa) + dec_grid*np.cos(pa)
    theta = np.zeros((npix, npix))
    r = np.zeros((npix, npix))
    mask = y_grid != 0
    incl = incl*math.pi/180.
    theta[mask] = 2*np.arctan((y_grid[mask]/np.cos(incl)) \
                  /(x_grid[mask] \
                  + np.sqrt(x_grid[mask]**2 + (y_grid[mask]/np.cos(incl))**2)))
    r = np.sqrt(x_grid**2 + (y_grid/np.cos(incl))**2)*dist # AU

    # calculate line peak intensity 
    peakint = np.zeros((npix, npix))
    mask = r != 0
    peakint[mask] = 10*pow(r[mask]/(r0*dist), qd)

    # calculate disk velocity 
    vtheta = np.zeros((npix,npix))
    vtheta[mask] = (vr0sin/np.sin(incl))*pow(r[mask]/(r0*dist), qr)

    # calculate projected velocity along the line of sight 
    vproj = np.zeros((npix,npix))
    vproj[mask] = vtheta[mask]*np.cos(theta[mask])*np.sin(incl)
    vproj = vproj + vlsr

    # Compute the synthetic datacube 
    mask = r != 0
    sigma = np.zeros((npix, npix))
    sigma = linewidth/(2*np.sqrt(2*np.log(2)))
    intensity = np.zeros((nchan, npix, npix)) 
    intensity = peakint[np.newaxis,:,:] \
                *np.exp(-(vel[:,np.newaxis,np.newaxis] - vproj)**2/(2*sigma**2))

    # inner and outer radius
    for i in range(np.size(intensity,0)):
        zeroint = np.zeros((npix,npix))
        intensity[i,:,:] = zeroint*(alg(r-ro*dist)+1)/2 \
                           + intensity[i,:,:]*(1-alg(r-ro*dist))/2
        intensity[i,:,:] = intensity[i,:,:]*(alg(r-ri*dist)+1)/2 \
                           + zeroint*(1-alg(r-ri*dist))/2

    # convolution 
    x_sigma = bmaj/pixsize / (2*np.sqrt(2*np.log(2)))
    y_sigma = bmin/pixsize / (2*np.sqrt(2*np.log(2)))
    bpa_r = bpa*math.pi/180.
    beam = Gaussian2DKernel(x_stddev=x_sigma, y_stddev=y_sigma, theta=bpa_r)
    
    for i in range(np.size(intensity,0)):
        if npix % 2 == 0:
            intensity[i,:,:] = scipy_convolve(intensity[i,:,:],beam,mode='same')
        else:
            intensity[i,:,:] = convolve(intensity[i,:,:],beam) #deal with singularity
    
    # scaling intensity 
    alpha = np.nansum(intensity*data)/np.nansum(pow(intensity,2))
    intensity = alpha*intensity
    
    return intensity

