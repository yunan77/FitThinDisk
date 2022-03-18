#!/usr/local/bin/python3
import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares
from scipy import linalg, stats
from sub_package.thindisk import thindisk
from sub_package.getparams import ReadCube,Parameters

#################################################
## please setup input parameters in params.ini ##
#################################################

## read parameters in params.ini ##
params = Parameters('params.ini')
params.getinitial()
xleft = int(params.fixed['imgrange'][0])
xright = int(params.fixed['imgrange'][1])+1
ydown = int(params.fixed['imgrange'][2])
yup = int(params.fixed['imgrange'][3])+1

## define residual between data and model ##
def residual(para,data,npix,nchan,velwidth,velstart,rms):
    i=0
    for name in params.fitname:
        exec(name+'='+str(para[i]),globals())
        i += 1
    for key,value in params.nofit.items():
        if key == 'ri' or key == 'ro':
            exec(key+'='+'False',globals())
        exec(key+'='+str(value),globals())

    resid = np.array([])
    for i in range(len(nchan)):
        imgpar = data[i],npix[i],nchan[i],velwidth[i],velstart[i]
        intensity = thindisk(vlsr,incl,pa,vr0,linewidth,ri,ro,imgpar)
        resid = np.append(resid,((data[i][:,ydown:yup,xleft:xright]-intensity[:,ydown:yup,xleft:xright])/rms[i][:,ydown:yup,xleft:xright]).reshape(-1)) 
   
    return resid
        
## read image parameters ##
chi = []
chf = []
data = []
npix = []
nchan = []
velwidth = []
velstart = []
rms = []
hd = []

for i in range(len(params.chanrange)):
    for j in range(len(params.chanrange[i])):
        if params.chanrange[i][j] == '~':
            chi.append(int(params.chanrange[i][0:j]))
            chf.append(int(params.chanrange[i][j+1:len(params.chanrange[i])+1]))

totalchan = 0
for i in range(len(params.imglist)):
    cube =  ReadCube(params.imgpath + params.imglist[i])
    cube.calrms()
    cube.select_chans(first=chi[i],last=chf[i])
    data.append(cube.data)
    npix.append(cube.npix)
    nchan.append(cube.nchan)
    velwidth.append(cube.velwidth)
    velstart.append(cube.velstart)
    rms.append(cube.rms)
    hd.append(cube.hd)
    totalchan += cube.nchan

## find minimun of residual ##
res = least_squares(residual, params.para0, method='trf',diff_step=1e-5,args=(data,npix,nchan,velwidth,velstart,rms))

## compute uncertainty for fit parameters ##
U, s, Vh = linalg.svd(res.jac, full_matrices=False)
tol = np.finfo(float).eps*s[0]*max(res.jac.shape)
w = s > tol
cov = (Vh[w].T/s[w]**2) @ Vh[w]  # robust covariance matrix
s_sq = np.sum(res.fun**2)/(res.fun.size - res.x.size)
cov *= s_sq
nbeam = params.fixed['bmaj']*params.fixed['bmin']/params.fixed['pixsize']**2
effpix = (xright-xleft)*(yup-ydown)*totalchan/nbeam
chi2 = stats.chi2.ppf(0.683,params.df+1)
perr = np.sqrt(np.diag(cov))*np.sqrt(nbeam)*np.sqrt(chi2)

## print fitting result ##
out = []
add = out.append
if res.success is True:
    add('Fitting success!')
    add(f'reduced chi square = {s_sq}')
    add(f'optimized with {effpix} effective pixels')
    add('[Variables]')
    for i in range(params.df):
        init = '(init = {})'.format(params.para0[i])
        add(f'   {params.fitname[i]:12}={res.x[i]:10.3f} +/- {perr[i]:<6.3f} {init:15} (bounds = [{params.lowbound[i]},{params.upbound[i]}])')
    for key,value in params.nofit.items():
        add(f"   {key:12}={value:10.3f}            (fixed)")
else:
    print('Fitting does not converge. Please try different initial values or bounds.')

print('\n'.join(out))

## output model fits ##
i=0
for name in params.fitname:
    exec(name + '=res.x[i]')
    i += 1

for i in range(len(params.imglist)):
    imgpar = data[i],npix[i],nchan[i],velwidth[i],velstart[i]
    inten = thindisk(vlsr,incl,pa,vr0,linewidth,ri,ro,imgpar)
    model = fits.PrimaryHDU(inten)
    hdmodel = model.header
    hdmodel["CTYPE1"] = hd[i]["CTYPE1"]
    hdmodel["CTYPE2"] = hd[i]["CTYPE2"]
    hdmodel["CTYPE3"] = hd[i]["CTYPE3"]
    hdmodel["CUNIT1"] = hd[i]["CUNIT1"]
    hdmodel["CUNIT2"] = hd[i]["CUNIT2"]
    hdmodel["CUNIT3"] = hd[i]["CUNIT3"]
    hdmodel["CDELT1"] = hd[i]["CDELT1"]
    hdmodel["CDELT2"] = hd[i]["CDELT2"]
    hdmodel["CDELT3"] = hd[i]["CDELT3"]
    hdmodel["CRPIX1"] = hd[i]["CRPIX1"]
    hdmodel["CRPIX2"] = hd[i]["CRPIX2"]
    hdmodel["CRPIX3"] = hd[i]["CRPIX3"]
    hdmodel["CRVAL1"] = hd[i]["CRVAL1"]
    hdmodel["CRVAL2"] = hd[i]["CRVAL2"]
    hdmodel["CRVAL3"] = hd[i]["CRVAL3"]
    hdmodel["BUNIT"] = hd[i]["BUNIT"]
    hdmodel["RESTFRQ"] = hd[i]["RESTFRQ"]
    hdmodel["VELO-LSR"] = vlsr
    model.writeto('modelfits/'+params.imglist[i]+'.model', overwrite = True)
