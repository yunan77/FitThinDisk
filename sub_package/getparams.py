#!/usr/local/bin/python3
import sys
import numpy as np
from astropy.io import fits
import configparser
import json

class ReadCube:
    def __init__(self,imgname):
        image = fits.open(imgname)
        self.data = image[0].data
        self.hd = image[0].header
        self.npix = self.hd['NAXIS1']
        self.nchan = self.hd['NAXIS3']
        self.restfreq = self.hd['RESTFRQ']
        if self.hd['CTYPE3'] == 'VRAD':
            self.velwidth = self.hd['CDELT3']
            self.velstart = self.hd['CRVAL3']
            self.velcrpix = self.hd['CRPIX3']
        elif self.hd['CTYPE3'] == 'FREQ':
            self.freqwidth = self.hd['CDELT3']
            self.freqstart = self.hd['CRVAL3']
            self.freqcrpix = self.hd['CRPIX3']

    def calrms(self):
        rms = np.zeros((self.nchan,self.npix,self.npix))
        for i in range(self.nchan):
            rms[i] = np.sqrt(np.nanmean(np.power(self.data[i],2)))
        self.rms = rms

    def select_chans(self,first=0,last=99):
        if len(self.data.shape) == 4:
            if np.size(self.data,0) != 1:
                self.data = self.data[first:last+1,0,:,:]
                self.rms = self.rms[first:last+1,0,:,:]
            else:
                self.data = self.data[0,first:last+1,:,:]
                self.rms = self.rms[0,first:last+1,:,:]
        else:
            self.data = self.data[first:last+1,:,:]
            self.rms = self.rms[first:last+1,:,:]
        self.nchan = last-first+1
        if self.hd['CTYPE3'] == 'VRAD':
            self.velstart = self.velstart + self.velwidth*first
            self.hd['CRVAL3'] = self.velstart
        elif self.hd['CTYPE3'] == 'FREQ':
            self.freqstart = self.freqstart + self.freqwidth*first
            self.hd['CRVAL3'] = self.freqstart

class Parameters:
    def __init__(self,inifile):
        para = configparser.ConfigParser()
        para.read(inifile)
        self.imgpath = para.get('file','imgpath')
        try:
            self.imglist = json.loads(para.get('file','imglist'))
        except:
            imglist = para.get('file','imglist').replace("'",'"')
            self.imglist = json.loads(imglist)
        try:
            self.chanrange = json.loads(para.get('file','chanrange'))
        except:
            chanrange = para.get('file','chanrange').replace("'",'"')
            self.chanrange = json.loads(chanrange)
        pixsize = para.getfloat('fixed', 'pixsize')
        samplesize = para.getfloat('fixed', 'samplesize')
        imgrange = para.get('fixed', 'imgrange').split(',')
        dist = para.getfloat('fixed', 'dist')
        qr = para.getfloat('fixed', 'qr')
        qd = para.getfloat('fixed', 'qd')
        r0 = para.getfloat('fixed', 'r0')
        bmaj = para.getfloat('fixed', 'bmaj')
        bmin = para.getfloat('fixed', 'bmin')
        bpa = para.getfloat('fixed', 'bpa')
        x0 = para.getfloat('fixed', 'x0')
        y0 = para.getfloat('fixed', 'y0')
        vlsr = para.get('free', 'vlsr').split(',')
        incl = para.get('free', 'incl').split(',')
        pa = para.get('free', 'pa').split(',')
        vr0 = para.get('free', 'vr0').split(',')
        linewidth = para.get('free', 'linewidth').split(',')
        ri = para.get('free', 'ri').split(',')
        ro = para.get('free', 'ro').split(',')
        self.fixed = {'pixsize':pixsize,'samplesize':samplesize,'imgrange':imgrange,
                      'dist':dist,'qr':qr,'qd':qd,'r0':r0,'bmaj':bmaj,'bmin':bmin,
                      'bpa':bpa,'x0':x0,'y0':y0
                     }
        try:
            self.free = {
            'vlsr':{'value':float(vlsr[0]),'vary':vlsr[1],'bound':[float(vlsr[2]),float(vlsr[3])]},
            'incl':{'value':float(incl[0]),'vary':incl[1],'bound':[float(incl[2]),float(incl[3])]},
            'pa':{'value':float(pa[0]),'vary':pa[1],'bound':[float(pa[2]),float(pa[3])]},
            'vr0':{'value':float(vr0[0]),'vary':vr0[1],'bound':[float(vr0[2]),float(vr0[3])]},
            'linewidth':{'value':float(linewidth[0]),'vary':linewidth[1],'bound':[float(linewidth[2]),float(linewidth[3])]},
            'ri':{'value':float(ri[0]),'vary':ri[1],'bound':[float(ri[2]),float(ri[3])]},
            'ro':{'value':float(ro[0]),'vary':ro[1],'bound':[float(ro[2]),float(ro[3])]},
                         }
        except IndexError:
            print('Please specify initial value, variation, minimum value and maximum value')

    def getinitial(self):
        self.para0 = []
        self.lowbound = []
        self.upbound = []
        self.fixval = []
        self.fitname = []
        self.nofit = {}
        for key, value in self.free.items():
            if value['vary'] == 'True':
                self.para0.append(value['value'])
                self.lowbound.append(value['bound'][0])
                self.upbound.append(value['bound'][1])
                self.fitname.append(key)
            elif value['vary'] == 'False':
                self.nofit['{}'.format(key)] = value['value']
            else:
                sys.exit("InputError: 'variation' should be 'True' or 'False'")
        self.df = len(self.para0)
