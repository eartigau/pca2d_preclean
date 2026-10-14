#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 18:53:18 2026

@author: carmonan
"""
import astropy
from astropy.io import fits
from astropy.io import ascii
from astropy.stats import sigma_clip
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import copy
from scipy.optimize import curve_fit


def read_fits(file):
    hdu    = fits.open(file)
    data = hdu[1].data
    return data

def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

tcorr_file = sys.argv[1]
wave_file = sys.argv[2]
sky_mask_file = sys.argv[3]
file_out = sys.argv[4]

#tcorr_file="TEST_nanMasking/2539802o_pp_e2dsff_tcorr_AB.fits"
#wave_file="TEST_nanMasking/6420F52BD3a_pp_e2dsff_AB_wave_night_AB.fits"
#level=0.5
#file_out='e2dsff_masked.fits'
#sky_mask_file="3334061o_pp_e2dsff_AB_sky_mask.fits"



QC_PLOTS="NO"

print(f'tcorr file: {tcorr_file}')
print(f'wave file: {wave_file}')
print(f'file_out: {file_out}')

#READ THE FILE
spec2D=read_fits(tcorr_file)
mask_sky2D=read_fits(sky_mask_file)

spec2Dmasked=copy.deepcopy(spec2D)
mask2D_sky=copy.deepcopy(spec2D)*0.+1.



for order in range(0,48):
    
    spec=spec2D[order,:]
    spec_masked=spec2Dmasked[order,:]
    
    #mask sky
    mask_sky_number=mask_sky2D[order,:]
    mask_sky= mask_sky_number >= 1
    
    spec_masked[mask_sky]=float('NaN')
    
    if QC_PLOTS=="YES":                   
        fig = plt.figure()
        fig.set_size_inches(10, 8)
        plt.rcParams["font.size"] = (10)
        plt.subplots_adjust(left=0.1,
                            bottom=0.1,
                            right=0.9,
                            top=0.9,
                            wspace=0.1,
                            hspace=0.3)
        
        plt.plot(spec,label='spec',color='red')
        plt.plot(spec_masked,label='spec masked sky',color='blue')
        plt.title(f'order {order:.0f}')
    
    spec2Dmasked[order,:]=spec_masked
    
    
    if QC_PLOTS=="YES": 
        plt.plot(spec_masked,label='spec masked sky ',color='chartreuse')
   
    
    if QC_PLOTS=="YES":
        plt.legend(loc="upper right")
        plt.show()
    
   
    
    #continue_str=input("Continue? ")
    
#save
hdul         = fits.open(tcorr_file)
hdul[1].data = spec2Dmasked
hdul.writeto(file_out,overwrite=True)
print(file_out,"CREATED")
