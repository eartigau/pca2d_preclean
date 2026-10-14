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
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


#n_sigma_vec=[1,2,3,4,5,6,7,8,9,10,15,20]
n_sigma_vec=[3]
QC_PLOTS="YES"
sky_file="sky_files/3334061o_pp_e2dsff_AB.fits"

for n_sigma in n_sigma_vec:
    
    file_out="masks/3334061o_pp_e2dsff_AB_sky_mask_"+f'{n_sigma:0.0f}'+"sigma.fits"
    
    
    print(f'sky file: {sky_file}')
    print(f'file_out: {file_out}')
    
    #READ THE FILE
    sky2D=read_fits(sky_file)
    
    mask2D_sky=copy.deepcopy(sky2D)*0.+1.
    
    
    for order in range(0,48):
    
    
        if QC_PLOTS=="YES":
            fig = plt.figure()
            fig.set_size_inches(10, 8)
            plt.rcParams["font.size"] = (10)
            plt.subplots_adjust(left=0.1,
                                bottom=0.1,
                                right=0.9,
                                top=0.9,
                                wspace=0.1,
                                hspace=0.5)
            
            
         
        #read the spectrum
        sky_spec=copy.deepcopy(sky2D[order,:])
        
        
        #continuum normalization
        XX = np.arange(0,len(sky_spec),1.)
        flux_median=np.nanmedian(sky_spec)
        flux_std=np.nanstd(sky_spec)
        
        mask=sky_spec < 2*flux_median
        coefficients = np.polyfit(XX[mask], sky_spec[mask], 5.)
        # Create polynomial function
        p = np.poly1d(coefficients)     
        #print("Linear Fit Coefficients:", coefficients)
        continuum=p(XX)
        
        if QC_PLOTS=="YES":
            plt.subplot(3,1,1)
            plt.plot(XX,sky_spec,color='red',label='sky')
            #plt.ylim(0,2000)    
            plt.plot(XX, continuum,color='chartreuse',label='sky')
            plt.title(f'{sky_file} order {order:.0f}')
        
        sky_spec=sky_spec-continuum
                
                
        # find median and std of the sky
        flux_median=np.nanmedian(sky_spec)
        flux_std=np.nanstd(sky_spec)
        
        #define 
        mask=(sky_spec < flux_median + 5*flux_std)
        
        
        #make histogram using the flux masked
        #counts, bin_edges, patches = plt.hist(sky_spec[mask], density=True, bins=100)  # density=False would make counts
        counts, bin_edges = np.histogram(sky_spec[mask], density=True,bins=100)
        
        #construc vectors to fit a gaussian to the continuum flux
        XX=[]
        YY=counts
        
        #construct XX vector from histrogram
        for pixel in range(len(bin_edges)-1):
            x_local=0.5*(bin_edges[pixel]+bin_edges[pixel+1])
            XX.append(x_local)
        
        #all works using numpy arrays
        XX=np.array(XX)
        YY=np.array(YY)
        
        
        
        
        #fit the gaussian
        n = len(XX)
        mean = sum(XX*YY)/n                   #note this correction
        sigma = 50 #sum(YY*(XX-mean)**2)/n
        popt,pcov = curve_fit(gaus,XX,YY,p0=[1,mean,sigma])
        YFIT=gaus(XX,*popt)
        
        
        
        a=popt[0]
        x0=popt[1]
        sigma=np.abs(popt[2])
        
        # Define the upper value of masking
        # print("a,x0,sigma",popt)
        upper_value=x0+n_sigma*sigma
        
        
        #check the fit
        if QC_PLOTS=="YES":
            # # Histogram of the 3 sigma kippling sky
            plt.subplot(3,1,2)
            plt.stairs(counts, bin_edges,color='blue',fill=True)
            #plt.axvline(flux_median , linestyle='--', color='black')
            #plt.axvline(flux_median + 3*flux_std , linestyle='--', color='black')
            plt.plot(XX,YY,color='red')
            plt.ylabel('Probability')
            plt.xlabel('flux [ADU/s]')
            plt.title(f'order {order:.0f} histogram sky flux $\sigma$={n_sigma:.0f}')
            plt.plot(XX,YFIT,color='chartreuse')
            plt.axvline(x0, linestyle='--', color='yellow')
            plt.axvline(x0+sigma, linestyle='--', color='yellow')
            plt.axvline(upper_value, linestyle='--', color='chartreuse')
            #xmax=max([1000,x0+(n_sigma+1)*sigma])
            plt.xlim(x0-4*sigma,x0+(n_sigma+2)*sigma)
            #plt.xlim(flux_median-n_sigma*flux_std,flux_median+n_sigma*flux_std)
        
        
        #*****
        # perform the masking
        #****
        
        mask_sky=(sky_spec > upper_value)
        
        #we load it again because we put out the continuum
        sky_spec=copy.deepcopy(sky2D[order,:])
        
        sky_masked=copy.deepcopy(sky2D[order,:])
        sky_masked[mask_sky]=float('NaN')
        
        #store the mask in a 2D array
        mask2D_sky[order,:]=mask_sky
        
        if QC_PLOTS=="YES":
            plt.subplot(3,1,3)
            plt.plot(sky_spec,color='red',label='sky')
            plt.plot(sky_masked,color='blue',label='sky masked')
            plt.title(f'order {order:.0f} sky spectrum filtered $\sigma$={n_sigma:.0f}')
            #plt.axhline(upper_value, linestyle='--', color='chartreuse')
            #plt.axhline(x0, linestyle='--', color='yellow')
            plt.plot(continuum+upper_value, linestyle='--', color='chartreuse')
            plt.plot(continuum+x0, linestyle='--', color='yellow')
            plt.legend(loc="upper right")
            #plt.ylim(x0-3*sigma,x0+10*sigma)
            flux_median=np.nanmedian(sky_spec)
            plt.ylim(0.5*flux_median,np.max([500,2*flux_median]))
            #plt.ylim(0.5*flux_median,1.5*flux_median)
        
        plt.show()
        
        #continue_str=input("Continue? ")
        
        
    #save the mask
    hdul         = fits.open(sky_file)
    hdul[1].data = mask2D_sky
    hdul.writeto(file_out,overwrite=True)
    print(file_out,"CREATED")    
