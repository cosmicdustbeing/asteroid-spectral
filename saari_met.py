#!/usr/bin/env python3.9
# coding: utf-8

import numpy as np
import re
#import numba
#from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from loess.loess_1d import loess_1d
from astropy.io import ascii

#import statsmodels.nonparametric.smoothers_lowess

def band_params(wavx,refly,xnew):
    
    xout, ysmooth, weights = loess_1d(wavx, refly, xnew, frac=0.1, degree=2)

    plt.clf()
    plt.scatter(wavx, refly, c='limegreen', label='Data',s=5)
    plt.plot(xout, ysmooth, 'red', label='LOESS')
    plt.legend()
    plt.pause(1)
    plt.show(block=True)
    #plt.clf()
    #plt.plot(xout, BII_slopes)
    #plt.pause(1)

    # BAND I

    # where to look for the BI blue edge
    BI_bedge = np.where([(xnew > 0.63) & (xnew < 0.9)], ysmooth,0)

    # band I blue edge (first guess)
    BIb = [np.argmax(BI_bedge),np.max(BI_bedge)]

    # calculate slopes from blue edge
    BI_slopes = np.divide(BIb[1]-ysmooth,xnew[BIb[0]]-xnew, # avoid divide by zero
                          out=np.zeros_like(BIb[1]-ysmooth),where=xnew[BIb[0]]-xnew!=0)

    # where to look for the BI red edge
    BI_redge = np.where([(xnew > 1.25) & (xnew < 1.75)], BI_slopes,0)

    # band I red edge
    BIr = [np.argmax(BI_redge),ysmooth[np.argmax(BI_redge)]]

    # re-calculate slopes
    BI_slopes = np.divide(ysmooth-BIr[1],xnew-xnew[BIr[0]],
                          out=np.zeros_like(ysmooth-BIr[1]),where=xnew-xnew[BIr[0]]!=0)

    # where to look for the BI red edge 
    BI_bedge = np.where([(xnew > 0.63) & (xnew < 0.9)], BI_slopes,99.)

    # new band I blue edge
    BIb = [np.argmin(BI_bedge),ysmooth[np.argmin(BI_bedge)]]

    # calculate slope for output
    BIs = (BIr[1]-BIb[1])/(xnew[BIr[0]]-xnew[BIb[0]])/BIb[1]

    # wavelength, continuum and division of reflectance
    BI_wav = xnew[BIb[0]:BIr[0]]
    BI_cont = np.linspace(BIb[1], BIr[1], BIr[0]-BIb[0])
    BI_div = ysmooth[BIb[0]:BIr[0]]/BI_cont

    # band I area
    BIa = np.trapz(BI_div,BI_wav)

    # band I center
    BIc = BI_wav[np.argmin(BI_div)]
    BId = 1.-np.min(BI_div)

    # band I midpoint
    BIm = (np.max(np.where([BI_div < 1.-0.5*BId], BI_wav,0))+
           np.min(np.where([BI_div < 1.-0.5*BId], BI_wav,99.)))/2.

    # band I full width at half minimum
    BIfwhm = (np.max(np.where([BI_div < 1.-0.5*BId], BI_wav,0)) - BIm)*2.

    ##### BAND II

    # band I red edge
    BIIr = [np.argmax(xnew),ysmooth[np.argmax(xnew)]]

    # calculate slopes
    BII_slopes = np.divide(ysmooth-BIIr[1],xnew-xnew[BIIr[0]],
                           out=np.zeros_like(ysmooth-BIIr[1]),where=xnew-xnew[BIIr[0]]!=0)

    # where to look for the BII red edge
    BII_bedge = np.where([(xnew > 1.3) & (xnew < 1.9)], BII_slopes,99.)

    # new band I blue edge
    BIIb = [np.argmin(BII_bedge),ysmooth[np.argmin(BII_bedge)]]

    # calculate slope for output
    BIIs = (BIIr[1]-BIIb[1])/(xnew[BIIr[0]]-xnew[BIIb[0]])/BIIb[1]

    # wavelength, continuum and division of reflectance
    BII_wav = xnew[BIIb[0]:BIIr[0]]
    BII_cont = np.linspace(BIIb[1], BIIr[1], BIIr[0]-BIIb[0]) 
    BII_div = ysmooth[BIIb[0]:BIIr[0]]/BII_cont

    # band II area
    BIIa = np.trapz(BII_div,BII_wav)

    # band II center
    BIIc = BII_wav[np.argmin(BII_div)]
    BIId = 1.-np.min(BII_div)

    # band II midpoint
    BIIm = (np.max(np.where([BII_div < 1.-0.5*BIId], BII_wav,0))+
            np.min(np.where([BII_div < 1.-0.5*BIId], BII_wav,99.)))/2.

    # band II full width at half minimum
    BIIfwhm = (np.max(np.where([BII_div < 1.-0.5*BIId], BII_wav,0)) - BIIm)*2.

    output=[BIc,BIm,BIs,BIa,BIIc,BIIm,BIIs,BIIa]

    return output

# main 
if __name__ == '__main__':

# read in spectrum
    info = ascii.read('00080.EMMspex.txt')

    wavx  = np.array(info['col1'].data)
    refly = np.array(info['col2'].data)
#    err  = np.array(info['col3'].data)

    xnew = np.linspace(np.min(wavx), 2.45, 5000)

    # initialize arrays
#    N = 10000
#    mcerr = np.random.normal(0.,err,(N,len(wavx)))
#    mcBIc = np.zeros(N)
#    mcBIm = np.zeros(N)
#    mcBIs = np.zeros(N)
#    mcBIa = np.zeros(N)
#    mcBIIc = np.zeros(N)
#    mcBIIm = np.zeros(N)
#    mcBIIs = np.zeros(N)
#    mcBAR = np.zeros(N)

    # Monte Carlo loop
#    for i in range(N):

BIc,BIm,BIs,BIa,BIIc,BIIm,BIIs,BIIa = band_params(wavx,refly,xnew)

#results = Parallel(n_jobs=3)(delayed(band_params)(wavx,refly+mcerr[i,:],xnew,1./err^2.) for i in range(N))

#for j in range(N):
#    mcBIc[j] = results[j][0]
#    mcBIm[j] = results[j][1]
#    mcBIs[j] = results[j][2]
#    mcBIIc[j] = results[j][4]
#    mcBIIm[j] = results[j][5]
#    mcBIIs[j] = results[j][6]
#    mcBAR[j] = results[j][3]/results[j][7]

#        print(BIc,BIm,BIs,BIIc,BIIm,BIIs,BIa/BIIa)

print(BIc)
print(BIm)
print(BIs)
print(BIIc)
print(BIIm)
print(BIIs)
print(BIIa/BIa)

