# coding: utf-8

import numpy as np
import argparse
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from loess.loess_1d import loess_1d
from astropy.io import ascii

def band_params(wavx, refly, xnew):
    xout, ysmooth, weights = loess_1d(wavx, refly, xnew, frac=0.1, degree=2)

    # BAND I

    # first guess at blue edge (maximum reflectance)
    BI_bedge = np.where((xnew > 0.63) & (xnew < 0.9), ysmooth, 0)

    # index and smoothed reflectance of BI blue edge
    BIb = [np.argmax(BI_bedge), np.max(BI_bedge)]

    # calculate slopes from blue edge
    BI_slopes = np.divide(BIb[1] - ysmooth, xnew[BIb[0]] - xnew,
                          out=np.zeros_like(BIb[1] - ysmooth), where=xnew[BIb[0]] - xnew != 0)

    # where to look for the BI red edge
    BI_redge = np.where((xnew > 1.) & (xnew < 1.75), BI_slopes, 0)

    # band I red edge
    BIr = [np.argmax(BI_redge), ysmooth[np.argmax(BI_redge)]]

    # re-calculate slopes
    BI_slopes = np.divide(ysmooth - BIr[1], xnew - xnew[BIr[0]],
                          out=np.zeros_like(ysmooth - BIr[1]), where=xnew - xnew[BIr[0]] != 0)

    # where to look for the BI red edge
    BI_bedge = np.where((xnew > 0.63) & (xnew < 0.9), BI_slopes, 99.)

    # new band I blue edge
    BIb = [np.argmin(BI_bedge), ysmooth[np.argmin(BI_bedge)]]

    # calculate slope for output
    BIs = (BIr[1] - BIb[1]) / (xnew[BIr[0]] - xnew[BIb[0]]) / BIb[1]

    # wavelength, continuum and division of reflectance
    BI_wav = xnew[BIb[0]:BIr[0]+1]
    BI_cont = np.linspace(BIb[1], BIr[1], BIr[0] - BIb[0] + 1)
    BI_div = ysmooth[BIb[0]:BIr[0]+1] / BI_cont

    # band I area
    BIa = np.trapz(BI_div, BI_wav)

    # band I center
    BIc = BI_wav[np.argmin(BI_div)]
    BId = 1. - np.min(BI_div) # band depth

    # band I midpoint
    BIm = (np.max(np.where(BI_div < 1. - 0.5 * BId, BI_wav, 0)) +
           np.min(np.where(BI_div < 1. - 0.5 * BId, BI_wav, 99.))) / 2.

    # band I full width at half minimum
    BIfwhm = (np.max(np.where(BI_div < 1. - 0.5 * BId, BI_wav, 0)) - BIm) * 2.

    # BAND II
    BII_redge = np.argmax(xnew)
    BIIr = [BII_redge, ysmooth[BII_redge]]
    BII_slopes = np.divide(ysmooth - BIIr[1], xnew - xnew[BIIr[0]],                              
                           out=np.zeros_like(ysmooth - BIIr[1]), where=xnew - xnew[BIIr[0]] != 0)
    BII_bedge = np.where((xnew > 1.3) & (xnew < 1.7), BII_slopes, 99.)
    BIIb = [np.argmin(BII_bedge), ysmooth[np.argmin(BII_bedge)]]
    BIIs = (BIIr[1] - BIIb[1]) / (xnew[BIIr[0]] - xnew[BIIb[0]]) / BIIb[1]

    BII_wav = xnew[BIIb[0]:BIIr[0]+1]
    BII_cont = np.linspace(BIIb[1], BIIr[1], BIIr[0] - BIIb[0] + 1)
    BII_div = ysmooth[BIIb[0]:BIIr[0]+1] / BII_cont
    BIIa = np.trapz(BII_div, BII_wav)
    BIIc = BII_wav[np.argmin(BII_div)]
    BIId = 1. - np.min(BII_div)
    BIIm = (np.max(np.where(BII_div < 1. - 0.5 * BIId, BII_wav, 0)) +
            np.min(np.where(BII_div < 1. - 0.5 * BIId, BII_wav, 99.))) / 2.
    BIIfwhm = (np.max(np.where(BII_div < 1. - 0.5 * BIId, BII_wav, 0)) - BIIm) * 2.

    # band are ratio
    BAR = BIIa/BIa
    
    output = [xnew[BIb[0]], xnew[BIr[0]], BIc, BIm, BIs, BId, BIa, BIfwhm,
              xnew[BIIb[0]], xnew[BIIr[0]], BIIc, BIIm, BIIs, BIId, BIIa, BIIfwhm, BAR]

    return output

# main 
if __name__ == '__main__':
    # read in spectrum
    parser = argparse.ArgumentParser(description="Process meteoroid spectrum data.")
    parser.add_argument('input_file', type=str, help="Path to the input text file containing spectrum data.")
    parser.add_argument('output_file', type=str, help="Path to the output text file to save results.")
    args = parser.parse_args()

    # read in spectrum
    info = ascii.read(args.input_file)
    
    wavx = np.array(info['col1'].data)
    refly = np.array(info['col2'].data)
    xnew = np.linspace(np.min(wavx), 2.45, 5000)

    # get the parameters
    results = band_params(wavx, refly, xnew)

    BIb = results[0]
    BIr = results[1]
    BIc = results[2]
    BIm = results[3]
    BIs = results[4]
    BId = results[5]
    BIa = results[6]
    BIw = results[7]
    BIIb = results[8]
    BIIr = results[9]
    BIIc = results[10]
    BIIm = results[11]
    BIIs = results[12]
    BIId = results[13]
    BIIa = results[14]
    BIIw = results[15]
    BAR = results[16]

    # Save results to a text file
    with open(args.output_file, 'w') as f:
        f.write("BIc\tBIm\tBIs\tBIa\tBIIc\tBIIm\tBIIs\tBAR\n")
        f.write(f"{BIc:.3g}\t{BIm:.3g}\t{BIs:.3g}\t{BAR:.3g}\t{BIIc:.3g}\t{BIIm:.3g}\t{BIIs:.3g}\t{BAR:.3g}\n")
    
