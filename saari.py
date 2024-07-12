# coding: utf-8

import numpy as np
import argparse
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from loess.loess_1d import loess_1d
from astropy.io import ascii

def band_params(wavx, refly, xnew, wt):
    xout, ysmooth, weights = loess_1d(wavx, refly, xnew, frac=0.1, degree=2, sigy=wt)

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

    # band area ratio
    BAR = BIIa/BIa
    
    output = [xnew[BIb[0]], xnew[BIr[0]], BIc, BIm, BIs, BId, BIa, BIfwhm,
              xnew[BIIb[0]], xnew[BIIr[0]], BIIc, BIIm, BIIs, BIId, BIIa, BIIfwhm, BAR]

    return output

# main 
if __name__ == '__main__':
    # read in spectrum
    parser = argparse.ArgumentParser(description="Process asteroid spectrum data.")
    parser.add_argument('input_file', type=str, help="Path to the input text file containing spectrum data.")
    parser.add_argument('output_file', type=str, help="Path to the output text file to save results.")
    args = parser.parse_args()

    # read in spectrum
    info = ascii.read(args.input_file)
    
    wavx = np.array(info['col1'].data)
    refly = np.array(info['col2'].data)
    err = np.array(info['col3'].data)
    xnew = np.linspace(np.min(wavx), 2.45, 5000)

    # initialize arrays
    N = 1000
    mcerr = np.random.normal(0., err, (N, len(wavx)))
    
    # Parallel execution
    results = Parallel(n_jobs=-2)(delayed(band_params)(wavx, refly + mcerr[i, :], xnew, 1. / err) for i in range(N))

    mcBIb = np.zeros(N)
    mcBIr = np.zeros(N)
    mcBIc = np.zeros(N)
    mcBIm = np.zeros(N)
    mcBIs = np.zeros(N)
    mcBId = np.zeros(N)
    mcBIa = np.zeros(N)
    mcBIw = np.zeros(N)
    mcBIIb = np.zeros(N)
    mcBIIr = np.zeros(N)
    mcBIIc = np.zeros(N)
    mcBIIm = np.zeros(N)
    mcBIIs = np.zeros(N)
    mcBIId = np.zeros(N)
    mcBIIa = np.zeros(N)
    mcBIIw = np.zeros(N)
    mcBAR = np.zeros(N)

    for j in range(N):
        mcBIb[j] = results[j][0]
        mcBIr[j] = results[j][1]
        mcBIc[j] = results[j][2]
        mcBIm[j] = results[j][3]
        mcBIs[j] = results[j][4]
        mcBId[j] = results[j][5]
        mcBIa[j] = results[j][6]
        mcBIw[j] = results[j][7]
        mcBIIb[j] = results[j][8]
        mcBIIr[j] = results[j][9]
        mcBIIc[j] = results[j][10]
        mcBIIm[j] = results[j][11]
        mcBIIs[j] = results[j][12]
        mcBIId[j] = results[j][13]
        mcBIIa[j] = results[j][14]
        mcBIIw[j] = results[j][15]
        mcBAR[j] = results[j][16]

    # Save results to a text file with 3 significant figures
    with open(args.output_file, 'w') as f:
        f.write("BIb\tBIr\tBIc\tBIm\tBIs\tBId\tBIw\tBIIb\tBIIc\tBIIm\tBIIs\tBIId\tBIIw\tBAR\n")
        for i in range(N):
            f.write(f"{mcBIb[i]:.3g}\t{mcBIr[i]:.3g}\t{mcBIc[i]:.3g}\t{mcBIm[i]:.3g}\t{mcBIIs[i]:.3g}\t{mcBId[i]:.3g}\t{mcBIw[i]:.3g}\t{mcBIIb[i]:.3g}\t{mcBIIc[i]:.3g}\t{mcBIIm[i]:.3g}\t{mcBIIs[i]:.3g}\t{mcBIId[i]:.3g}\t{mcBIIw[i]:.3g}\t{mcBAR[i]:.3g}\n")

        
    print(f"BI center: {np.mean(mcBIc):.3g} ± {np.std(mcBIc):.3g}")
    print(f"BI slope: {np.mean(mcBIs):.3g} ± {np.std(mcBIs):.3g}")
    print(f"BI depth: {np.mean(mcBId):.3g} ± {np.std(mcBId):.3g}")
    print(f"BII midpoint: {np.mean(mcBIIm):.3g} ± {np.std(mcBIIm):.3g}")
    print(f"BII slope: {np.mean(mcBIIs):.3g} ± {np.std(mcBIIs):.3g}")
    print(f"BII depth: {np.mean(mcBIId):.3g} ± {np.std(mcBIId):.3g}")
    print(f"B.A.R.: {np.mean(mcBAR):.3g} ± {np.std(mcBAR):.3g}")

    # Plot histograms of all parameters
    fig, axs = plt.subplots(3, 3, figsize=(8, 8))
    fig.suptitle('Histograms of Band Parameters')

    axs[0, 0].hist(mcBIc, bins=50, color='blue', edgecolor='black')
    axs[0, 0].set_title('BI Center (mcBIc)')

    axs[0, 1].hist(mcBIs, bins=50, color='blue', edgecolor='black')
    axs[0, 1].set_title('BI Slope (mcBIm)')

    axs[0, 2].hist(mcBId, bins=50, color='blue', edgecolor='black')
    axs[0, 2].set_title('BI Depth')

    axs[1, 0].hist(mcBIIm, bins=50, color='blue', edgecolor='black')
    axs[1, 0].set_title('BII Midpoint (mcBIIc)')

    axs[1, 1].hist(mcBIIs, bins=50, color='blue', edgecolor='black')
    axs[1, 1].set_title('BII Slope (mcBIIm)')

    axs[1, 2].hist(mcBIId, bins=50, color='blue', edgecolor='black')
    axs[1, 2].set_title('BII Depth (mcBIIs)')

    axs[2, 0].hist(mcBAR, bins=50, color='blue', edgecolor='black')
    axs[2, 0].set_title('Band Area Ratio (mcBAR)')

    axs[2, 1].hist(mcBIw, bins=50, color='blue', edgecolor='black')
    axs[2, 1].set_title('BI Width (mcBAR)')

    axs[2, 2].hist(mcBIIw, bins=50, color='blue', edgecolor='black')
    axs[2, 2].set_title('BII Width (mcBAR)')
    
    # Remove empty subplots
    #fig.delaxes(axs[2, 1])
    #fig.delaxes(axs[2, 2])

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

