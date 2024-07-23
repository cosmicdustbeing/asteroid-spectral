import argparse
import pandas as pd
import numpy as np
from astropy.io import ascii
from scipy.integrate import simpson
from astropy.io import fits

#def read_spectrum(filename):
#    # Load the data
#    data = ascii.read(filename, comment='#',
#                      names=['wavelength','flux','err','flag'])
##    data = pd.read_csv(filename, delim_whitespace=True, comment='#', header=None)
##    wavelength = data[0].values  # in micrometers
##    flux = data[2].values  # flux values
##    wavelength = data['wavelength'].value  # in micrometers
##    flux = data['flux'].value  # flux values       
#
#    return data['wavelength'].value, data['flux'].value

def read_spectrum_from_fits(filename):
    with fits.open(filename) as hdul:
        data = hdul[0].data
        wavelength = data[0,]
        flux = data[1,]
        error = data[2,]
    return wavelength, flux, error

def read_header_parameters_from_fits(filename):
    parameters = {}
    with fits.open(filename) as hdul:
        header = hdul[0].header
        parameters['mjd'] = header.get('MJD_OBS', 'N/A')
        parameters['hourangle'] = header.get('TCS_HA', 'N/A')
        parameters['azimuth'] = header.get('TCS_AZ', 'N/A')
        parameters['elevation'] = header.get('TCS_EL', 'N/A')
        parameters['airmass'] = header.get('AIRMASS', 'N/A')
        parameters['humidity'] = header.get('HUMIDITY', 'N/A')
        parameters['position_angle'] = header.get('POSANGLE', 'N/A')
        parameters['airtemp'] = header.get('TEMPERAT', 'N/A')

    return parameters

def read_filter(filename):
    # Load the filter transmission data

    filter_data = ascii.read('j_filter.txt', names=['wavelength','transmission'])
#    filter_data = pd.read_csv(filename, delim_whitespace=True, comment='#', header=None, names=['wavelength', 'transmission'])
#    filter_data['wavelength'] /= 10000  # Convert wavelength to micrometers
    return filter_data

def calculate_magnitudes(wavelength, flux, filter_data, zp):
    transmission = np.interp(wavelength, filter_data['wavelength'], filter_data['transmission'])
    integrated_flux = simpson(flux * transmission, x=wavelength)
    magnitude = -2.5 * np.log10(integrated_flux / zp )
    return magnitude

def main():
    parser = argparse.ArgumentParser(description="Calculate 2MASS magnitudes and colors from a spectrum file")
    parser.add_argument('inputfile', type=str, help="Input spectrum file")
    parser.add_argument('outputfile', type=str, help="Output file for magnitudes and colors")
    args = parser.parse_args()

    # Read spectrum data
    wavelength, flux, error = read_spectrum_from_fits(args.inputfile)
    
    # Read airmass from the header
    parameters = read_header_parameters_from_fits(args.inputfile)
    
    # Load filter transmission data
    j_filter = read_filter('j_filter.txt')
    h_filter = read_filter('h_filter.txt')
    k_filter = read_filter('k_filter.txt')

    # 2MASS zero points
    zp_j = 1594
    zp_h = 1024
    zp_k = 666.7

    # Calculate magnitudes
    m_j = calculate_magnitudes(wavelength, flux, j_filter, zp_j)
    m_h = calculate_magnitudes(wavelength, flux, h_filter, zp_h)
    m_k = calculate_magnitudes(wavelength, flux, k_filter, zp_k)

    # Calculate colors
    color_jh = m_j - m_h
    color_hk = m_h - m_k
    color_jk = m_j - m_k
    
    output_data = {
        'MJD' : [parameters.get('mjd', 'N/A')],
        'HourAngle' : [parameters.get('hourangle', 'hr')],
        'Azimuth' : [parameters.get('azimuth', 'deg')],
        'Elevation' : [parameters.get('elevation', 'deg')],
        'Airmass': [parameters.get('airmass', 'N/A')],
        'Temperature': [parameters.get('airtemp', 'C')],
        'Humidity': [parameters.get('humidity', '%')],
        'J_magnitude': [m_j],
        'H_magnitude': [m_h],
        'K_magnitude': [m_k],
        'J_H_color': [color_jh],
        'H_K_color': [color_hk]
        'J_K_color': [color_jk]
    }

    df = pd.DataFrame(output_data)
    df.to_csv(args.outputfile, index=False)
    
if __name__ == "__main__":
    main()
