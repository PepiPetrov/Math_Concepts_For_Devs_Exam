import numpy as np
import matplotlib.pyplot as plt

from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as u
import astropy.constants as const

import ipyaladin.aladin_widget as ipyal


def read_spec(filename):
    '''Read a UVES spectrum from the ESO pipeline

    Parameters
    ----------
    filename : string
       name of the fits file with the data

    Returns
    -------
    wavelength : np.ndarray
        wavelength (in Ang)
    flux : np.ndarray
        flux (in erg/s/cm**2)
    date_obs : string
        time of observation
    '''
    sp = fits.open(filename)
    header = sp[0].header
    wcs = WCS(header)
    # make index array
    index = np.arange(header['NAXIS1'])

    wavelength = wcs.wcs_pix2world(index[:, np.newaxis], 0)
    wavelength = wavelength.flatten()
    flux = sp[0].data

    date_obs = header['Date-OBS']
    return wavelength, flux, date_obs


wavelength, flux, date_obs = read_spec('data/r.UVES.2011-08-11T234521.856-A01_0000.fits')

WIEN_DISPLACEMENT_CONSTANT = 2.9 * 10**7
STEFAN_BOLTZMAN_CONSTANT = 5.67 * 10 ** -8


def find_temperature(wavelength):
    '''
    Finds temperature of star by using its wavelength

    Args:
        wavelength: numpy.ndarray
    Returns:
        temperature: The temperature of the star (in Kelvin)
    '''

    temperature = WIEN_DISPLACEMENT_CONSTANT / np.max(wavelength)

    return temperature


find_temperature(wavelength)


def find_luminosity(wavelength, radius, temperature=None):
    '''
    Finds luminosity of star by using its wavelength

    Args:
        wavelength: numpy.ndarray
        radius: float - The radius of the star
        temperature: float | None - The temperature of the star, if none, it will be calculated
    Returns:
        luminosity: The luminosity of the star (in )
    '''

    if temperature is None:
        temperature = find_temperature(wavelength)

    star_surface = 4 * np.pi * radius ** 2
    luminosity_per_square_meter = STEFAN_BOLTZMAN_CONSTANT * temperature ** 4

    return star_surface * luminosity_per_square_meter

def calculate_redshift(velocity):
    '''
    Calculates redshift
    
    Args:
    velocity: The velocity (in Km/h)
    
    Returns:
    
    z: The redshift
    '''
    v_divided_by_speed_of_light = velocity / 1.07925285 * 10^9
    
    z =  v * (np.sqrt(1 - v_divided_by_speed_of_light) / np.sqrt(1 + v_divided_by_speed_of_light))