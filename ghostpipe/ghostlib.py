# -*- coding: utf-8 -*-
"""
Core utilities for handling Gemini/GHOST spectra.

This module provides the low-level building blocks used across the GHOST
Pipeline:

* ``red_order_limits`` / ``blue_order_limits``: adopted wavelength boundaries
  between adjacent echelle orders for the red and blue arms, used to stitch
  orders into a continuous 1D spectrum without overlap.
* ``write_1dspectrum_to_fits`` / ``write_spectrum_to_fits``: writers for the
  s1d FITS products (WAVE/FLUX/FLUXERR image extensions).
* ``read_template_product``: reader for template products created by
  ``ghost_template_s1d``.
* ``load_spectrum``: s1d spectrum loader with wavelength/flux filtering and
  a median-gradient outlier (spike) rejection.
* ``set_berv_bjd``: computes the Barycentric Julian Date (BJD) and the
  Barycentric Earth Radial Velocity (BERV) at Gemini South for a given
  observation and stores them in the FITS header.

Created on Nov 16, 2021.

Author
------
Eder Martioli <emartioli@lna.br>
Laboratório Nacional de Astrofísica (LNA/MCTI), Brazil
Institut d'Astrophysique de Paris, France

License
-------
This file is part of the GHOST Pipeline, distributed under the terms of
the GNU General Public License v3.0. See the LICENSE file at the root of
this project or <https://www.gnu.org/licenses/gpl-3.0.html> for details.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from copy import deepcopy

from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

#from reduc_lib import fit_continuum

def red_order_limits() :
    """Return the adopted order-boundary wavelengths for the RED arm.

    Returns
    -------
    list of float
        Wavelengths (nm), in decreasing order, marking the adopted
        boundaries between adjacent echelle orders of the GHOST red
        camera. These are used to stitch orders into a continuous 1D
        spectrum without overlap.
    """
    order_limits = [1040.88,
                    1026.09,
                    996.21,
                    968.16,
                    941.45,
                    916.33,
                    892.59,
                    870.16,
                    848.83,
                    828.68,
                    809.16,
                    790.61,
                    772.95,
                    756.17,
                    740.01,
                    724.44,
                    709.70,
                    695.35,
                    681.68,
                    668.48,
                    655.82,
                    643.57,
                    631.79,
                    620.43,
                    609.50,
                    598.88,
                    588.69,
                    578.81,
                    569.30,
                    560.04,
                    551.20,
                    542.57,
                    534.13]
    return order_limits

def blue_order_limits() :
    """Return the adopted order-boundary wavelengths for the BLUE arm.

    Returns
    -------
    list of float
        Wavelengths (nm), in decreasing order, marking the adopted
        boundaries between adjacent echelle orders of the GHOST blue
        camera. These are used to stitch orders into a continuous 1D
        spectrum without overlap.
    """
    order_limits = [543.80,
                    533.81,
                    525.89,
                    518.03,
                    510.29,
                    502.95,
                    495.68,
                    488.68,
                    481.82,
                    475.23,
                    468.75,
                    462.48,
                    456.37,
                    450.38,
                    444.57,
                    438.94,
                    433.41,
                    428.05,
                    422.86,
                    417.69,
                    412.71,
                    407.79,
                    403.06,
                    398.39,
                    393.83,
                    389.52,
                    385.24,
                    380.84,
                    376.65,
                    372.58,
                    368.71,
                    364.75,
                    360.96,
                    357.21]
    return order_limits

def write_1dspectrum_to_fits(filename, wave, flux, fluxerr, header=None) :
    """Save a 1D spectrum to a FITS file (BinTable format).

    Parameters
    ----------
    filename : str
        Output FITS file path (overwritten if it exists).
    wave : array-like
        Wavelength array (nm).
    flux : array-like
        Flux array.
    fluxerr : array-like
        Flux uncertainty array.
    header : astropy.io.fits.Header, optional
        Header to store in the primary HDU.
    """
    
    if header is None :
        header = fits.Header()
    
    header.set('TTYPE1', "WAVE")
    header.set('TUNIT1', "NM")
    header.set('TTYPE2', "FLUX")
    header.set('TUNIT2', "COUNTS")
    header.set('TTYPE3', "FLUXERR")
    header.set('TUNIT3', "COUNTS")

    primary_hdu = fits.PrimaryHDU(header=header)
    
    hdu_wl = fits.ImageHDU(data=wave, name="WAVE")
    hdu_flux = fits.ImageHDU(data=flux, name="FLUX")
    hdu_fluxerr = fits.ImageHDU(data=fluxerr, name="FLUXERR")

    listofhuds = [primary_hdu, hdu_wl, hdu_flux, hdu_fluxerr]

    mef_hdu = fits.HDUList(listofhuds)

    mef_hdu.writeto(filename, overwrite=True)


def read_template_product(filename) :
    """Read a template product created by ghost_template_s1d.

    Parameters
    ----------
    filename : str
        Input template FITS file path.

    Returns
    -------
    tuple
        ``(wave, flux, fluxerr, times, fluxes, fluxerrs, header)`` where
        ``wave``, ``flux`` and ``fluxerr`` describe the template spectrum,
        and ``times``, ``fluxes`` and ``fluxerrs`` contain the calibrated
        time series of spectra used to build the template.
    """
    hdulist = fits.open(filename)

    wave = hdulist["WAVE"].data
    flux = hdulist["FLUX"].data
    fluxerr = hdulist["FLUXERR"].data

    loc = {}
    loc["wl"] = wave
    loc["flux"] = flux
    loc["fluxerr"] = fluxerr

    return loc

def write_spectrum_to_fits(waves, fluxes, fluxerrs, filename, header=None):
    """Save a spectrum to a FITS file with WAVE/FLUX/FLUXERR image HDUs.

    Parameters
    ----------
    waves : array-like
        Wavelength array (nm).
    fluxes : array-like
        Flux array.
    fluxerrs : array-like
        Flux uncertainty array.
    filename : str
        Output FITS file path (overwritten if it exists).
    header : astropy.io.fits.Header, optional
        Header to store in the primary HDU.
    """
    
    if header is None :
        header = fits.Header()

    header.set('TTYPE1', "WAVE")
    header.set('TUNIT1', "NM")
    header.set('TTYPE2', "FLUX")
    header.set('TUNIT2', "COUNTS")
    header.set('TTYPE3', "FLUXERR")
    header.set('TUNIT3', "COUNTS")

    primary_hdu = fits.PrimaryHDU(header=header)
    hdu_wl = fits.ImageHDU(data=waves, name="WAVE")
    hdu_flux = fits.ImageHDU(data=fluxes, name="FLUX")
    hdu_err = fits.ImageHDU(data=fluxerrs, name="FLUXERR")
    mef_hdu = fits.HDUList([primary_hdu, hdu_wl, hdu_flux, hdu_err])

    mef_hdu.writeto(filename, overwrite=True)


def load_spectrum(filename, minwl=-1e20, maxwl=1e20, minflux=0, maxflux=1e30) :
    """Load an s1d spectrum from a FITS file, filtering out bad points.

    Points outside the given wavelength and flux ranges are rejected, as
    well as spectral spikes detected from the median absolute deviation of
    the flux-gradient differences between neighboring pixels.

    Parameters
    ----------
    filename : str
        Input s1d FITS file path (WAVE/FLUX/FLUXERR image HDUs).
    minwl, maxwl : float, optional
        Minimum and maximum wavelength (nm) to keep.
    minflux, maxflux : float, optional
        Minimum and maximum flux values to keep.

    Returns
    -------
    dict
        Dictionary with keys ``header``, ``wl``, ``flux`` and ``fluxerr``.
    """

    spectrum = {}
    
    hdu = fits.open(filename)
    hdr = hdu[0].header
    
    spectrum["header"] = deepcopy(hdr)

    wl = hdu["WAVE"].data
    
    threshold=100
    if wl[-1] > 700 :
        threshold=30
    
    flux = np.gradient(hdu["FLUX"].data) / wl
    #flux = hdu["FLUX"].data
    fluxerr = hdu["FLUXERR"].data

    tmp1before, tmp1after = np.full_like(hdu["FLUX"].data,np.nan), np.full_like(hdu["FLUX"].data,np.nan)
    tmp2before, tmp2after = np.full_like(hdu["FLUX"].data,np.nan), np.full_like(hdu["FLUX"].data,np.nan)
    
    # add array shifted by i pixels forward
    tmp1before[1:], tmp1after[:-1] = flux[:-1], flux[1:]
    tmp2before[2:], tmp2after[:-2] = flux[:-2], flux[2:]

    flux_array= []
    flux_array.append(np.abs(flux - tmp1before))
    flux_array.append(np.abs(flux - tmp2before))
    flux_array.append(np.abs(flux - tmp1after))
    flux_array.append(np.abs(flux - tmp2after))
    flux_array = np.array(flux_array)
    median_fluxes_diff = np.nanmedian(flux_array,axis=0)
 
    finite = np.isfinite(median_fluxes_diff)
    
    #continuum = np.full_like(median_fluxes_diff, np.nan)
    #continuum[finite] = fit_continuum(wl[finite], median_fluxes_diff[finite], function="polynomial", order=8, nit=10, rej_low=3., rej_high=1.0, grow=1, med_filt=1, percentile_low=0., percentile_high=100.,min_points=100, xlabel="wavelength", ylabel="flux", plot_fit=False, silent=True)
    #median_fluxes_diff /= continuum
    
    mad_flux_diff = stats.median_abs_deviation(median_fluxes_diff[finite])
    median_flux_diff = np.nanmedian(median_fluxes_diff)

    #print(median_flux_diff, mad_flux_diff)

    keep = (hdu["WAVE"].data > minwl) & (hdu["WAVE"].data < maxwl)
    keep &= np.isfinite(hdu["FLUX"].data)
    keep &= (hdu["FLUX"].data > minflux) & (hdu["FLUX"].data < maxflux)
    keep &= median_fluxes_diff < median_flux_diff + threshold * mad_flux_diff
    
    spectrum["wl"] = deepcopy(hdu["WAVE"].data[keep])
    spectrum["flux"] = deepcopy(hdu["FLUX"].data[keep])
    spectrum["fluxerr"] = deepcopy(hdu["FLUXERR"].data[keep])
    
    
    """
    
    plt.plot(hdu["WAVE"].data[keep],hdu["FLUX"].data[keep]/np.nanmedian(hdu["FLUX"].data[keep]),'ko')
    plt.plot(hdu["WAVE"].data[~keep],hdu["FLUX"].data[~keep]/np.nanmedian(hdu["FLUX"].data[keep]),'ro')
    
    #plt.plot(hdu["WAVE"].data[keep],flux[keep]/np.nanmedian(flux[keep]),'ko')
    #plt.plot(hdu["WAVE"].data[~keep],flux[~keep]/np.nanmedian(flux[keep]),'ro')
    
    norm_factor = np.nanmedian(median_fluxes_diff)
    
    plt.plot(hdu["WAVE"].data,median_fluxes_diff/norm_factor,'g-')
    
    ymedian = np.full_like(hdu["FLUX"].data,median_flux_diff)
    ythreshold = np.full_like(hdu["FLUX"].data,median_flux_diff + threshold * mad_flux_diff)
    
    plt.plot(hdu["WAVE"].data,ymedian/norm_factor,'b--')
    plt.plot(hdu["WAVE"].data,ythreshold/norm_factor,'b-')

    plt.show()
    exit()

    """
    hdu.close()
    
    return spectrum


def set_berv_bjd(hdr, use_jpl_ephemeris=False) :
    """Compute BJD and BERV at Gemini South and store them in the header.

    The Barycentric Julian Date (BJD, TDB) at mid-exposure and the
    Barycentric Earth Radial Velocity (BERV) toward the target coordinates
    are computed with astropy and added to the header under the keywords
    ``BJD`` and ``BERV``, along with ``JD``, ``MJD``, ``HJD``, the target
    coordinates in degrees (``RA_DEG``/``DEC_DEG``), the observatory
    geographic coordinates, and the ISOT date at mid-exposure.

    Parameters
    ----------
    hdr : astropy.io.fits.Header
        FITS header containing the target coordinates (``RA``/``DEC``),
        observation date (``DATE-OBS``/``UTSTART``) and exposure time.
    use_jpl_ephemeris : bool, optional
        If True, use the JPL DE430 ephemeris (requires download); if
        False, use astropy's built-in ephemeris.

    Returns
    -------
    astropy.io.fits.Header
        The updated header, with ``BJD`` and ``BERV`` (km/s) keywords set.
    """

    # set Gemini South geographic coordinates
    longitude = -70.7353463919
    latitude = -30.2377823822
    altitude = 2722*u.m #hdr['OBSALT']

    hdr.set("OBSLONG",longitude,"East longitude")
    hdr.set("OBSLAT",longitude,"North latitude")
    hdr.set("OBSALT",2722,"elevation above sea level [m]")

    observatory_location = EarthLocation.from_geodetic(lat=latitude, lon=longitude, height=altitude)

    # set equinox to 2000 or get it from header if it exists
    equinox="J{:.1f}".format(hdr["EQUINOX"])

    # set source observed
    source = SkyCoord(hdr['RA'], hdr['DEC'], unit=(u.hourangle, u.deg), frame='icrs', equinox=equinox)

    radeg, decdeg = source.ra.value, source.dec.value

    halfexptime = TimeDelta((hdr["EXPTIME"]/2.) * u.s)
    
    timestr = "{}T{}".format(hdr["DATE-OBS"],hdr["UTSTART"])

    obstime = Time(timestr, format='isot', scale='utc', location=observatory_location) + halfexptime

    jd = obstime.jd
    mjd = obstime.mjd
    
    # Set light travel time for source observed
    if use_jpl_ephemeris :
        # below is for more precise time
        ltt_bary = obstime.light_travel_time(source, ephemeris='jpl')
    else :
        ltt_bary = obstime.light_travel_time(source)
        
    bjd = obstime.tdb.jd + ltt_bary
    
    barycorr = source.radial_velocity_correction(obstime=obstime)
    berv = barycorr.to(u.km/u.s).value
    #### HJD
    ltt_helio = obstime.light_travel_time(source, 'heliocentric') ### para o HJD
    hjd = obstime.tdb.jd + ltt_helio

    hdr.set("JD",jd,"Julian date at mid of exposure")
    hdr.set("MJD",mjd,"Modified Julian date at mid of exposure")
    hdr.set("BJD",bjd.value,"Barycentric Julian date at mid of exposure")
    hdr.set("HJD",hjd.value,"Heliocentric Julian date at mid of exposure")
    hdr.set("BERV",berv,"Barycentric Earth radial velocity [km/s]")
    hdr.set("RA_DEG",radeg,"Target right ascension in degree")
    hdr.set("DEC_DEG",decdeg,"Target declination in degree")
    hdr.set("ISOTDATE",obstime.isot,"UT data at mid of exposure [isot]")

    #sidereal = obstime.sidereal_time('apparent')
    #hdr.set("ST",sidereal,"Sidereal time")
    #hdr.set("SD",sidereal,"Sidereal time")
    
    # calculate airmass
    #airmass = source.transform_to(AltAz(obstime=obstime,location=observatory_location)).secz
    #hdr.set("AIRMASS",airmass.value,"Airmass at start of exposure")

    return hdr
