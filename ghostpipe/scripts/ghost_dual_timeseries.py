# -*- coding: utf-8 -*-
"""
Photometric flux-ratio time series from dual-IFU GHOST observations.

Given two matched sets of s1d spectra (e.g., target in IFU1 and a
comparison star in IFU2), this tool computes the time series of the flux
ratio between the two objects, both integrated and in wavelength channels,
producing a low-resolution differential (spectro-)photometric light curve.
A transit model computed with ``batman`` can be over-plotted for
comparison.

Usage example
-------------
::

    ghost_dual_timeseries --input1="*red*_TARGET_s1d.fits" \
                          --input2="*red*_COMPARISON_s1d.fits" -p

Created on Feb 25, 2024.

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

from optparse import OptionParser
import os,sys
import glob
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np

try:
    import batman
except ImportError:
    # Optional dependency, only required for the transit model overlay.
    # Install it with: pip install batman-package
    batman = None


def batman_model(time, per, t0, a, inc, rp, u0, u1=0., ecc=0., w=90.) :
    """Compute a transit model light curve with batman.

    Requires the optional dependency ``batman-package``.
    """
    if batman is None :
        raise ImportError("The transit model requires the optional dependency "
                          "'batman-package'. Install it with: pip install batman-package")
    
    """
        Function for computing transit models for the set of 8 free paramters
        x - time array
        """
    params = batman.TransitParams()
    
    params.per = per
    params.t0 = t0
    params.inc = inc
    params.a = a
    params.ecc = ecc
    params.w = w
    params.rp = rp
    params.u = [u0,u1]
    params.limb_dark = "quadratic"       #limb darkening model
    
    m = batman.TransitModel(params, time)    #initializes model
    
    flux_m = m.light_curve(params)          #calculates light curve
    
    return np.array(flux_m)


def wavelength_to_rgb(wavelength, gamma=0.8):
    ''' taken from http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    Additionally alpha value set to 0.5 outside range
    '''
    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 750:
        A = 1.
    else:
        A=0.5
    if wavelength < 380:
        wavelength = 380.
    if wavelength >750:
        wavelength = 750.
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R,G,B,A)



def main():
    """Entry point: parse command-line options and compute flux-ratio time series."""
    parser = OptionParser()
    parser.add_option("-1", "--input1", dest="input1", help="Input spectral s1d FITS data for object 1",type='string',default="")
    parser.add_option("-2", "--input2", dest="input2", help="Input spectral s1d FITS data for object 2 (comparison)",type='string',default="")
    parser.add_option("-p", action="store_true", dest="plot", help="plot", default=False)

    try:
        options,args = parser.parse_args(sys.argv[1:])
    except SystemExit:
        # let optparse handle --help and usage errors with its own exit code
        raise
    except:
        print("Error: check usage with -h ghost_dual_timeseries.py")
        sys.exit(1)

    print('Input Object 1 data pattern: ', options.input1)
    print('Input Object 2 data pattern: ', options.input2)

    ghost_dir = os.path.dirname(__file__) + '/'

    inputfiles1 = sorted(glob.glob(options.input1))
    inputfiles2 = sorted(glob.glob(options.input2))

    fluxratio = np.array([])
    bjd = np.array([])

    redchannel = True
    if "blue" in inputfiles1[0].lower() :
        redchannel = False

    calibrated = True

    #channel_size = 0.5
    channel_size = 5

    if redchannel :
        channels = np.arange(535,1050,channel_size)
        #channels = np.arange(650,800,channel_size)
    else :
        channels = np.arange(360,540,channel_size)
        #channels = np.arange(360,380,channel_size)

    print(channels)
    fluxratio_channels = []
    wlc = []
    for j in range(len(channels)-1) :
        fluxratio_channels.append(np.array([]))
        wlc.append((channels[j] + channels[j+1])/2)

    airmass=np.array([])

    for i in range(len(inputfiles1)) :
        print("Loading 1D spectra for {}/{} -> {} and {}".format(i+1,len(inputfiles1),os.path.basename(inputfiles1[i]),os.path.basename(inputfiles2[i])))

        hdu1, hdu2 = fits.open(inputfiles1[i]), fits.open(inputfiles2[i])
        hdr1, hdr2 = hdu1[0].header, hdu2[0].header

        wl1, wl2 = hdu1["WAVE"].data, hdu2["WAVE"].data
        flux1, fluxerr1 = hdu1["FLUX"].data, hdu1["FLUXERR"].data
        flux2, fluxerr2 = hdu2["FLUX"].data, hdu2["FLUXERR"].data

        fluxratio = np.append(fluxratio, np.nansum(flux1)/np.nansum(flux2))
        bjd = np.append(bjd, hdr1["BJD"])

        for j in range(len(channels)-1) :
            keep = (wl1 > channels[j]) & (wl1 < channels[j+1])
            f = np.nansum(flux1[keep])/np.nansum(flux2[keep])
            fluxratio_channels[j] = np.append(fluxratio_channels[j],f)

        airmass = np.append(airmass,hdr1["AIRMASS"])

    # Set initial values of transit parameters from TESS fit
    per = 2.6755504163
    t0 = 1597.0390602989 + 2457000
    a = 7.0714843363
    inc = 89.2782060751
    rp = 0.1116823681
    u0 = 0.3798792560
    u1 = 0.0676448537
    ecc = 0.0
    w = 90.0
    # Definir o array de tempos do modelo.
    t = np.arange(bjd[0],bjd[-1],0.001)
    # calcular o fluxo utilizando o biblioteca Batman.
    magmodel = 1+2.5*np.log10(batman_model(t, per, t0, a, inc, rp, u0, u1, ecc, w))


    obs_tr_model = 1+2.5*np.log10(batman_model(bjd, per, t0, a, inc, rp, u0, u1, ecc, w))

    median_magdiff = -2.5*np.log10(fluxratio)/np.nanmedian(-2.5*np.log10(fluxratio))

    calib = median_magdiff / obs_tr_model

    if not calibrated :
        plt.plot(bjd, median_magdiff, 'k-', lw=2, alpha=0.6, zorder=2, label="median")
    #plt.plot(bjd, calib, 'r--', lw=2, zorder=2, label="systematics")
    #plt.plot(bjd, airmass, '-', lw=2, zorder=2, label="airmass")

    sigma = np.array([])
    offset = 0.025

    for j in range(len(channels)-1) :
        magdiff = -2.5*np.log10(fluxratio_channels[j])/np.nanmedian(-2.5*np.log10(fluxratio_channels[j]))
        sig = np.nanstd(magdiff / calib - obs_tr_model)
        print("sigma = {:.4f} % or {:.4f} mmag".format(sig*100, sig*1000))

        sigma = np.append(sigma,sig)

        calib_mag = magdiff / calib + j*offset

        if wlc[-1] > 750 :
            color = wavelength_to_rgb(wlc[j] - (1030-750), gamma=0.8)
        else :
            color = wavelength_to_rgb(wlc[j], gamma=0.8)

        if calibrated :
            plt.plot(bjd, calib_mag, '.',color=color, zorder=1)
            plt.text(bjd[-1]+0.02, 1 + j*offset, r"$\lambda=${:.0f} nm".format(wlc[j]), verticalalignment='bottom', horizontalalignment='left', color=color, fontsize=15)
            plt.plot(t, magmodel+ j*offset, '-', color="k", lw=2)
        else :
            plt.plot(bjd, magdiff, '.',color=color, alpha=0.6, zorder=1)

    plt.xlabel("BJD",fontsize=16)
    plt.ylabel(r"$\Delta$mag",fontsize=16)
    #plt.legend(fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.show()

    if redchannel :
        plt.plot(wlc, sigma*100, "ro", label=r"Red channel, $\Delta\lambda={:.2f}$ nm".format(channel_size))
    else :
        plt.plot(wlc, sigma*100, "bo", label=r"Blue channel, $\Delta\lambda={:.2f}$ nm".format(channel_size))

    plt.ylabel(r"$\sigma$ (%)",fontsize=22)
    plt.xlabel(r"$\lambda$ (nm)",fontsize=22)
    plt.legend(fontsize=20)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.show()


if __name__ == "__main__":
    main()
