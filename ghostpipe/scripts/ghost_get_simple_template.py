# -*- coding: utf-8 -*-
"""
Inspect blue+red GHOST template products and export a simple template.

This tool loads a blue and a red template FITS product (as produced by
``ghost_template_s1d``), plots the template spectra along with the
individual spectra in the time series, and optionally combines the two
arms into a single simple s1d spectrum saved to an output FITS file.

Usage example
-------------
::

    ghost_get_simple_template --inputBlue=target_ghost_blue_template.fits \
        --inputRed=target_ghost_red_template.fits -p

Created on May 30, 2023.

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
import sys

import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from ghostpipe import ghostlib


def main():
    """Entry point: parse command-line options, plot and export the template."""
    parser = OptionParser()
    parser.add_option("-b", "--inputBlue", dest="inputBlue", help="Input Blue template fits file",type='string',default="")
    parser.add_option("-r", "--inputRed", dest="inputRed", help="Input Red template fits file",type='string',default="")
    parser.add_option("-o", "--output", dest="output", help="Output simple template fits file",type='string',default="")
    parser.add_option("-p", action="store_true", dest="plot", help="plot", default=False)
    parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

    try:
        options,args = parser.parse_args(sys.argv[1:])
    except SystemExit:
        # let optparse handle --help and usage errors with its own exit code
        raise
    except:
        print("Error: check usage with  -h ghost_get_simple_template.py")
        sys.exit(1)

    if options.verbose:
        print('Input Blue template fits file: ', options.inputBlue)
        print('Input Red template fits file: ', options.inputRed)
        print('Output simple template fits product: ', options.output)

    plot_snr = False
    max_nspc = 5

    hdulBlue = fits.open(options.inputBlue)
    hdrBlue = hdulBlue[0].header

    hdulRed = fits.open(options.inputRed)
    hdrRed = hdulRed[0].header

    print(repr(hdrBlue))
    print(repr(hdrRed))


    waveBlue, waveRed = hdulBlue['WAVE'].data, hdulRed['WAVE'].data
    fluxBlue, fluxRed = hdulBlue['TEMPLATE_FLUX'].data, hdulRed['TEMPLATE_FLUX'].data
    fluxerrBlue, fluxerrRed = hdulBlue['TEMPLATE_FLUXERR'].data, hdulRed['TEMPLATE_FLUXERR'].data

    Bluefluxes, Redfluxes = hdulBlue['FLUXES'].data, hdulRed['FLUXES'].data
    Bluefluxerrs, Redfluxerrs = hdulBlue['FLUXERRS'].data, hdulRed['FLUXERRS'].data

    if options.plot :
        nspcBlue, nspcRed = len(Bluefluxes), len(Redfluxes)
        print("Nblue={} NRed={}".format(nspcBlue, nspcRed))
        if nspcBlue > max_nspc :
            nspcBlue = max_nspc
        if nspcRed > max_nspc :
            nspcRed = max_nspc

        for i in range(nspcBlue) :
            if not plot_snr :
                plt.errorbar(waveBlue, Bluefluxes[i], yerr=Bluefluxerrs[i], fmt='.', color='darkblue', alpha=0.3, zorder=1)

        for i in range(nspcRed) :
            if not plot_snr :
                plt.errorbar(waveRed, Redfluxes[i], yerr=Redfluxerrs[i], fmt='.', color='darkred', alpha=0.3, zorder=1)

        if plot_snr :
            plt.plot(waveBlue, fluxBlue/fluxerrBlue, 'b-', alpha=0.8, lw=2, zorder=2)
            plt.plot(waveRed, fluxRed/fluxerrRed, 'r-', alpha=0.8, lw=2, zorder=2)
        else :
            plt.plot(waveBlue, fluxBlue, 'b-', alpha=0.8, lw=2, zorder=2)
            plt.plot(waveRed, fluxRed, 'r-', alpha=0.8, lw=2, zorder=2)

        plt.tick_params(axis='x', labelsize=14)
        plt.tick_params(axis='y', labelsize=14)
        plt.minorticks_on()
        plt.tick_params(which='minor', length=3, width=0.7, direction='in',bottom=True, top=True, left=True, right=True)
        plt.tick_params(which='major', length=7, width=1.2, direction='in',bottom=True, top=True, left=True, right=True)
        plt.xlabel(r"wavelength [nm]",fontsize=20)
        if plot_snr :
            plt.ylabel(r"SNR",fontsize=20)
        else :
            plt.ylabel(r"flux",fontsize=20)
        plt.show()


    if options.output != "" :
        # Combine the blue and red templates into a single simple spectrum,
        # sorted by wavelength, and save it as an s1d FITS product.
        wave = np.concatenate((waveBlue, waveRed))
        flux = np.concatenate((fluxBlue, fluxRed))
        fluxerr = np.concatenate((fluxerrBlue, fluxerrRed))
        sortmask = np.argsort(wave)
        ghostlib.write_spectrum_to_fits(wave[sortmask], flux[sortmask], fluxerr[sortmask], options.output, header=hdrBlue)


if __name__ == "__main__":
    main()
