# -*- coding: utf-8 -*-
"""
Quick-look plot of GHOST spectra in DRAGONS format.

Plots the order-by-order flux of a blue and/or red GHOST spectrum reduced
with DRAGONS ("dragons" format, one array per echelle order).

Usage example
-------------
::

    ghost_plot_dragons --blue=S20251128S0087_blue001_dragons.fits \
                       --red=S20251128S0087_red001_dragons.fits

Created on Dec 1, 2025.

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
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt


def main():
    """Entry point: parse command-line options and plot the spectra."""
    parser = OptionParser()
    parser.add_option("-b", "--blue", dest="blue", help="Input blue spectrum",type='string',default="")
    parser.add_option("-r", "--red", dest="red", help="Input red spectrum",type='string',default="")
    try:
        options,args = parser.parse_args(sys.argv[1:])
    except SystemExit:
        # let optparse handle --help and usage errors with its own exit code
        raise
    except:
        print("Error: check usage with -h ghost_plot_dragons.py")
        sys.exit(1)

    print('Plotting GHOST spectra: ')
    print('\t Blue: ', options.blue)
    print('\t Red: ', options.red)

    ext = 1

    blue = fits.open(options.blue)
    red = fits.open(options.red)

    wl_blue = blue[ext].header['CRVAL1'] + blue[ext].header['CD1_1']*(np.arange(len(blue[ext].data))-blue[ext].header['CRPIX1'])
    wl_red = red[ext].header['CRVAL1'] + red[ext].header['CD1_1']*(np.arange(len(red[ext].data))-red[ext].header['CRPIX1'])

    plt.plot(wl_blue, blue[ext].data, "-", color="darkblue", label="Blue channel")
    plt.plot(wl_red, red[ext].data, "-", color="darkred", label="Red channel")

    plt.tick_params(axis='x', labelsize=10)
    plt.tick_params(axis='y', labelsize=10)
    plt.minorticks_on()
    plt.tick_params(which='minor', length=3, width=0.7, direction='in',bottom=True, top=True, left=True, right=True)
    plt.tick_params(which='major', length=7, width=1.2, direction='in',bottom=True, top=True, left=True, right=True)
    plt.xlabel(r"wavelength [nm]",fontsize=20)
    plt.ylabel(r"flux",fontsize=20)
    plt.legend(fontsize=18)
    plt.show()


if __name__ == "__main__":
    main()
