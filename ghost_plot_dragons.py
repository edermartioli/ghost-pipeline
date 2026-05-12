# -*- coding: iso-8859-1 -*-
"""
    Created on Dec 1 2025
    
    Description: This routine plot the GHOST spectra in DRAGONS format
    
    @author: Eder Martioli <emartioli@lna.br>
    Laboratório Nacional de Astrofísica - LNA/MCTI
    
    python ghost_plot_dragons.py --blue=/Users/eder/Desktop/q-204/S20251128S0087_blue001_dragons.fits --red=/Users/eder/Desktop/q-204/S20251128S0087_red001_dragons.fits
    
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-b", "--blue", dest="blue", help="Input blue spectrum",type='string',default="")
parser.add_option("-r", "--red", dest="red", help="Input red spectrum",type='string',default="")
try:
    options,args = parser.parse_args(sys.argv[1:])
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

