# -*- coding: iso-8859-1 -*-
"""
    Created on May 30 2023
    
    Description: Get simple template
    
    @author: Eder Martioli <emartioli@lna.br>
    
    Laboratório Nacional de Astrofísica, Brazil.
    Institut d'Astrophysique de Paris, France.
    
    Simple usage example:
    
    python ghost_get_simple_template.py --inputBlue=/Users/eder/Science/EtaCarinae/ghost_templates/etacar_ghost_blue.fits --inputRed=/Users/eder/Science/EtaCarinae/ghost_templates/etacar_ghost_red.fits -p

    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
import glob

import astropy.io.fits as fits
import matplotlib.pyplot as plt
import ghostlib

parser = OptionParser()
parser.add_option("-b", "--inputBlue", dest="inputBlue", help="Input Blue template fits file",type='string',default="")
parser.add_option("-r", "--inputRed", dest="inputRed", help="Input Red template fits file",type='string',default="")
parser.add_option("-o", "--output", dest="output", help="Output simple template fits file",type='string',default="")
parser.add_option("-p", action="store_true", dest="plot", help="plot", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
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
    ghostlib.write_spectrum_to_fits(wave, flux, fluxerr, options.output, header=hdr)

