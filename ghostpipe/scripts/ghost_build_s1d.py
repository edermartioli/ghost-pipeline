# -*- coding: utf-8 -*-
"""
Build s1d (1D spectrum) products from GHOST spectra reduced with DRAGONS.

This tool stitches the echelle orders of ``*_calibrated.fits`` products
into a continuous 1D spectrum per arm (red/blue), using the adopted order
boundaries in :mod:`ghostpipe.ghostlib`. It also computes the BJD and BERV
and stores them in the header of the output s1d product. Fluxes from IFU2
and from the sky fiber can be included when available.

Usage examples
--------------
For data reduced with DRAGONS with sky subtraction turned off, when the
object was observed also with IFU2::

    ghost_build_s1d --input="/path/to/*_calibrated.fits" -2

For single-IFU observations::

    ghost_build_s1d --input="/path/to/*_calibrated.fits"

Created on May 26, 2023.

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
from ghostpipe import ghostlib
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from scipy import constants


def main():
    """Entry point: parse command-line options and build s1d products."""
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input", help="Input spectral FITS data pattern",type='string',default="*_calibrated.fits")
    parser.add_option("-2", action="store_true", dest="hasifu2", help="has IFU2", default=False)
    parser.add_option("-b", action="store_true", dest="applyberv", help="apply barycentric correction", default=False)
    parser.add_option("-s", action="store_true", dest="hassky", help="has sky", default=False)
    parser.add_option("-p", action="store_true", dest="plot", help="plot", default=False)

    try:
        options,args = parser.parse_args(sys.argv[1:])
    except SystemExit:
        # let optparse handle --help and usage errors with its own exit code
        raise
    except:
        print("Error: check usage with -h ghost_build_s1d.py")
        sys.exit(1)

    print('Input GHOST data pattern: ', options.input)

    ghost_dir = os.path.dirname(__file__) + '/'

    inputfiles = sorted(glob.glob(options.input))

    for i in range(len(inputfiles)) :

        print("Generating 1D spectra for {}/{} -> {}".format(i+1,len(inputfiles),inputfiles[i]))
        hdul = fits.open(inputfiles[i])
        hdr = hdul[0].header

        hdr = ghostlib.set_berv_bjd(hdr, use_jpl_ephemeris=True)

        lower_limit, upper_limit = 0,1e30
        if hdr['CAMERA'] == 'RED' :
            lower_limit = 530.0
            norders = 33
            wl_limits = ghostlib.red_order_limits()

        elif hdr['CAMERA'] == 'BLUE' :
            upper_limit = 533.9
            wl_limits = ghostlib.blue_order_limits()

        norders = len(wl_limits)

        wl = np.array([],dtype=float)
        flux1, fluxerr1 = np.array([],dtype=float),np.array([],dtype=float)
        flux2, fluxerr2 = np.array([],dtype=float),np.array([],dtype=float)
        skyflux, skyfluxerr = np.array([],dtype=float),np.array([],dtype=float)

        for order in range(norders) :

            #wave = np.array(hdul['WAVL'].data[order]/10,dtype=float)
            wave = np.array(hdul['AWAV'].data[order],dtype=float)
            keep = (wave > lower_limit) & (wave < upper_limit)

            if order == norders-1 :
                keep &= wave < wl_limits[order]
            else :
                keep &= (wave >= wl_limits[order+1]) & (wave < wl_limits[order])

            #print(order, wave[0], wave[-1],"->", wl_limits[order])

            wl = np.append(wl,wave[keep])
            #print("processing order:{} -> array size:{}".format(order, len(wave)))

            #flux = np.array(hdul[1].data[order,:,0],dtype=float)
            flux = np.array(hdul[1].data[order,:],dtype=float)
            flux1 = np.append(flux1,flux[keep])
            #fluxerr1 = np.append(fluxerr1,np.sqrt(hdul[2].data[order,:,0][keep]))
            fluxerr1 = np.append(fluxerr1,np.sqrt(hdul[2].data[order,:][keep]))
            #plt.plot(wave,flux/np.median(flux))
            #plt.plot(wave,np.sqrt(hdul[2].data[order,:,0]))

            if options.hasifu2 :
                #flux = np.array(hdul[1].data[order,:,1],dtype=float)
                flux = np.array(hdul[5].data[order,:],dtype=float)
                flux2 = np.append(flux2,flux[keep])
                #fluxerr2 = np.append(fluxerr2,np.sqrt(hdul[2].data[order,:,1][keep]))
                fluxerr2 = np.append(fluxerr2,np.sqrt(hdul[6].data[order,:][keep]))

            #plt.plot(wave,0.5 + flux/np.median(flux))
            if options.hassky :
                #flux = np.array(hdul[1].data[order,:,2],dtype=float)
                flux = np.array(hdul[5].data[order,:],dtype=float)
                skyflux = np.append(skyflux,flux[keep])
                #skyfluxerr = np.append(skyfluxerr,np.sqrt(hdul[2].data[order,:,2][keep]))
                skyfluxerr = np.append(skyfluxerr,np.sqrt(hdul[6].data[order,:][keep]))

        #plt.show()
        sort = np.argsort(wl)
        #print("Final array size: {}".format(len(wl[sort])))

        outwl = deepcopy(wl[sort])
        if options.applyberv :
            speed_of_light_in_kps = constants.c / 1000.
            vel_shift = hdr['BERV']
            outwl = outwl / (1.0 + vel_shift / speed_of_light_in_kps)

        object1 = hdr["OBJECT"].replace(" ","")
        if hdr["RESOLUT"] == "Standard" :
            object1 = hdr["SRIFU1"].replace(" ","")

        output1 = inputfiles[i].replace("_calibrated.fits","_{}_s1d.fits".format(object1))
        print("Saving 1D spectrum for IFU1 object: {} to file -> {}".format(hdr["SRIFU1"],output1))
        hdr["OBJECT"] = object1
        ghostlib.write_1dspectrum_to_fits(output1, outwl, flux1[sort], fluxerr1[sort], header=hdr)
        if options.plot :
            plt.plot(outwl, flux1[sort])

        if options.hasifu2 :
            object2 = hdr["SRIFU2"].replace(" ","")
            output2 = inputfiles[i].replace("_calibrated.fits","_{}_s1d.fits".format(object2))
            print("Saving 1D spectrum for IFU2 object: {} to file -> {}".format(hdr["SRIFU2"],output2))
            hdr["OBJECT"] = object2
            ghostlib.write_1dspectrum_to_fits(output2, outwl, flux2[sort], fluxerr2[sort], header=hdr)
            if options.plot :
                plt.plot(outwl, flux2[sort])

        if options.hassky :
            sky_output = inputfiles[i].replace("_calibrated.fits","_sky_s1d.fits")
            print("Saving 1D spectrum for Sky to file -> {}".format(sky_output))
            hdr["OBJECT"] = "Sky"
            ghostlib.write_1dspectrum_to_fits(sky_output, outwl, skyflux[sort], skyfluxerr[sort], header=hdr)
            if options.plot :
                plt.plot(outwl, skyflux[sort])

        if options.plot :
            plt.show()


if __name__ == "__main__":
    main()
