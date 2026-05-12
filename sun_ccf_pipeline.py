# -*- coding: iso-8859-1 -*-
"""
    Created on Feb 20 2023
    
    Description: Calculate CCF for a set of Solar IAG Atlas spectra
    
    @author: Eder Martioli <emartioli@lna.br>
    
    Laboratório Nacional de Astrofísica, Brazil.
    Institut d'Astrophysique de Paris, France.
    
    Simple usage example:
    
    python /Users/eder/ghost-pipeline/sun_ccf_pipeline.py --ccf_mask=/Users/eder/ghost-pipeline/masks/G2_nm.mas --input=/Users/eder/Data/IAGSolarAtlas/solarspectrum_mu*.fits --source_rv=83.74638 --output_ccfs_file="sun_ccfs.fits" -pv

    python /Users/eder/ghost-pipeline/sun_ccf_pipeline.py --ccf_mask=/Users/eder/Science/WASP-108/SYNTHETIC_SPECTRUM/WASP-108_excl_tel_ranges.mas --input=/Users/eder/Data/IAGSolarAtlas/solarspectrum_mu*.fits --source_rv=83.74638 --output_ccfs_file="sun_ccfs.fits" -pv



    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
import glob

import matplotlib.pyplot as plt
import ghostlib
import reduc_lib
import ccf_lib
from copy import deepcopy

import numpy as np
import spectrallib
from scipy import constants

import astropy.io.fits as fits

ghost_ccf_dir = os.path.dirname(__file__)


def load_array_of_solar_spectra(inputdata, wl0=0, wlf=1e30, verbose=False) :

    loc = {}
    loc["input"] = inputdata

    spectra = []

    for i in range(len(inputdata)) :
        
        hdul = fits.open(inputdata[i])
        hdr = hdul[0].header

        spectrum = {}
        
        # set source RVs
        spectrum['FILENAME'] = inputdata[i]
        
        spectrum['OBJECT'] = "Sun"
        spectrum['MU'] = hdr["MU"]
        spectrum['EXPTIME'] = hdr["EXPTIME"]
        spectrum['NOBS'] = hdr["NR_OBS"]

        keep = (hdul[0].data[0]>wl0*10) & (hdul[0].data[0]<wlf*10)

        spectrum["wl"] = hdul[0].data[0][keep] / 10.
        spectrum["flux"] = hdul[0].data[1][keep]
        spectrum["fluxerr"] = spectrum["flux"] * 0.001

        spectrum["weights"] = np.ones_like(spectrum["wl"])

        spectrum['SNR'] = np.nanmedian(spectrum["flux"]) / np.nanmedian(spectrum["fluxerr"])
        hdr['OBJECT'] = "Sun"
        hdr['SNR'] = spectrum['SNR']
        hdr['BJD'] = 0.
        
        if verbose :
            print("Spectrum ({0}/{1}): {2} OBJ={3} MU={4} EXPTIME={5} NOBS={6}".format(i+1,len(inputdata),os.path.basename(inputdata[i]),spectrum['OBJECT'],spectrum['MU'],spectrum['EXPTIME'],spectrum['NOBS']))
        
        spectrum['header'] = hdr
        spectra.append(spectrum)
    
    loc["spectra"] = spectra

    return loc


def save_ccfs_to_fits(template_ccf, mus, filename, header=None) :

    rv = template_ccf["wl"]
    tccf = template_ccf["flux"]
    ccfs = template_ccf["flux_arr_sub"] * template_ccf["flux"]
    residuals = template_ccf["flux_arr_sub"]
    nspc = len(template_ccf["flux_arr_sub"])

    if header is None :
        header = fits.Header()

    header.set('TTYPE1', "RV")
    header.set('TUNIT1', "KPS")
    header.set('TTYPE2', "MEANCCF")
    header.set('TUNIT2', "COUNTS")
    header.set('TTYPE3', "CCFS")
    header.set('TUNIT3', "COUNTS")
    header.set('TTYPE4', "RESIDUAL_CCFS")
    header.set('TUNIT4', "COUNTS")

    primary_hdu = fits.PrimaryHDU(header=header)
    hdu_rv = fits.ImageHDU(data=rv, name="RV")
    hdu_mccf = fits.ImageHDU(data=tccf, name="MEANCCF")
    hdu_ccfs = fits.ImageHDU(data=ccfs, name="CCFS")
    hdu_resccfs = fits.ImageHDU(data=residuals, name="RESIDUAL_CCFS")
    hdu_mus = fits.ImageHDU(data=mus, name="MU")

    mef_hdu = fits.HDUList([primary_hdu, hdu_rv, hdu_mccf, hdu_ccfs, hdu_resccfs, hdu_mus])

    mef_hdu.writeto(filename, overwrite=True)




parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Spectral *s1d_A.fits data pattern",type='string',default="*.fits")
parser.add_option("-m", "--ccf_mask", dest="ccf_mask", help="Input CCF mask",type='string',default="")
parser.add_option("-c", "--output_ccfs_file", dest="output_ccfs_file", help="Output CCFs file (fits format)",type='string',default="")
parser.add_option("-r", "--source_rv", dest="source_rv", help="Input source RV (km/s)",type='float',default=0.)
parser.add_option("-w", "--ccf_width", dest="ccf_width", help="CCF half width (km/s)",type='string',default="60")
parser.add_option("-s", action="store_true", dest="saveccfs", help="Save CCF to FITS files", default=False)
parser.add_option("-p", action="store_true", dest="plot", help="plot", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with  -h sun_ccf_pipeline.py")
    sys.exit(1)

if options.verbose:
    print('Spectral s1d.fits data pattern: ', options.input)
    if options.ccf_mask != "":
        print('Input CCF mask: ', options.ccf_mask)
    print('Initial CCF width (km/s): ', options.ccf_width)

# make list of tfits data files
if options.verbose:
    print("Creating list of s1d fits spectrum files...")
inputdata = sorted(glob.glob(options.input))

verbose = options.verbose

# First load spectra into a container
#array_of_spectra = load_array_of_solar_spectra(inputdata, wl0=400, wlf=700, verbose=verbose)
array_of_spectra = load_array_of_solar_spectra(inputdata, verbose=verbose)

ccf_width = options.ccf_width
ccf_mask = options.ccf_mask

# Start dealing with CCF related parameters and construction of a weighted mask
# load science CCF parameters
ccf_params = ccf_lib.set_ccf_params(ccf_mask)
# update ccf width with input value
ccf_params["CCF_WIDTH"] = float(ccf_width)

index = 12
spectrum = array_of_spectra["spectra"][index]
header = spectrum["header"]

ccfmask = ccf_lib.apply_weights_to_ccf_mask(ccf_params, spectrum["wl"], spectrum["flux"], spectrum["fluxerr"], spectrum["weights"], median=True, remove_lines_with_nans=True, source_rv=options.source_rv, verbose=False, plot=False)

sun_ccf = ccf_lib.run_ccf_eder(ccf_params, spectrum["wl"], spectrum["flux"], header, ccfmask, targetrv=options.source_rv, normalize_ccfs=True, plot=options.plot, verbose=options.verbose)

source_rv = sun_ccf["header"]['RV_OBJ']
ccf_params["SOURCE_RV"] = source_rv
ccf_params["CCF_WIDTH"] = 7 * sun_ccf["header"]['CCFMFWHM']

if verbose :
    print("RV={:.5f} km/s CCF FWHM={:.2f} km/s CCF window size={:.2f} km/s".format(sun_ccf["header"]['RV_OBJ'], sun_ccf["header"]['CCFMFWHM'], sun_ccf["CCF_WIDTH"]))


spectra = array_of_spectra['spectra']

mean_fwhm, calib_rv = [], []
rvccf_data, ccf_data = [], []
mu = []

for i in range(len(spectra)) :

    spectrum = spectra[i]
    basename = os.path.basename(spectrum['FILENAME'])
    
    if verbose :
        print("Running CCF on file {0}/{1} -> {2}".format(i,len(spectra)-1,basename))

    fluxes, waves_sf = spectrum['flux'], spectrum['wl']

    # run main routine to process ccf on science fiber
    header = spectrum["header"]

    # run an adpated version of the ccf codes using reduced spectra as input
    sci_ccf = ccf_lib.run_ccf_eder(ccf_params, waves_sf, fluxes, header, ccfmask, filename=spectrum['FILENAME'], targetrv=ccf_params["SOURCE_RV"], normalize_ccfs=True, output="", plot=False, verbose=False)

    mean_fwhm.append(sci_ccf["header"]['CCFMFWHM'])
    calib_rv.append(sci_ccf["header"]['RV_OBJ'] )
    mu.append(sci_ccf["header"]['MU'])

    if verbose :
        print("Spectrum: {0} DATE={1} Sci_RV={2:.5f} km/s".format(basename, sci_ccf["header"]["DATE"], sci_ccf["header"]['RV_OBJ']))
            
    rvccf_data.append(sci_ccf['RV_CCF'])
    ccf_data.append(sci_ccf['MEAN_CCF'])

mean_fwhm = np.array(mean_fwhm)
velocity_window = 1.5*np.nanmedian(mean_fwhm)

obj = sci_ccf["header"]["OBJECT"].replace(" ","")

# The solution is to work on a simpler version of ccf2rv routines, but it takes a bit of time for coding and testing. To do soon.
rvs, rverrs, template_ccf = ccf_lib.ccf_analysis(rvccf_data[0], np.array(ccf_data,dtype=float), calib_rv, nsig_clip=0, velocity_window=velocity_window, plot=True, verbose=True)

if options.output_ccfs_file != "" :
    save_ccfs_to_fits(template_ccf, mu, options.output_ccfs_file, header=header)

if options.plot :
    ccf_lib.plot_ccfs(template_ccf)

    plt.errorbar(mu, rvs, yerr=rverrs, fmt="o")
    plt.xlabel(r"$\mu$")
    plt.ylabel(r"Velocity [m/s]")
    plt.show()
