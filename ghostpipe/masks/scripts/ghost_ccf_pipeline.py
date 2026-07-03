# -*- coding: utf-8 -*-
"""
Cross-correlation (CCF) radial-velocity analysis of GHOST s1d spectra.

This is the main radial-velocity tool of the GHOST Pipeline. It loads a
time series of s1d spectra, builds a template, computes the CCF of each
spectrum with a weighted line mask, and measures radial velocities by a
template-matching analysis of the CCFs. It can also produce the CCF FITS
products and the time series of bisector span and FWHM.

Output time series are saved in rdb format, and CCFs in FITS format.

Usage examples
--------------
::

    ghost_ccf_pipeline --ccf_mask=masks/G2_nm.mas --nknots=40 \
        --input="*red*_s1d.fits" \
        --output_ccfs_file=target_ghost_red_ccfs.fits \
        --output_rv_file=target_ghost_red_ccfrv.rdb \
        --output_bis_file=target_ghost_red_ccfbis.rdb \
        --output_fwhm_file=target_ghost_red_ccffwhm.rdb \
        --output_obslog_file=target_ghost_red_log.txt -pv

    ghost_ccf_pipeline --ccf_mask=masks/G2_nm.mas --nknots=10 \
        --input="*blue*_s1d.fits" \
        --output_rv_file=target_ghost_blue_ccfrv.rdb -pv

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

import matplotlib.pyplot as plt
from ghostpipe import ghostlib
from ghostpipe import reduc_lib
from ghostpipe import ccf_lib
from copy import deepcopy

import numpy as np
from ghostpipe import spectrallib
from scipy import constants

ghost_ccf_dir = os.path.dirname(__file__)


def main():
    """Entry point: parse command-line options and run the CCF RV analysis."""
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input", help="Spectral *s1d_A.fits data pattern",type='string',default="*.fits")
    parser.add_option("-m", "--ccf_mask", dest="ccf_mask", help="Input CCF mask",type='string',default="")
    parser.add_option("-o", "--output_template", dest="output_template", help="Output template spectrum",type='string',default="")
    parser.add_option("-c", "--output_ccfs_file", dest="output_ccfs_file", help="Output CCFs file (fits format)",type='string',default="")
    parser.add_option("-t", "--output_rv_file", dest="output_rv_file", help="Output RV file (rdb format)",type='string',default="")
    parser.add_option("-b", "--output_bis_file", dest="output_bis_file", help="Output bisector file (rdb format)",type='string',default="")
    parser.add_option("-f", "--output_fwhm_file", dest="output_fwhm_file", help="Output FWHM file (rdb format)",type='string',default="")
    parser.add_option("-l", "--output_obslog_file", dest="output_obslog_file", help="Output observation log file",type='string',default="")
    parser.add_option("-z", "--rvfile", dest="rvfile", help="Input RV file",type='string',default="")
    parser.add_option("-r", "--source_rv", dest="source_rv", help="Input source RV (km/s)",type='float',default=0.)
    parser.add_option("-w", "--ccf_width", dest="ccf_width", help="CCF half width (km/s)",type='string',default="150")
    parser.add_option("-a", "--vel_sampling", dest="vel_sampling", help="Velocity sampling for the template spectrum (km/s)",type='float',default=1.8)
    parser.add_option("-k", "--nknots", dest="nknots", help="Number of knots for normalization",type='int',default=10)
    parser.add_option("-s", action="store_true", dest="saveccfs", help="Save CCF to FITS files", default=False)
    parser.add_option("-p", action="store_true", dest="plot", help="plot", default=False)
    parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

    try:
        options,args = parser.parse_args(sys.argv[1:])
    except SystemExit:
        # let optparse handle --help and usage errors with its own exit code
        raise
    except:
        print("Error: check usage with  -h ghost_ccf_pipeline.py")
        sys.exit(1)

    if options.verbose:
        print('Spectral s1d.fits data pattern: ', options.input)
        if options.ccf_mask != "":
            print('Input CCF mask: ', options.ccf_mask)
        if options.output_template != "":
            print('Output template spectrum: ', options.output_template)
        if options.source_rv != 0 :
            print('Input source RV (km/s): ', options.source_rv)
        print('Initial CCF width (km/s): ', options.ccf_width)
        print('Velocity sampling (km/s): ', options.vel_sampling)

    # make list of tfits data files
    if options.verbose:
        print("Creating list of s1d fits spectrum files...")
    inputdata = sorted(glob.glob(options.input))

    max_gap_size = 5.0
    min_window_size = 150.
    vel_sampling = 1.0
    verbose = options.verbose

    number_of_knots_for_normalization = options.nknots
    speed_of_light_in_kps = constants.c / 1000.

    # First load spectra into a container
    array_of_spectra = reduc_lib.load_array_of_ghost_spectra(inputdata, rvfile=options.rvfile, obslog=options.output_obslog_file, apply_berv=True, plot=options.plot, verbose=verbose)

    # Then load data into vector
    spectra = reduc_lib.get_spectral_data(array_of_spectra, verbose=verbose)

    # Use wide values to avoid too much clipping at this point. This will improve the noise model
    #spectra = reduc_lib.get_gapfree_windows(spectra, max_vel_distance=max_gap_size, min_window_size=min_window_size, fluxkey="fluxes", wavekey="waves_sf", verbose=verbose)

    # Set a common wavelength grid for all input spectra
    spectra = reduc_lib.set_common_wl_grid(spectra, vel_sampling=vel_sampling, verbose=verbose)

    # Interpolate all spectra to a common wavelength grid
    spectra = reduc_lib.resample_and_align_spectra(spectra, verbose=verbose, plot=False)
    #spectra["aligned_waves"]
    #spectra["sf_fluxes"],spectra["sf_fluxerrs"]
    #spectra["rest_fluxes"], spectra["rest_fluxerrs"]

    spectra, template = reduc_lib.reduce_spectra(spectra, nsig_clip=5.0, combine_by_median=False, subtract=True, fluxkey="sf_fluxes", fluxerrkey="sf_fluxerrs", wavekey="common_wl", update_spectra=True, plot=False, verbose=True)

    #plt.plot(template["wl"], template["flux"],'-')
    #plt.show()
    #exit()

    #spectra, template = reduc_lib.normalize_spectra(spectra, template, fluxkey="rest_fluxes", fluxerrkey="rest_fluxerrs", cont_function='spline3', polyn_order=40, med_filt=1, plot=True)
    spectra, template = reduc_lib.normalize_spectra(spectra, template, fluxkey="sf_fluxes", fluxerrkey="sf_fluxerrs", cont_function='spline3', polyn_order=number_of_knots_for_normalization, plot=options.plot)

    #plt.plot(template["wl"], template["flux"],'-')
    #plt.show()

    # Calculate statistical weights based on the time series dispersion 1/sig^2
    spectra = reduc_lib.calculate_weights(spectra, template, use_err_model=False, plot=True)

    # Initialize drift containers with zeros
    drifts = reduc_lib.get_zero_drift_containers(inputdata)

    ccf_width = options.ccf_width

    source_rv=options.source_rv
    output_template = options.output_template
    ccf_mask = options.ccf_mask

    # Start dealing with CCF related parameters and construction of a weighted mask
    # load science CCF parameters
    ccf_params = ccf_lib.set_ccf_params(ccf_mask)

    # update ccf width with input value
    ccf_params["CCF_WIDTH"] = float(ccf_width)

    ccfmask = ccf_lib.apply_weights_to_ccf_mask(ccf_params, template["wl"], template["flux"], template["fluxerr"], spectra["weights"], median=True, remove_lines_with_nans=True, source_rv=source_rv, verbose=False, plot=True)

    base_header = deepcopy(array_of_spectra["spectra"][0]["header"])

    template_ccf = ccf_lib.run_ccf_eder(ccf_params, template["wl"], template["flux"], base_header, ccfmask, targetrv=source_rv, normalize_ccfs=True, plot=True, verbose=False)

    source_rv = template_ccf["header"]['RV_OBJ']
    ccf_params["SOURCE_RV"] = source_rv
    ccf_params["CCF_WIDTH"] = 8 * template_ccf["header"]['CCFMFWHM']

    if verbose :
        print("Source RV={:.4f} km/s  CCF FWHM={:.2f} km/s CCF window size={:.2f} km/s".format(source_rv,template_ccf["header"]['CCFMFWHM'],ccf_params["CCF_WIDTH"]))
    # Apply weights to stellar CCF mask
    ccfmask = ccf_lib.apply_weights_to_ccf_mask(ccf_params, template["wl"], template["flux"], template["fluxerr"], spectra["weights"], median=True, remove_lines_with_nans=True, source_rv=source_rv, verbose=verbose, plot=False)

    # plot template
    # blue channel:
    #reduc_lib.plot_template_products_with_CCF_mask(template, ccfmask, source_rv=ccf_params["SOURCE_RV"],  wl0=487, wlf=491, pfilename="")
    # red channel :
    #reduc_lib.plot_template_products_with_CCF_mask(template, ccfmask, source_rv=ccf_params["SOURCE_RV"],  wl0=613, wlf=623, pfilename="")

    if output_template != "" :
        if verbose :
            print("Saving template spectrum to file: {0} ".format(output_template))
        ghostlib.write_spectrum_to_fits(template["wl"], template["flux"], template["fluxerr"], output_template, header=template_ccf["header"])

    ###### START CCF  #######
    save_output = True
    normalize_ccfs = True
    run_analysis = True
    #fluxkey="sf_fluxes"
    fluxkey="rest_fluxes"
    waveskey="aligned_waves"
    plot = options.plot

    if plot :
        templ_legend = "Template of {}".format(template_ccf["header"]["OBJECT"].replace(" ",""))
        plt.plot(template_ccf['RV_CCF'], template_ccf['MEAN_CCF'], "-", color='green', lw=2, label=templ_legend, zorder=2)

    calib_rv, drift_rv  = [], []
    mean_fwhm = []
    sci_ccf_file_list = []

    rvccf_data, ccf_data = [], []

    for i in range(spectra['nspectra']) :

        basename = os.path.basename(spectra['filenames'][i])

        if verbose :
            print("Running CCF on file {0}/{1} -> {2}".format(i,spectra['nspectra']-1,basename))

        rv_drifts = drifts[i]

        fluxes, waves_sf = spectra[fluxkey][i], spectra[waveskey][i]

        # run main routine to process ccf on science fiber
        header = array_of_spectra["spectra"][i]["header"]

        output_ccf_filename = ""
        if options.saveccfs :
            filedir = os.path.dirname(spectra['filenames'][i])
            maskbasename = os.path.basename(options.ccf_mask).split('.')[0]
            output_ccf_filename = os.path.join(filedir, "CCFTABLE_{}_{}".format(maskbasename,basename))

        # run an adapted version of the ccf codes using reduced spectra as input
        sci_ccf = ccf_lib.run_ccf_eder(ccf_params, waves_sf, fluxes, header, ccfmask, rv_drifts=rv_drifts, filename=spectra['filenames'][i], targetrv=ccf_params["SOURCE_RV"], normalize_ccfs=normalize_ccfs, output=output_ccf_filename, plot=False, verbose=False)

        if options.saveccfs :
            sci_ccf_file_list.append(os.path.abspath(sci_ccf["file_path"]))

        #rvtrue = sci_ccf["header"]['RV_OBJ'] * (1 - header["BERV"]/speed_of_light_in_kps)
        calib_rv.append(sci_ccf["header"]['RV_OBJ'] )
        mean_fwhm.append(sci_ccf["header"]['CCFMFWHM'])
        drift_rv.append(sci_ccf["header"]['RV_DRIFT'])

        if verbose :
            print("Spectrum: {0} DATE={1} Sci_RV={2:.5f} km/s RV_DRIFT={3:.5f} km/s".format(os.path.basename(spectra['filenames'][i]), sci_ccf["header"]["DATE"], sci_ccf["header"]['RV_OBJ'], sci_ccf["header"]["RV_DRIFT"]))

        if plot :
            if i == spectra['nspectra'] - 1 :
                scilegend = "{}".format(sci_ccf["header"]["OBJECT"].replace(" ",""))
            else :
                scilegend = None
            #plt.plot(esci_ccf['RV_CCF'],sci_ccf['MEAN_CCF']-esci_ccf['MEAN_CCF'], "--", label="spectrum")
            plt.plot(sci_ccf['RV_CCF'], sci_ccf['MEAN_CCF'], "-", color='#2ca02c', alpha=0.5, label=scilegend, zorder=1)

        rvccf_data.append(sci_ccf['RV_CCF'])
        ccf_data.append(sci_ccf['MEAN_CCF'])

    mean_fwhm = np.array(mean_fwhm)
    velocity_window = 1.5*np.nanmedian(mean_fwhm)

    calib_rv, median_rv = np.array(calib_rv), np.nanmedian(calib_rv)

    if plot :
        plt.xlabel('Velocity [km/s]')
        plt.ylabel('CCF')
        plt.legend()
        plt.show()

        plt.plot(spectra["bjds"], (calib_rv  - median_rv), 'o', color='#2ca02c', label="Sci RV = {0:.4f} km/s".format(median_rv))
        plt.plot(spectra["bjds"], (mean_fwhm  - np.nanmean(mean_fwhm)), '--', color='#2ca02c', label="Sci FWHM = {0:.4f} km/s".format(np.nanmean(mean_fwhm)))

        drift_rv = np.array(drift_rv)

        mean_drift, sigma_drift = np.nanmedian(drift_rv), np.nanstd(drift_rv)
        plt.plot(spectra["bjds"], drift_rv, '.', color='#ff7f0e', label="Inst. FP drift = {0:.4f}+/-{1:.4f} km/s".format(mean_drift,sigma_drift))

        plt.xlabel(r"BJD")
        plt.ylabel(r"Velocity [km/s]")
        plt.legend()
        plt.show()

    if verbose :
        print("Running CCF analysis: velocity_window = {0:.3f} km/s".format(velocity_window))

    obj = sci_ccf["header"]["OBJECT"].replace(" ","")

    # The solution is to work on a simpler version of ccf2rv routines, but it takes a bit of time for coding and testing. To do soon.
    rvs, rverrs, template_ccf = ccf_lib.ccf_analysis(rvccf_data[0], np.array(ccf_data,dtype=float), calib_rv, nsig_clip=0, velocity_window=velocity_window, plot=plot, verbose=verbose)

    if options.output_ccfs_file != "" :
        ccf_lib.save_ccfs_to_fits(template_ccf, spectra["bjds"], options.output_ccfs_file, header=base_header)

    fwhm, fwhmerr, bis, biserr = ccf_lib.ccf_fwhm_and_biss(template_ccf, plot=plot, verbose=verbose)

    if options.output_rv_file != "" :
        spectrallib.save_time_series(options.output_rv_file, spectra["bjds"]-2400000., rvs, rverrs, xlabel="rjd", ylabel="vrad", yerrlabel="svrad", write_header_rows=True)

    if options.output_bis_file != "" :
        spectrallib.save_time_series(options.output_bis_file, spectra["bjds"]-2400000., bis, biserr, xlabel="rjd", ylabel="biss", yerrlabel="bisserr", write_header_rows=True)

    if options.output_fwhm_file != "" :
        spectrallib.save_time_series(options.output_fwhm_file, spectra["bjds"]-2400000., fwhm, fwhmerr, xlabel="rjd", ylabel="fwhm", yerrlabel="fwhmerr", write_header_rows=True)

    source_rv = np.nanmedian(rvs)
    source_rverr = np.nanstd(rvs-source_rv)

    if plot :
        plt.errorbar(spectra["bjds"], 1000*(rvs-np.nanmedian(rvs)), 1000*rverrs, fmt='o', color='#2ca02c', label="Sci RV = {:.4f} +/- {:.4f} km/s".format(source_rv, source_rverr))
        plt.xlabel(r"BJD")
        plt.ylabel(r"Velocity [m/s]")
        plt.legend()
        plt.show()


if __name__ == "__main__":
    main()
