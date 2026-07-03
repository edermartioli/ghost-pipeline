# -*- coding: utf-8 -*-
"""
Build a template spectrum out of a time series of GHOST s1d spectra.

The input spectra are resampled onto a common wavelength grid, optionally
shifted to the stellar rest frame using an input RV time series, combined
into a high signal-to-noise template with iterative sigma-clipping, and
optionally normalized to the continuum. The result is saved as a template
FITS product containing the template spectrum and the calibrated time
series of spectra.

Usage examples
--------------
::

    ghost_template_s1d --nknots=40 --input="*red*_s1d.fits" \
        --rv_file=target_ghost_red_ccfrv.rdb \
        --output=target_ghost_red_template.fits -v

    ghost_template_s1d --nknots=10 --input="*blue*_s1d.fits" \
        --rv_file=target_ghost_blue_ccfrv.rdb \
        --output=target_ghost_blue_template.fits -v

Created on May 28, 2023.

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
import glob

import matplotlib.pyplot as plt
from ghostpipe import reduc_lib
from ghostpipe import spectrallib



def main():
    """Entry point: parse command-line options and build the template product."""
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input", help="Spectral *s1d.fits data pattern",type='string',default="*.fits")
    parser.add_option("-o", "--output", dest="output", help="Output time series",type='string',default="")
    parser.add_option("-r", "--rv_file", dest="rv_file", help="Input file with RVs (km/s)",type='string',default="")
    parser.add_option("-a", "--vel_sampling", dest="vel_sampling", help="Velocity sampling for the template spectrum (km/s)",type='float',default=1.8)
    parser.add_option("-n", action="store_true", dest="normalize", help="Normalize spectra to the continuum", default=False)
    parser.add_option("-k", "--nknots", dest="nknots", help="Number of knots for normalization",type='int',default=10)
    parser.add_option("-p", action="store_true", dest="plot", help="plot", default=False)
    parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

    try:
        options,args = parser.parse_args(sys.argv[1:])
    except SystemExit:
        # let optparse handle --help and usage errors with its own exit code
        raise
    except:
        print("Error: check usage with  -h ghost_template_s1d.py")
        sys.exit(1)

    if options.verbose:
        print('Spectral s1d_A.fits data pattern: ', options.input)
        if options.output != "":
            print('Output template product: ', options.output)
        if options.rv_file != "" :
            print('Input file with RVs (km/s) ', options.rv_file)
        print('Velocity sampling (km/s): ', options.vel_sampling)

    # make list of spectral data files
    if options.verbose:
        print("Creating list of s1d fits spectrum files...")
    inputdata = sorted(glob.glob(options.input))

    number_of_knots_for_normalization = options.nknots

    # First load spectra into a container
    array_of_spectra = reduc_lib.load_array_of_ghost_spectra(inputdata, rvfile=options.rv_file, apply_berv=True, plot_spectra=True, plot=True, verbose=options.verbose)

    # Then load data into a vector
    spectra = reduc_lib.get_spectral_data(array_of_spectra, verbose=options.verbose)
    # Set a common wavelength grid for all input spectra
    spectra = reduc_lib.set_common_wl_grid(spectra, vel_sampling=options.vel_sampling, verbose=options.verbose)
    # Interpolate all spectra to a common wavelength grid
    spectra = reduc_lib.resample_and_align_spectra(spectra, verbose=options.verbose, plot=False)
    #The following dict entries contain the spectra: spectra["aligned_waves"], spectra["sf_fluxes"],spectra["sf_fluxerrs"], spectra["rest_fluxes"], spectra["rest_fluxerrs"]
    # Reduce spectra (template matching, sigma clip, etc.) and generate template
    spectra, template = reduc_lib.reduce_spectra(spectra, nsig_clip=5.0, combine_by_median=True, subtract=True, fluxkey="sf_fluxes", fluxerrkey="sf_fluxes", wavekey="common_wl", update_spectra=True, plot=False, verbose=options.verbose)

    if options.normalize :
        # Detect continuum and normalize spectra
        spectra, template = reduc_lib.normalize_spectra(spectra, template, fluxkey="sf_fluxes", fluxerrkey="sf_fluxes", cont_function='spline3', polyn_order=number_of_knots_for_normalization, med_filt=1, plot=False)

    if options.plot :
        plt.plot(template["wl"], template["flux"],'r-', zorder=2)

    halfbinsize = 30
    fluxerrs = []
    for j in range(len(template["flux_arr"])) :
        if options.verbose :
            print("Calculating errors for spectrum {} of {}".format(j+1,len(template["flux_arr"])))

        ferr = spectrallib.spectral_errors_from_residuals(template["wl"], template["flux_residuals"][j], max_delta_wl=0.2, halfbinsize=halfbinsize, min_points_per_bin=3, use_mad=True)
        fluxerrs.append(ferr)

        if options.plot :
            plt.errorbar(template["wl"], template["flux_arr"][j], yerr=ferr, fmt='.', alpha=0.1, zorder=1)


    if options.plot :
        plt.tick_params(axis='x', labelsize=16)
        plt.tick_params(axis='y', labelsize=16)
        plt.minorticks_on()
        plt.tick_params(which='minor', length=3, width=0.7, direction='in', bottom=True, top=True, left=True, right=True)
        plt.tick_params(which='major', length=7, width=1.2, direction='in', bottom=True, top=True, left=True, right=True)
        plt.xlabel("wavelength (nm)",fontsize=18)
        plt.ylabel("flux",fontsize=18)

        plt.show()


    if options.output != "" :
        spectrallib.write_spectra_times_series_to_fits(options.output, template["wl"], template["flux"], template["fluxerr"], spectra["bjds"], template["flux_arr"], fluxerrs, header=spectra["header"][0])


if __name__ == "__main__":
    main()
