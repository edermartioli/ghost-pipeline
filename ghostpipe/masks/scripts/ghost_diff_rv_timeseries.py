# -*- coding: utf-8 -*-
"""
Differential radial-velocity time series between two objects.

Given two RV time series in rdb format (e.g., target and a simultaneous
comparison star observed in the second IFU), this tool subtracts the RV
drift measured on the comparison from the target RVs and saves the
corrected differential time series to a new rdb file (suffix ``_corr``).

Usage example
-------------
::

    ghost_diff_rv_timeseries --input1=target_ccfrv.rdb \
                             --input2=comparison_ccfrv.rdb -p

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
import os,sys
import matplotlib.pyplot as plt
import numpy as np
from ghostpipe import spectrallib

from astropy.io import ascii


def main():
    """Entry point: parse command-line options and compute differential RVs."""
    parser = OptionParser()
    parser.add_option("-1", "--input1", dest="input1", help="Input RV data for object 1",type='string',default="")
    parser.add_option("-2", "--input2", dest="input2", help="Input RV data for object 2 (comparison)",type='string',default="")
    parser.add_option("-p", action="store_true", dest="plot", help="plot", default=False)

    try:
        options,args = parser.parse_args(sys.argv[1:])
    except SystemExit:
        # let optparse handle --help and usage errors with its own exit code
        raise
    except:
        print("Error: check usage with -h ghost_diff_rv_timeseries.py")
        sys.exit(1)

    print('Input Object 1 data pattern: ', options.input1)
    print('Input Object 2 data pattern: ', options.input2)

    ghost_dir = os.path.dirname(__file__) + '/'

    tbl1 = ascii.read(options.input1, data_start=2)
    tbl2 = ascii.read(options.input2, data_start=2)

    rjd1, rv1, erv1 = np.array(tbl1['rjd']), np.array(tbl1['vrad'])*1000, np.array(tbl1['svrad'])*1000
    rjd2, rv2, erv2 = np.array(tbl2['rjd']), np.array(tbl2['vrad'])*1000, np.array(tbl2['svrad'])*1000

    sysrv1 = np.nanmedian(rv1)
    sysrv1err = np.nanstd(rv1)
    sysrv2 = np.nanmedian(rv2)

    rvdrift = rv2 - sysrv2

    rvcorr1 = rv1 - rvdrift

    rvcorr1err = np.sqrt(erv1*erv1 + erv2*erv2)

    plt.errorbar(rjd1,rvcorr1-sysrv1,yerr=rvcorr1err,fmt='.',label="Systemic RV={:.3f}+-{:.3f} km/s".format(sysrv1/1000,sysrv1err/1000))
    plt.xlabel(r"BJD-2400000",fontsize=18)
    plt.ylabel(r"$\Delta$RV [m/s]",fontsize=18)
    plt.legend(fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.show()

    output = (options.input1).replace(".rdb","_corr.rdb")
    spectrallib.save_time_series(output, rjd1, rvcorr1/1000, rvcorr1err/1000, xlabel="rjd", ylabel="vrad", yerrlabel="svrad", write_header_rows=True)


if __name__ == "__main__":
    main()
