# -*- coding: iso-8859-1 -*-
"""
    Created on May 28 2023
    
    Description: This routine calculate differential RV time series
    
    @author: Eder Martioli <emartioli@lna.br>
    Laboratório Nacional de Astrofísica - LNA/MCTI
    
    python ghost_diff_rv_timeseries.py --input1="/Volumes/Samsung_T5/Science/WASP-108/WASP-108_ghost_red_ccfrv.rdb" --input2="/Volumes/Samsung_T5/Science/WASP-108/TYC8254-1616-1_ghost_red_ccfrv.rdb"
    
    python ghost_diff_rv_timeseries.py --input1="/Volumes/Samsung_T5/Science/WASP-108/WASP-108_ghost_blue_ccfrv.rdb" --input2="/Volumes/Samsung_T5/Science/WASP-108/TYC8254-1616-1_ghost_blue_ccfrv.rdb"
    
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
import glob
import astropy.io.fits as fits
import ghostlib
import matplotlib.pyplot as plt
import numpy as np
import spectrallib

from astropy.io import ascii

parser = OptionParser()
parser.add_option("-1", "--input1", dest="input1", help="Input RV data for object 1",type='string',default="")
parser.add_option("-2", "--input2", dest="input2", help="Input RV data for object 2 (comparison)",type='string',default="")
parser.add_option("-p", action="store_true", dest="plot", help="plot", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
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
