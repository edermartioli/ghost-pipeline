# -*- coding: iso-8859-1 -*-
"""
    Created on May 27 2023
    
    Description: This routine calculate doppler tomography from ccf data
    
    @author: Eder Martioli <emartioli@lna.br>
    Laboratório Nacional de Astrofísica - LNA/MCTI
    
    python ghost_rm_tomography.py --input=/Volumes/Samsung_T5/Science/WASP-108/WASP-108_ghost_blue_ccfs.fits --inputrv=/Volumes/Samsung_T5/Science/WASP-108/WASP-108_ghost_blue_ccfrv.rdb
    python ghost_rm_tomography.py --input=/Volumes/Samsung_T5/Science/WASP-108/WASP-108_ghost_red_ccfs.fits --inputrv=/Volumes/Samsung_T5/Science/WASP-108/WASP-108_ghost_red_ccfrv.rdb


    python ghost_rm_tomography.py --input=/Users/eder/Science/WASP-108/WASP-108_OLD_CCF/WASP-108_ghost_blue_ccfs.fits --inputrv=/Users/eder/Science/WASP-108/WASP-108_OLD_CCF/WASP-108_ghost_blue_ccfrv.rdb
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
from scipy.optimize import curve_fit
from astropy.io import ascii
from scipy import interpolate
from copy import deepcopy
from astropy.convolution import Gaussian1DKernel, convolve

def plot_2d(x, y, z, model=[], LIM=None, LAB=None, z_lim=None, transit=None, title="", pfilename="", cmap="gist_heat"):
    """
        Use pcolor to display 2D map of sequence of spectra
    
    Inputs:
    - x:        x array of the 2D map (if x is 1D vector, then meshgrid; else: creation of Y)
    - y:        y 1D vector of the map
    - z:        2D array (sequence of spectra; shape: (len(x),len(y)))
    - LIM:      list containing: [[lim_inf(x),lim_sup(x)],[lim_inf(y),lim_sup(y)],[lim_inf(z),lim_sup(z)]]
    - LAB:      list containing: [label(x),label(y),label(z)] - label(z) -> colorbar
    - title:    title of the map
    - **kwargs: **kwargs of the matplolib function pcolor
    
    Outputs:
    - Display 2D map of the sequence of spectra
    
    """
    
    if len(np.shape(x))==1:
        X,Y  = np.meshgrid(x,y)
    else:
        X = x
        Y = []
        for n in range(len(x)):
            Y.append(y[n] * np.ones(len(x[n])))
        Y = np.array(Y,dtype=float)
    Z = z

    if LIM == None :
        x_lim = [np.min(X),np.max(X)] #Limits of x axis
        y_lim = [np.min(Y),np.max(Y)] #Limits of y axis
        if z_lim == None :
            z_lim = [np.min(Z),np.max(Z)]
        LIM   = [x_lim,y_lim,z_lim]

    if LAB == None :
        ### Labels of the map
        x_lab = r"$Wavelength$ [nm]"   #Wavelength axis
        y_lab = r"Time [BJD]"         #Time axis
        z_lab = r"Flux"     #Intensity (exposures)
        LAB   = [x_lab,y_lab,z_lab]

    fig = plt.figure()
    plt.rcParams["figure.figsize"] = (10,7)
    ax = plt.subplot(111)

    if transit :
        ax.hlines(y=transit["tini"],xmin=LIM[0][0],xmax=LIM[0][1],ls=':',color='k', lw=2)
        ax.hlines(y=transit["tcen"],xmin=LIM[0][0],xmax=LIM[0][1],ls='--',color='k', lw=2)
        ax.hlines(y=transit["tend"],xmin=LIM[0][0],xmax=LIM[0][1],ls=':',color='k', lw=2)

    if len(model) :
        #print("Input model:", model)
        ax.plot(model, Y, ls='--',color='k', lw=2)

    cc = ax.pcolor(X, Y, Z, vmin=LIM[2][0], vmax=LIM[2][1], cmap=cmap)
    cb = plt.colorbar(cc,ax=ax)
    
    ax.set_xlim(LIM[0][0],LIM[0][1])
    ax.set_ylim(LIM[1][0],LIM[1][1])
    
    ax.set_xlabel(LAB[0])
    ax.set_ylabel(LAB[1],labelpad=15)
    cb.set_label(LAB[2],rotation=270,labelpad=30)

    ax.set_title(title,pad=35)

    if pfilename=="" :
        plt.show()
    else :
        plt.savefig(pfilename, format='png')
    plt.clf()
    plt.close()




def plot_ccfs(rv,tccf,ccfs,residuals) :

    nspc = len(ccfs)

    fig,ax = plt.subplots(nrows = 2, ncols = 1, sharex=True)
    
    for i in range(nspc):
        color = [i/nspc,1-i/nspc,1-i/nspc]
        ax[0].plot(rv, tccf, color = "green", lw=2, label="median CCF")
        ax[0].plot(rv, ccfs[i], color = color, alpha = 0.2)
        ax[1].plot(rv, residuals[i], color = color,alpha = 0.2)

    ax[0].set_title('Mean CCFs', fontsize=20)
    ax[0].set_xlabel('Velocity [km/s]', fontsize=20)
    ax[0].set_ylabel('CCF depth', fontsize=20)
    ax[0].tick_params(axis='both', which='major', labelsize=16)
    ax[0].tick_params(axis='both', which='minor', labelsize=12)

    ax[1].set_title('Residual CCFs', fontsize=20)
    ax[1].set_xlabel('Velocity [km/s]', fontsize=20)
    ax[1].set_ylabel('CCF residual depth', fontsize=20)
    ax[1].tick_params(axis='both', which='major', labelsize=16)
    ax[1].tick_params(axis='both', which='minor', labelsize=12)

    plt.tight_layout()
    #plt.legend()
#    plt.xticks(fontsize=16)
#    plt.yticks(fontsize=16)
    
    plt.show()

def func(x, a, b) :
    return a + b*x

def gauss(v,v0,ew,zp,amp):
    # gaussian with a constant offset. As we know that the ccfs are negative structures, amp will be negative
    return zp+amp*np.exp( -0.5*(v-v0)**2/ew**2)

parser = OptionParser()
parser.add_option("-1", "--input", dest="input", help="Input spectral CCF FITS data",type='string',default="")
parser.add_option("-r", "--inputrv", dest="inputrv", help="Input RV data rdb file",type='string',default="")
parser.add_option("-p", action="store_true", dest="plot", help="plot", default=False)
parser.add_option("-c", action="store_true", dest="convolve", help="convolve ccfs", default=False)
parser.add_option("-v", action="store_true", dest="verbose", help="verbose", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with -h ghost_rm_tomography.py")
    sys.exit(1)

print('Input data file: ', options.input)
print('Input RV data file: ', options.inputrv)

ghost_dir = os.path.dirname(__file__) + '/'

transit = {}
# WASP-108
#transit["tini"] = 2460076.613 - 180.31/(60*24*2)
#transit["tcen"] = 2460076.613
#transit["tend"] = 2460076.613 + 180.31/(60*24*2)

# WASP-101
transit["tini"] = 2460682.6336
transit["tcen"] = 2460682.680
transit["tend"] = 2460682.7247


# WASP-87
transit["tini"] = 2460701.7040
transit["tcen"] = 2460701.7040 + 180. / (60*24*2)
transit["tend"] = 2460701.7040 + 180. / (60*24)

tbl = ascii.read(options.inputrv, data_start=2)
bjd, rv, erv = np.array(tbl['rjd'])+2400000., np.array(tbl['vrad']), np.array(tbl['svrad'])

keep_for_fit = (bjd <= transit["tini"]) | (bjd >= transit["tend"])
popt, pcov = curve_fit(func, bjd[keep_for_fit], rv[keep_for_fit])

rv_model = func(bjd, *popt)

rvfulloffset = np.abs(rv_model[-1] - rv_model[0])
print("RVs full offset: {:.4f} km/s".format(rvfulloffset))

if options.plot :
    plt.plot(bjd-2400000, rv, 'ko', label='RV data')
    plt.plot(bjd[keep_for_fit]-2400000, rv[keep_for_fit], '.', ms=20, alpha=0.4, label='out-of-transit data')
    plt.plot(bjd-2400000, rv_model, 'r--', label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
    plt.xlabel(r"BJD-2400000",fontsize=18)
    plt.ylabel(r"RV [km/s]",fontsize=18)
    plt.legend(fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.show()

rmeffect = rv - rv_model

hdul = fits.open(options.input)
"""
print(hdul.info())
  1  RV            1 ImageHDU         7   (824,)   float64
  2  MEANCCF       1 ImageHDU         7   (824,)   float64
  3  CCFS          1 ImageHDU         8   (824, 28)   float64
  4  RESIDUAL_CCFS    1 ImageHDU         8   (824, 28)   float64
  5  TIMES         1 ImageHDU         7   (28,)   float64
"""

rvs = hdul["RV"].data
tccf = hdul["MEANCCF"].data
ccfs = hdul["CCFS"].data
residuals = hdul["RESIDUAL_CCFS"].data
bjds = hdul["TIMES"].data
nspc = len(ccfs)

if options.plot :
    plot_ccfs(rvs,tccf,ccfs,residuals)

# subtract RVs and stack CCF data in the tpl_rvs grid
outccfs = []
for i in range(nspc):
    if (bjds[i] <= transit["tini"]) or (bjds[i] >= transit["tend"]) :
        outccfs.append(ccfs[i])
        #color = [i/nspc,1-i/nspc,1-i/nspc]
        #plt.plot(tpl_rvs,f(tpl_rvs),lw=1,alpha=0.5,color=color)
mccf = np.nanmean(outccfs,0)

if options.convolve :
    tp0 = [rvs[np.argmin(mccf)],1,1,-0.1]
    tfit, tpcov = curve_fit(gauss, rvs, mccf, p0 = tp0)
    tfit_err = np.sqrt(np.diag(tpcov))
    
    fwhm = 2*np.sqrt(2*np.log(2)) * tfit[1]
    fwhmerr = 2*np.sqrt(2*np.log(2)) * tfit_err[1]
    rvstep = np.median(np.abs(rvs[1:]-rvs[:-1]))

    # degrade resolution in X%
    degradation = 1.1
    sigma_kernel = (tfit[1]/2) * degradation / rvstep

    #gausskernel = Gaussian1DKernel(stddev=kernel_size)
    gauss_1D_kernel = Gaussian1DKernel(sigma_kernel)
    mccf = convolve(mccf, gauss_1D_kernel, boundary='fill', fill_value=1.0)
    #plt.plot(rvs, mccf, 'o', color="green", zorder=1)
    #plt.plot(rvs,conv_mccf,'o',color="blue",zorder=1)
    #plt.plot(rvs,gauss(rvs,*tfit),'-',color="red",zorder=2)
    #plt.show()

newresiduals = []

for i in range(nspc) :
    rv_corr = rvs - rmeffect[i]
    f = interpolate.interp1d(rv_corr, deepcopy(ccfs[i]), kind='cubic', fill_value='extrapolate')
    if options.convolve :
        conv_ccf = convolve(f(rvs), gauss_1D_kernel, boundary='fill', fill_value=1.0)
        resccf = conv_ccf / mccf
    else :
        resccf = f(rvs) / mccf
    newresiduals.append(resccf)
    
newresiduals = np.array(newresiduals,dtype=float)

sig = []
for i in range(len(newresiduals)) :
    sigma = np.std(newresiduals)
    sig.append(sigma)
mediansig = np.nanmedian(sig)
medianrv = np.nanmedian(newresiduals.flatten())

for i in range(len(newresiduals)) :
    color = [i/nspc,1-i/nspc,1-i/nspc]
    plt.plot(rvs,newresiduals[i]+6*mediansig*i,lw=2,color=color)
plt.show()

x_lab = r"Barycentric RV [m/s]"   #Wavelength axis
y_lab = r"Time [BJD]"         #Time axis
z_lab = r"CCF"     #Intensity (exposures)
LAB   = [x_lab,y_lab,z_lab]

plot_2d(rvs, bjds, newresiduals, model=[], LIM=None, LAB=LAB,  z_lim = [medianrv-5*mediansig,medianrv+5*mediansig], transit=transit, title="", pfilename="", cmap="gist_heat")
