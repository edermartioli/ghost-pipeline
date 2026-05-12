# -*- coding: iso-8859-1 -*-
"""
    Created on May 26 2023
    
    Description: This routine to reduce the GHOST spectra
    
    @author: Eder Martioli <emartioli@lna.br>
    Laboratório Nacional de Astrofísica - LNA/MCTI
    
    python ghost_log.py --input="/Volumes/Samsung_T5/Data/GHOST/"
    
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
import glob
import astropy.io.fits as fits

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Input spectral FITS data pattern",type='string',default="*.fits")
try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with -h ghost_log.py")
    sys.exit(1)

print('GHOST data pattern: ', options.input)

ghost_dir = os.path.dirname(__file__) + '/'

inputfiles = sorted(glob.glob(options.input))

print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format("IDX","FILENAME","OBSTYPE","OBJECT","DATE","UT","RESOLUT","TARGETM","READRED","READBLU","REDCCDS","BLUCCDS"))
for i in range(len(inputfiles)) :
    hdr = fits.getheader(inputfiles[i],0)

    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(i,os.path.basename(inputfiles[i]),hdr["OBSTYPE"],hdr["OBJECT"],hdr["DATE"],hdr["UT"],hdr["RESOLUT"],hdr["TARGETM"],hdr["READRED"],hdr["READBLU"],hdr["REDCCDS"],hdr["BLUCCDS"]))
