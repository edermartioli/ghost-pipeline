# -*- coding: utf-8 -*-
"""
Print an observation log for a set of raw GHOST FITS files.

For each input file, this tool prints one row with the main header
keywords (OBSTYPE, OBJECT, DATE, UT, resolution and read modes), which is
useful to inspect the contents of a GHOST data package before reduction.

Usage example
-------------
::

    ghost_log --input="/path/to/data/*.fits"

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


def main():
    """Entry point: parse command-line options and print the observation log."""
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input", help="Input spectral FITS data pattern",type='string',default="*.fits")
    try:
        options,args = parser.parse_args(sys.argv[1:])
    except SystemExit:
        # let optparse handle --help and usage errors with its own exit code
        raise
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


if __name__ == "__main__":
    main()
