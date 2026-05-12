# -*- coding: iso-8859-1 -*-
"""
    Created on April 26 2023
    
    Description: This routine runs several DRAGONS routines to reduce the GHOST spectra
    
    @author: Eder Martioli <emartioli@lna.br>
    Laboratório Nacional de Astrofísica - LNA/MCTI
    
    Simple usage examples:

    First install DRAGONS-3.0 and GHOST routines following installation instructions on:
        https://ghost-drtutorial.readthedocs.io/en/release-3.0.x/sv_installation.html

    In ~/.geminidr/, create or edit the configuration file rsys.cfg as follow:
    [calibs]
    standalone = True
    database_dir = <path_to_my_data>

    Then save the file ghost_pipeline.py in a given $PATH_TO_GHOST_PIPELINE_TOOL

    conda activate ghost-sv
    cd $PATH_TO_GHOST_PIPELINE_TOOL

    ATTENTION! the option "-r" runs the command "caldb init --wipe", which can erase previous information in the calibration database.
    
    Two-step example to reduce first a standard "CD -32 9927" and then a scientific target "XX Oph":
    ------------------------------------------------------------------------------------------------
    python -W"ignore" ghost_pipeline.py --input=/Volumes/Samsung_T5/Gemini/GHOST-SV/DRS/ghost_tutorial/playdata/example1/*.fits \
                                        --reduced_data_dir=/Volumes/Samsung_T5/Gemini/GHOST-SV/DRS/ghost_tutorial/playground \
                                        --object="CD -32 9927" -s \
                                        -bfao1r
    
    python -W"ignore" ghost_pipeline.py --input=/Volumes/Samsung_T5/Gemini/GHOST-SV/DRS/ghost_tutorial/playdata/example1/*.fits \
                                        --reduced_data_dir=/Volumes/Samsung_T5/Gemini/GHOST-SV/DRS/ghost_tutorial/playground \
                                        --object="XX Oph" \
                                        --red_std_spec=S20230416S0073_red001_standard.fits \
                                        --blue_std_spec=S20230416S0073_blue001_standard.fits \
                                        -o1

    Example for a full reduction of a science target without standards :
    --------------------------------------------------------------------
    python  -W"ignore" ghost_pipeline.py --input=/Volumes/Samsung_T5/Gemini/GHOST-SV/DRS/ghost_tutorial/playdata/example1/*.fits \
                                         --reduced_data_dir=/Volumes/Samsung_T5/Gemini/GHOST-SV/DRS/ghost_tutorial/playground \
                                         --object="XX Oph" \
                                         -bfao1r
                                         
    Example for a full reduction of all science targets without standards :
    ----------------------------------------------------------------------
    python  -W"ignore" ghost_pipeline.py --input=/Volumes/Samsung_T5/Gemini/GHOST-SV/DRS/ghost_tutorial/playdata/example1/*.fits \
                                         --reduced_data_dir=/Volumes/Samsung_T5/Gemini/GHOST-SV/DRS/ghost_tutorial/playground \
                                         -bfao1r
                                         

    Example for the reduction of SV data :
    ----------------------------------------------------------------------
    # To reduce WASP-108 transit data
    python -W"ignore" /Volumes/Samsung_T5/Gemini/ghost-pipeline/ghost_pipeline.py --input=/Volumes/Samsung_T5/Data/GHOST/wasp108/*.fits --reduced_data_dir=/Volumes/Samsung_T5/Data/GHOST/reduced_WASP-108/ --object="Sky" --detector_x_bin=2 --detector_y_bin=4 -bfaort

    # To reduce 18 Sco data
    python -W"ignore" /Volumes/Samsung_T5/Gemini/ghost-pipeline/ghost_pipeline.py --input=/Volumes/Samsung_T5/Data/GHOST/18Sco/*.fits --reduced_data_dir=/Volumes/Samsung_T5/Data/GHOST/reduced_18Sco/ --object="18 Sco" --detector_x_bin=1 --detector_y_bin=1 -bfaor
    
    # To reduce WASP-189 HR data
    python -W"ignore" /Volumes/Samsung_T5/Gemini/ghost-pipeline/ghost_pipeline.py --input=/Volumes/Samsung_T5/Data/GHOST/wasp189/*.fits --reduced_data_dir=/Volumes/Samsung_T5/Data/GHOST/reducedWASP189/ --object="WASP-189" --detector_x_bin=1 --detector_y_bin=4 -bfaor
    """

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os,sys
import glob
import astropy.io.fits as fits

def match_object(inputdata, object_name="") :
    outputdata = []
    for i in range(len(inputdata)) :
        hdr = fits.getheader(inputdata[i],0)
        #print(i, inputdata[i], hdr["OBSTYPE"], hdr["OBJECT"], object_name)
        if hdr["OBSTYPE"] == 'OBJECT':
            if hdr["OBJECT"] == object_name or object_name == "":
                print(i,hdr["OBJECT"],inputdata[i])
                outputdata.append(inputdata[i])
    return outputdata

parser = OptionParser()
parser.add_option("-i", "--input", dest="input", help="Input spectral FITS data pattern",type='string',default="*.fits")
parser.add_option("-O", "--object", dest="object", help="Object name",type='string',default="")
parser.add_option("-d", "--reduced_data_dir", dest="reduced_data_dir", help="Reduced data dir",type='string',default="./")
parser.add_option("-R", "--red_std_spec", dest="red_std_spec", help="Red standard FITS spectrum",type='string',default="")

parser.add_option("-x", "--detector_x_bin", dest="detector_x_bin", help="CCD spectral x-binning",type='int',default=1)
parser.add_option("-y", "--detector_y_bin", dest="detector_y_bin", help="CCD spatial y-binning",type='int',default=1)

parser.add_option("-B", "--blue_std_spec", dest="blue_std_spec", help="Blue standard FITS spectrum",type='string',default="")
parser.add_option("-s", action="store_true", dest="standard", help="Object is standard", default=False)
parser.add_option("-v", action="store_true", dest="display", help="display", default=False)
parser.add_option("-b", action="store_true", dest="run_bias", help="run bias reduction", default=False)
parser.add_option("-f", action="store_true", dest="run_flat", help="Run flat reduction", default=False)
parser.add_option("-a", action="store_true", dest="run_arc", help="Run arc reduction", default=False)
parser.add_option("-o", action="store_true", dest="run_object", help="Run object reduction", default=False)
parser.add_option("-n", action="store_true", dest="no_sky_subtraction", help="No sky subtraction", default=False)
parser.add_option("-t", action="store_true", dest="object_in_ifu2", help="Object in IFU2", default=False)
parser.add_option("-1", action="store_true", dest="run_onedspec", help="Generate 1D spectra", default=False)
parser.add_option("-r", action="store_true", dest="reset_caldb", help="reset caldb", default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with -h ghost_pipeline.py")
    sys.exit(1)

print('GHOST data pattern: ', options.input)
print('Input object name: ', options.object)
print('Reduced data dir: ', options.reduced_data_dir)

ghost_dir = os.path.dirname(__file__) + '/'
reduced_data_dir = options.reduced_data_dir

object_name = options.object
inputdata = sorted(glob.glob(options.input))

# get list of object data files
objectdata = match_object(inputdata, object_name=object_name)

# set input control switches
display_data_contents = options.display
run_bias = options.run_bias
run_flat = options.run_flat
run_arc = options.run_arc
run_object = options.run_object
run_onedspec = options.run_onedspec
obj_is_standard = options.standard

delete_intermediate_products = False # do not turn this off.
# If delete_intermediate_products=False, some intermediate products won't be deleted and they
# can mix up with others and mess the reduction. It probably needs better wild cards than just *.fits to select files

#############################
##### Access reduced data dir
#############################
if not os.path.exists(reduced_data_dir) :
    os.mkdir(reduced_data_dir)

os.chdir(reduced_data_dir)

#############################
##### Display data contents ####
#############################
if options.reset_caldb :
    command = "caldb init --wipe"
    print("Running: ",command)
    os.system(command)
    
#############################
##### Display data contents ####
#############################
if display_data_contents :
    command = "showd --adpkg=ghost_instruments {} -d object,detector_x_bin,detector_y_bin,read_mode".format(options.input)
    print("Running: ",command)
    os.system(command)

#############################
##### RUN BIAS reduction ####
#############################
if run_bias :

    command = "dataselect --adpkg=ghost_instruments {} --tags BIAS -o biasbundles.lis".format(options.input)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @biasbundles.lis"
    print("Running: ",command)
    os.system(command)

    command = "dataselect --adpkg=ghost_instruments {}/*_slit.fits --tags BIAS,SLIT -o biasslit.lis".format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @biasslit.lis"
    print("Running: ",command)
    os.system(command)

    command = 'dataselect --adpkg=ghost_instruments {}/*_red*.fits --tags BIAS,RED --expr="detector_x_bin=={} and detector_y_bin=={}" -o biasredsci.lis'.format(reduced_data_dir, options.detector_x_bin, options.detector_y_bin)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @biasredsci.lis"
    print("Running: ",command)
    os.system(command)

    command = 'dataselect --adpkg=ghost_instruments {}/*_blue*.fits --tags BIAS,BLUE --expr="detector_x_bin=={} and detector_y_bin=={}" -o biasbluesci.lis'.format(reduced_data_dir, options.detector_x_bin, options.detector_y_bin)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @biasbluesci.lis"
    print("Running: ",command)
    os.system(command)

    command = 'dataselect --adpkg=ghost_instruments {}/*_red*.fits --tags BIAS,RED --expr="detector_x_bin==1 and detector_y_bin==1" -o biasredflatarc.lis'.format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @biasredflatarc.lis"
    print("Running: ",command)
    os.system(command)

    command = 'dataselect --adpkg=ghost_instruments {}/*_blue*.fits --tags BIAS,BLUE --expr="detector_x_bin==1 and detector_y_bin==1" -o biasblueflatarc.lis'.format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @biasblueflatarc.lis"
    print("Running: ",command)
    os.system(command)

    command = "caldb add {}/calibrations/processed_bias/*.fits".format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    if delete_intermediate_products :
        command = "rm {}/*.fits".format(reduced_data_dir)
        print("Running: ",command)
        os.system(command)


#############################
##### RUN FLAT reduction ####
#############################
if run_flat :

    command = "dataselect --adpkg=ghost_instruments {} --tags FLAT -o flatbundles.lis".format(options.input)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @flatbundles.lis"
    print("Running: ",command)
    os.system(command)

    command = "dataselect --adpkg=ghost_instruments {}/*_slit.fits --tags SLITFLAT -o slitflat.lis".format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @slitflat.lis"
    print("Running: ",command)
    os.system(command)

    command = "caldb add {}/calibrations/processed_slitflat/*.fits".format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    command = 'dataselect --adpkg=ghost_instruments {}/*_red*.fits --tags FLAT,RED -o flatred.lis'.format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @flatred.lis"
    print("Running: ",command)
    os.system(command)

    command = 'dataselect --adpkg=ghost_instruments {}/*_blue*.fits --tags FLAT,BLUE -o flatblue.lis'.format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @flatblue.lis"
    print("Running: ",command)
    os.system(command)

    command = "caldb add {}/calibrations/processed_flat/*.fits".format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    if delete_intermediate_products :
        command = "rm {}/*.fits".format(reduced_data_dir)
        print("Running: ",command)
        os.system(command)

#############################
##### RUN ARC reduction ####
#############################
if run_arc :

    command = "dataselect --adpkg=ghost_instruments {} --tags ARC -o arcbundles.lis".format(options.input)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @arcbundles.lis"
    print("Running: ",command)
    os.system(command)

    command = "dataselect --adpkg=ghost_instruments {}/*_slit.fits --tags ARC,SLIT | head -n 1 > arcslit.lis".format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @arcslit.lis"
    print("Running: ",command)
    os.system(command)

    command = "caldb add {}/calibrations/processed_slit/*.fits".format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    command = 'dataselect --adpkg=ghost_instruments {}/*_red*.fits --tags ARC,RED -o arcred.lis'.format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @arcred.lis"
    print("Running: ",command)
    os.system(command)

    command = 'dataselect --adpkg=ghost_instruments {}/*_blue*.fits --tags ARC,BLUE -o arcblue.lis'.format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @arcblue.lis"
    print("Running: ",command)
    os.system(command)

    command = "caldb add {}/calibrations/processed_arc/*.fits".format(reduced_data_dir)
    print("Running: ",command)
    os.system(command)

    if delete_intermediate_products :
        command = "rm {}/*.fits".format(reduced_data_dir)
        print("Running: ",command)
        os.system(command)

#############################
##### RUN OBJECT reduction ##
#############################
if run_object :

    for i in range(len(objectdata)) :

        basename_with_fitsext = os.path.basename(objectdata[i])
        basename = basename_with_fitsext.split(".")[0]
        
        flags = ""
        if obj_is_standard :
            flags = " -r reduceStandard"
        # the flag below can be used if an object falls at IFU2 instead of sky.
        if options.object_in_ifu2 :
            flags += " -p extractProfile:ifu2=object"

        stdredflag, stdblueflag = "", ""
        if options.red_std_spec != "" :
            stdredflag = " -p standard={}".format(options.red_std_spec)
        if options.blue_std_spec != "" :
            stdblueflag = " -p standard={}".format(options.blue_std_spec)

        command = 'dataselect --adpkg=ghost_instruments {} --expr="{}" -o objbundles.lis'.format(objectdata[i], "object=='{}'".format(object_name))
        print("Running: ",command)
        os.system(command)

        command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @objbundles.lis"
        print("Running: ",command)
        os.system(command)
    
        command = "dataselect --adpkg=ghost_instruments {}/{}_slit.fits --tags SLIT -o objslit.lis".format(reduced_data_dir,basename)
        print("Running: ",command)
        os.system(command)

        command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr @objslit.lis"
        print("Running: ",command)
        os.system(command)
    
        command = "caldb add {}/calibrations/processed_slit/{}_slit_*.fits".format(reduced_data_dir,basename)
        print("Running: ",command)
        os.system(command)

        command = 'dataselect --adpkg=ghost_instruments {}/{}_red*.fits --tags RED -o objred.lis'.format(reduced_data_dir,basename)
        print("Running: ",command)
        os.system(command)
        
        if options.no_sky_subtraction :
            command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr{} @objred.lis -p extractProfile:sky_subtract=False -p barycentricCorrect:correction_factor=1{}".format(flags, stdredflag)
        else :
            command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr{} @objred.lis -p extractProfile:sky_subtract=True -p barycentricCorrect:correction_factor=1{}".format(flags, stdredflag)

        print("Running: ",command)
        os.system(command)
    
        command = 'dataselect --adpkg=ghost_instruments {}/{}_blue*.fits --tags BLUE -o objblue.lis'.format(reduced_data_dir,basename)
        print("Running: ",command)
        os.system(command)

        if options.no_sky_subtraction :
            command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr{} @objblue.lis -p extractProfile:sky_subtract=False -p barycentricCorrect:correction_factor=1{}".format(flags, stdblueflag)
        else :
            command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr{} @objblue.lis -p extractProfile:sky_subtract=True -p barycentricCorrect:correction_factor=1{}".format(flags, stdblueflag)
        print("Running: ",command)
        os.system(command)
    
        if delete_intermediate_products :

            basename_with_fitsext = os.path.basename(objectdata[i])
            basename = basename_with_fitsext.split(".")[0]

            command = "rm {}/{}*_slit.fits".format(reduced_data_dir,basename)
            print("Running: ",command)
            os.system(command)
        
            command = "rm {}/{}*_red???.fits".format(reduced_data_dir,basename)
            print("Running: ",command)
            os.system(command)
        
            command = "rm {}/{}*_blue???.fits".format(reduced_data_dir,basename)
            print("Running: ",command)
            os.system(command)

#############################
##### RUN 1D SPECTRUM #######
#############################
if run_onedspec :
    for i in range(len(objectdata)) :
        basename_with_fitsext = os.path.basename(objectdata[i])
        basename = basename_with_fitsext.split(".")[0]

        inputredproducts = sorted(glob.glob("{}/{}_red*_dragons.fits".format(reduced_data_dir,basename)))

        for j in range(len(inputredproducts)) :
            command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr -r write1DSpectra {}".format(inputredproducts[j])
            print("Running: ",command)
            os.system(command)

        inputblueproducts = sorted(glob.glob("{}/{}_blue*_dragons.fits".format(reduced_data_dir,basename)))

        for j in range(len(inputblueproducts)) :
            command = "reduce --adpkg=ghost_instruments --drpkg=ghostdr -r write1DSpectra {}".format(inputblueproducts[j])
            print("Running: ",command)
            os.system(command)
