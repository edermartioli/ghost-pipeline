#!/bin/bash

## Part 1 - BIAS

caldb init -w 

dataselect ../*.fits --tags BIAS -o biasbundles.lis

reduce @biasbundles.lis

dataselect *.fits --tags BIAS,SLIT -o biasslit.lis

reduce @biasslit.lis

#OS BIAS ABAIXO DEVEM TER BINNING CORRESPONDENTE AO BINNING DOS FRAMES DE CIÊNCIA
dataselect *.fits --tags BIAS,RED --expr="binning=='1x2'" -o biasredsci.lis
dataselect *.fits --tags BIAS,BLUE --expr="binning=='1x2'" -o biasbluesci.lis

reduce @biasredsci.lis
reduce @biasbluesci.lis

#OS BIAS PARA REDUZIR OS FLATS DEVEM SEMPRE SER 1X1, PORQUE OS FLATS SÃO SEMPRE 1X1
dataselect *.fits --tags BIAS,RED --expr="binning=='1x1'" -o biasredflatarc.lis
dataselect *.fits --tags BIAS,BLUE --expr="binning=='1x1'" -o biasblueflatarc.lis

reduce @biasredflatarc.lis
reduce @biasblueflatarc.lis

caldb add calibrations/processed_bias/*.fits

rm *.fits

## Part 2 - FLAT -> OS FLATS SÃO BINADOS PARA OS BINS DA CIÊNCIA DEPOIS

dataselect ../*.fits --tags FLAT -o flatbundles.lis

reduce @flatbundles.lis

dataselect *.fits --tags SLITFLAT -o slitflat.lis

reduce @slitflat.lis
caldb add calibrations/processed_slitflat/*.fits

dataselect *.fits --tags FLAT,RED -o flatred.lis
dataselect *.fits --tags FLAT,BLUE -o flatblue.lis

# add -p smoothing=6 for out of focus data
reduce @flatred.lis
reduce @flatblue.lis

caldb add calibrations/processed_flat/*.fits

rm *.fits

# Part 3 - ARC

dataselect ../*.fits --tags ARC -o arcbundles.lis
reduce @arcbundles.lis

dataselect *.fits --tags ARC,SLIT | head -n 1 > arcslit.lis
reduce @arcslit.lis

caldb add calibrations/processed_slit/*.fits

dataselect *.fits --tags ARC,RED -o arcred.lis
dataselect *.fits --tags ARC,BLUE -o arcblue.lis

reduce @arcred.lis
reduce @arcblue.lis
caldb add calibrations/processed_arc/*.fits

rm *.fits

## Part 4 - Sci bundle

dataselect ../*.fits  --xtags CAL -o scibundles.lis
reduce @scibundles.lis

dataselect *.fits --tags SLIT -o scislit.lis
reduce @scislit.lis

caldb add calibrations/processed_slit/*fits

dataselect *.fits --tags  RED  --xtags CAL -o scired.lis
dataselect *.fits --tags BLUE  --xtags CAL -o sciblue.lis
#
### Part 5 - Sci reduce
# TO REDUCE WITHOUT A SPECPHOTOMETRIC STANDARD ADD:
#-p fluxCalibrate:do_cal=skip TO THE REDUCE INSTANCE
reduce -p fluxCalibrate:do_cal=skip extractSpectra:sky_subtract=False combineOrders:stacking_mode=none @scired.lis
reduce -p fluxCalibrate:do_cal=skip extractSpectra:sky_subtract=False combineOrders:stacking_mode=none @sciblue.lis
