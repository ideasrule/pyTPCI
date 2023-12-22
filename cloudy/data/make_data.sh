#!/bin/bash
#
# This script contains commands to create some of the data files
# present in this directory (data/).  Users are advised to treat
# these commands as templates for recreating the existent files,
# or for creating variations of them.  Either way, consulting
# Hazy 1 for details is strongly encouraged.
#


## this will use default of 10 bins
#../source/cloudy.exe -e "compile agb-silicate ism grains"

## this will do 1 large bin
#../source/cloudy.exe -e "compile agb-silicate ism grains 1"

## this will do 1 large bin for pahs
#../source/cloudy.exe -e "compile pah ab08 grains 1"

## AB08 size distribution, charged PAH, 10 size bins
#../source/cloudy.exe -e 'compile grain 10 bins "ph2c.rfi" ab08'

## AB08 size distribution, neutral PAH, 10 size bins
#../source/cloudy.exe -e 'compile grain 10 bins "ph2n.rfi" ab08'

## this will compile all grain species
../source/cloudy.exe -e "compile grains"


## this will recreate he_iso_recomb.dat
../source/cloudy.exe -e "test; compile recomb coeff he-like"

## this will recreate h_iso_recomb.dat
../source/cloudy.exe -e "test; compile recomb coeff h-like"


## this will recompile the stellar atmospheres
../source/cloudy.exe -e "compile stars"
