#!/bin/bash

# Calculate GDT_TS, GDT_HA and GDT_SC for two molecules with a known amino acid correspondence (targets vs models/predictions)

workdir=$1
mol2process=$2
lga='/ebio/abt1_share/software/LGA_package_src/lga'

cd $workdir

# Create TMP file if it does not exist
mkdir -p TMP

# Align models based on the parameters in the LGA files in the casp tarfiles
$lga -4  -ie  -o1  -sia  -d:4  -gdc_sc  -swap $mol2process > $mol2process.res
# Use LGA alignment records to select residue-residue correspondences for GDT calculations
cat MOL2/$mol2process > MOL2/$mol2process.copy
grep "^LGA " $mol2process.res >> MOL2/$mol2process.copy

# Calculate structural alignment for GDC_TS and GDC_HA
$lga -3 -sia -al $mol2process.copy > $mol2process.gdc_ts_ha
cat  $mol2process.gdc_ts_ha
cat $mol2process.gdc_ts_ha | grep "GDT PERCENT_AT"
# Calculate structural alignment for GDC_SC
$lga -3 -sia -gdc_sc -al $mol2process.copy > $mol2process.gdc_sc
cat $mol2process.gdc_sc | grep "GDT PERCENT_AT"


# calculate GDT_HA and GDT_TS
cat $mol2process.gdc_ts_ha | grep "GDT PERCENT_AT" | awk '{ V=($3+$4+$6+$10)/4.0; printf "GDT_HA = %6.2f\n",V; }'
cat $mol2process.gdc_ts_ha | grep "GDT PERCENT_AT" | awk '{ V=($4+$6+$10+$18)/4.0; printf "GDT_TS = %6.2f\n",V; }'
rm $mol2process.gdc_ts_ha
cat $mol2process.gdc_sc | grep "GDT PERCENT_AT" | awk '{ V=(2*(10.0*$2+9.0*$3+8.0*$4+7.0*$5+6.0*$6+5.0*$7+4.0*$8+3.0*$9+2.0*$10+$11))/(11.0*10.0); printf "GDC_SC = %6.2f\n",V; }'
rm $mol2process.gdc_sc