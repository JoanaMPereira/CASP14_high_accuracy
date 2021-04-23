#!/bin/bash

# Superimpose two molecules with an unknown amino acid correspondence (targets vs models/predictions)

workdir=$1
mol2process=$2
lga='/ebio/abt1_share/software/LGA_package_src/lga'

cd $workdir

# Create TMP file if it does not exist
mkdir -p TMP
# Align models based on the parameters in the LGA files in the casp tarfiles
$lga -4  -ie  -o1  -sia  -d:4  -gdc_sc  -swap $mol2process > $mol2process.res

cat $mol2process.res
rm $mol2process.res