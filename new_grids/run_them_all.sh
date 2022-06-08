#!/bin/bash

# Script to generate all pineappl tables for NNPDF Theory 400

# NOTE: pineappl optimize will fail if the target grid already exists
# that's a feature and not a bug, either remove it yourself or use a different name

exe=../build/Vrap
kinfolder=input_kinematics
output_placeholder=test.pineappl.lz4

${exe} inputE605nlo.dat input_kinematics/E605.dat
pineappl optimize ${output_placeholder} E605nlo.pineappl.lz4

## E886P
${exe} inputE866nlo.dat  input_kinematics/E866P.dat
pineappl optimize ${output_placeholder} E866nlo.pineappl.lz4

## E886R
${exe} inputE866nlo.dat  input_kinematics/E866R.dat
pineappl optimize ${output_placeholder} E866Rnlo.pineappl.lz4

${exe} inputE866deutnlo.dat  input_kinematics/E866R.dat
pineappl optimize ${output_placeholder} E866deutRnlo.pineappl.lz4

## E906
./run_E906.sh inputE906nlo.dat
./run_E906.sh inputE906deutnlo.dat
