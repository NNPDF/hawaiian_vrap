#!/bin/bash

# Script to run the very complicated E906R

kinfile=input_kinematics/E906R.dat

generate_fake_kinematics() {
    fake_file=${1}
    mass=${2}
    rapidity=${3}

    if [[ ${mass} == 7.53 ]]
    then
        echo 7.53$RANDOM 0 >> ${fake_file}
    else
        echo ${mass} ${rapidity} >> ${fake_file}
    fi
}

i=0
bin=0
tmp_kinfile=/tmp/e906_kinfile.dat
runcard=${1} # ./inputE906nlo.dat
rm -f ${tmp_kinfile}


if [[ $runcard =~ "deut" ]]
then
    add_deut="deut"
fi

while read -r line
do

    generate_fake_kinematics ${tmp_kinfile} ${line}
    i=$(( i+1 ))
    if [[ i -eq 6 ]]
    then
        echo "Run and reset!"
        ../build/Vrap ${runcard}  ${tmp_kinfile}
        mv test.pineappl.lz4 E906${add_deut}nlo_bin_0${bin}.pineappl.lz4
        rm -f ${tmp_kinfile}
        i=0
        bin=$(( bin+1 ))
    fi

done < ${kinfile}

echo "Probably finished!"
