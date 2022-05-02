#!/bin/bash
#
# Runner for the regression tests for vrap
# note that the results are not exactly compatible with the cfactors from NNPDF for several reasons:
# 1) Simplified kinematics (in order to run faster)
# 2) No extra factors added (NNPDF cfactors had some extra flux normalization added to them)
# 3) Different PDFs

vrap=../build/Vrap
outputfolder=output
verbose=n # y/n

run_for_kinematics() {
    runcard=$1
    kinfile=$2
    outfile=$3

    echo " > Running vrap for $runcard"
    echo "# mass rapidity     xs" > ${outfile}

    while read -r line
    do
        if [[ $verbose == "y" ]]
        then
            res=$( ${vrap} ${runcard} ${line} | tee /dev/tty )
        else
            echo " > > m, y = " ${line}
            res=$( ${vrap} ${runcard} ${line} )
        fi
        val=$(echo $res | awk '{print $NF}')
        echo "     res=" ${val}
        echo ${line}     ${val} >> ${outfile}
    done < $kinfile
}

compare_two() {
    atol=1e-7
    rtol=1e-4
    python -c "
import numpy as np
a = np.loadtxt('${1}')
b = np.loadtxt('${2}')
np.testing.assert_allclose(a,b, atol=${atol}, rtol=${rtol})
    "
}

mkdir -p ${outputfolder}
unset failures

for dataset in E605 E866P E866deut
do
    for order in nlo nnlo
    do
        input_file=${dataset}${order}.dat
        output_file=${outputfolder}/${input_file}
        run_for_kinematics input${input_file} input_kinematics/${dataset}.dat ${output_file}
        compare_two ${output_file} regression_data/${input_file} || failures="${failures} ${input_file}"
    done
done

if [ -v ${failures} ]
then
    echo "Everything works, hooray!"
else
    echo "The following runcards have failed:" ${failures}
    exit -1
fi
