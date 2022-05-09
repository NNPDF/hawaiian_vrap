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
vrapout="results.out"

run_for_kinematics() {
    runcard=$1
    kinfile=$2
    outfile=$3

    echo " > Running vrap for $runcard"
    echo "# mass rapidity     xs" > ${outfile}

    if [[ $verbose == "y" ]]
    then
        ${vrap} ${runcard} ${kinfile}
    else
        echo "Running for ${kinfile}, this could take some minutes..."
        echo "If you want to follow the results do ~$ tail -f ${vrapout} in a different terminal"
        ${vrap} ${runcard} ${kinfile} >/dev/null
    fi
    cat ${vrapout} >> ${outfile}
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
