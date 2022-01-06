#!/bin/bash

ligand=ethanol.mol
outbasename=ethanol
ff=cgenff
charge=am1bcc #hf, am1bcc-pol

export OMP_NUM_THREADS=4

module load gcc/6.3.0

. /home/app/intel/intel2018_up3/bin/compilervars.sh intel64
. /home/app/gromacs/2020.6/cpu/bin/GMXRC.bash

module load g16/a03

export PerlChemistry=${HOME}/scripts/MATCH_RELEASE/PerlChemistry
export MATCH=${HOME}/scripts/MATCH_RELEASE/MATCH

export PATH=${PATH}:${MATCH}/scripts

export PATH=${PATH}:/home/app/a/GAMESS-20170420

stage3path=${HOME}/work/kyoto_pro-lig_md_2021q4/stage3/stage3
templatepath=${HOME}/work/kyoto_pro-lig_md_2021q4/stage3/templates


rm -rf ./${outbasename}_${ff}

python ${stage3path}/stage.py -i ${ligand} -o ${outbasename} -v --forcefields ${ff} -q ${charge}

#cd ${outbasename}_${ff}

#gmx_mpi grompp -f ${templatepath}/premin_template.mdp -c ../${outbasename}.gro -r ../${outbasename}.gro -p ${outbasename}.top -o min.tpr -po mdout.mdp

#cd ../
