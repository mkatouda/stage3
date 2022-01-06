#!/bin/bash

ligand=ethanol.mol
outbasename=etanol-wat
ff=cgenff #gaff
charge=am1bcc #am1bcc-pol, hf
water=tip3p

export OMP_NUM_THREADS=4

module load gcc/6.3.0
. /home/app/intel/intel2018_up3/bin/compilervars.sh intel64
. /home/app/gromacs/2020.6/cpu/bin/GMXRC.bash

export PerlChemistry=${HOME}/scripts/MATCH_RELEASE/PerlChemistry
export MATCH=${HOME}/scripts/MATCH_RELEASE/MATCH

export PATH=${PATH}:${MATCH}/scripts

stage3path=${HOME}/work/kyoto_pro-lig_md_2021q4/stage3/stage3
templatepath=${HOME}/work/kyoto_pro-lig_md_2021q4/stage3/templates


rm -rf ./${outbasename}_${ff}

python ${stage3path}/stage.py -i ${ligand} -o ${outbasename} -v --forcefields ${ff} -q ${charge} -w ${water}

cd ${outbasename}_${ff}

gmx_mpi grompp -f ${templatepath}/minim.mdp -c solvated.gro -r solvated.gro -p ${outbasename}.top -o min.tpr -po mdout.mdp

cd ../

