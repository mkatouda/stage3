#!/bin/bash

ligand=ethanol.mol
outbasename=ethanol-wat
ffligand=gaff #gaff cgenff
charge=am1bcc #hf, am1bcc-pol
water=tip3p

export OMP_NUM_THREADS=4

module load gcc/8.4.0
GMX_CMD=gmx
GMX_BIN=${HOME}/data/bin/x86_64/gromacs/2021.4/cpu/gcc/bin
. ${GMX_BIN}/GMXRC.bash

#export PerlChemistry=${HOME}/scripts/MATCH_RELEASE/PerlChemistry
#export MATCH=${HOME}/scripts/MATCH_RELEASE/MATCH
#export PATH=${PATH}:${MATCH}/scripts
#export PATH=${PATH}:/home/app/a/GAMESS-20170420
#export GMXLIB=${HOME}/param/gromacs/top

#templatepath=${HOME}/work/kyoto_pro-lig_md_2021q4/stage3/templates

rm -rf ./${outbasename}_${ff}

stage3 -l ${ligand} -o ${outbasename} -v --ffligand ${ffligand} -q ${charge} -w ${water}

#cd ${outbasename}_${ff}

#${GMX_CMD} grompp -f ${templatepath}/minim.mdp -c solvated.gro -r solvated.gro -p ${outbasename}.top -o min.tpr -po mdout.mdp -maxwarn 10

#cd ../

