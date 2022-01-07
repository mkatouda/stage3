#!/bin/bash

ligand=jz4.mol
protein=3HTB_protein.pdb
outbasename=3HTB-jz4-wat
ff=cgenff #cgenff #gaff
charge=am1bcc #am1bcc am1bcc-pol hf 
ffprotein=charmm27 #charmm27 charmm36-feb2021 amber99sb-ildn
water=tip3p

export OMP_NUM_THREADS=4

module load gcc/6.3.0

. /home/app/intel/intel2018_up3/bin/compilervars.sh intel64
. /home/app/gromacs/2020.6/cpu/bin/GMXRC.bash

#module load g16/a03

export PerlChemistry=${HOME}/scripts/MATCH_RELEASE/PerlChemistry
export MATCH=${HOME}/scripts/MATCH_RELEASE/MATCH

export PATH=${PATH}:${MATCH}/scripts

#export PATH=${PATH}:/home/app/a/GAMESS-20170420

export GMXLIB=${HOME}/param/gromacs/top

stage3path=${HOME}/work/kyoto_pro-lig_md_2021q4/stage3/stage3
templatepath=${HOME}/work/kyoto_pro-lig_md_2021q4/stage3/templates

rm -rf ./${outbasename}_${ff}

python ${stage3path}/stage.py -i ${ligand} -o ${outbasename} -v --forcefields ${ff} -q ${charge} -c ${protein} -w ${water} --ffprotein ${ffprotein}
#python ${stage3path}/stage.py -i ${ligand} -o ${outbasename} -v --forcefields ${ff} -q ${charge} -c ${protein} -w ${water} --ffprotein ${ffprotein} --pname SOD --nname CLA

cd ${outbasename}_${ff}

gmx_mpi grompp -f ${templatepath}/minim.mdp -c solvated_ionised.gro -r solvated_ionised.gro -p ${outbasename}.top -o min.tpr -po mdout.mdp -maxwarn 10

cd ../

