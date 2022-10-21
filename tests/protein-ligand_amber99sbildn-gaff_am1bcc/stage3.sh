#!/bin/bash

# MD input settings
ligand=jz4.mol
protein=3HTB_protein.pdb
outbasename=3HTB-jz4-wat
ff_ligand=gaff #gaff cgenff
charge_model=am1bcc #am1bcc am1bcc-pol hf 
ff_protein=amber99sb-ildn #amber99sb-ildn charmm27
water_model=tip3p
box_type=dodecahedron
box_buff=1.0
pname=NA
nname=CL
rdd=1
maxwarn=10
mdppath=`pwd`/mdp

# Gromacs settings
module load gcc/8.4.0
GMX_CMD=gmx
GMX_BIN=${HOME}/data/bin/x86_64/gromacs/2021.4/cpu/gcc/bin
. ${GMX_BIN}/GMXRC.bash

# MATCH settings
#export PerlChemistry=${HOME}/data/scripts/MATCH_RELEASE/PerlChemistry
#export MATCH=${HOME}/data/scripts/MATCH_RELEASE/MATCH
#export PATH=${PATH}:${MATCH}/scripts

if [ -d ./${outbasename}_${ff_ligand} ]; then
    rm -rf ./${outbasename}_${ff_ligand}
fi

# Run STaGE3
stage3 -l ${ligand} -o ${outbasename} --ffligand ${ff_ligand} -q ${charge_model} -c ${protein} -b ${box_type} -d ${box_buff} -w ${water_model} --ffprotein ${ff_protein} --pname ${pname} --nname ${nname} -v

cd ${outbasename}_${ff_ligand}

if [ -e solvated_ionised.gro ]; then
    gro=solvated_ionised.gro
elif [ -e solvated.gro ]; then
    gro=solvated.gro
else
    gro=${outbasebame}.gro
fi
${GMX_CMD} grompp -f ${mdppath}/em.mdp -c ${gro} -r ${gro} -p ${outbasename}.top -o min.tpr -po mdout.mdp -maxwarn ${maxwarn}

cd ../

