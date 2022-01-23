#!/bin/bash
# ITO-A
# PJM -L rscunit=ito-a
# PJM -L rscgrp=ito-ss
# PJM -L rscgrp=ito-ss-dbg
# PJM -L vnode=1
# PJM -L vnode-core=36
# ITO-B
# PJM -L rscunit=ito-b
# PJM -L rscgrp=ito-g-1
# PJM -L vnode=1
# PJM -L vnode-core=36
# Flow-fx
# PJM -L rscunit=fx
# PJM -L rscgrp=fx-debug
# PJM -L node=1
# PJM --mpi proc=4
# Flow-cx
#PJM -L rscunit=cx
#PJM -L rscgrp=cx-share
#PJM -L gpu=1
# Flow-cloud
# PJM -L rscunit=cl
# PJM -L rscgrp=cl-share
# PJM -L node=1
# PJM --mpi proc=20
# Common settings
#PJM -L elapse=1:00:00
#PJM -j
#PJM -S

# MD input settings
ligand=jz4.mol
protein=3HTB_protein.pdb
outbasename=3HTB-jz4-wat
ff_ligand=cgenff #cgenff gaff
charge_model=am1bcc #am1bcc am1bcc-pol hf 
ff_protein=charmm36-feb2021 #charmm27 charmm36-feb2021 amber99sb-ildn
water_model=tip3p
box_type=dodecahedron
box_buff=1.0
pname=NA
nname=CL
rdd=1
maxwarn=10
#mdppath=${HOME}/data/work/kyoto_pro-lig_md_2021q4/stage3/templates
mdppath=../mdp

# computing system settings
system=flow-cx # flow-fx flow-cloud ito-a ito-b
conda_venv=py38-shortmd
quesystem=none
USE_GPU=false
GPU_ID="0"
MPI_lib=none #intel OpenMPI
NUM_NODES=1
NUM_CORES=4
NUM_THREADS=4
LANG=C

# software envieonment settings
# Groamcs settings
GMX_CMD=gmx
GMX_BIN=${HOME}/data/gromacs-2021.4/bin
export GMXLIB=${HOME}/data/param/gromacs/top

# computing system dependent settings

if [ ${system,,} = 'ito-a' ]; then
    module load gcc/6.3.0
    . /home/app/intel/intel2018_up3/bin/compilervars.sh intel64
    GMX_CMD=gmx_mpi
    GMX_BIN=/home/app/gromacs/2020.6/cpu/bin
    . ${GMX_BIN}/GMXRC.bash
    module load g16/a03
    export PATH=${PATH}:/home/app/a/GAMESS-20170420
    quesystem=pjm
    MPI_lib=intel
    USE_GPU=false
    NUM_NODES=1
    NUM_CORES=36
    NUM_THREADS=36
elif [ ${system,,} = 'ito-b' ]; then
    module load gcc/6.3.0 cuda/9.1
    . /home/app/intel/intel2018_up3/bin/compilervars.sh intel64
    GMX_BIN=/home/app/gromacs/2020.6/gpu/bin
    GMX_CMD=gmx_mpi
    . ${GMX_BIN}/GMXRC.bash
    module load g16/a03
    export PATH=${PATH}:/home/app/a/GAMESS-20170420
    quesystem=pjm
    USE_GPU=true
    GPU_ID="0"
    NUM_NODES=1
    NUM_CORES=9
    NUM_THREADS=9
elif [ ${system,,} = 'flow-fx' ]; then
    module load fftw/3.3.9-tune
    module load gromacs/2021.2
    GMX_CMD=gmx_mpi
    quesystem=pjm
    MPI_lib=tcs
    USE_GPU=false
    NUM_NODES=1
    NUM_CORES=40
    NUM_THREADS=10
elif [ ${system,,} = 'flow-cx' ]; then
    module load gcc/8.4.0 cuda/11.2.1 openmpi_cuda/4.0.4 nccl/2.8.4
    GMX_CMD=gmx
    GMX_BIN=${HOME}/data/bin/x86_64/gromacs/2021.4/gpu/gcc/bin
    . ${GMX_BIN}/GMXRC.bash
    quesystem=pjm
    MPI_lib=none
    USE_GPU=true
    GPU_ID="0"
    NUM_NODES=1
    NUM_CORES=10
    NUM_THREADS=10
elif [ ${system,,} = 'flow-cloud' ]; then
    module load gcc/8.4.0
    GMX_CMD=gmx
    GMX_BIN=${HOME}/data/bin/x86_64/gromacs/2021.4/gcc/bin
    . ${GMX_BIN}/GMXRC.bash
    quesystem=pjm
    MPI_lib=OpenMPI
    USE_GPU=false
    NUM_NODES=1
    NUM_CORES=20
    NUM_THREADS=20
else
    . ${GMX_BIN}/GMXRC.bash
fi

# MATCH settings
export PerlChemistry=${HOME}/data/scripts/MATCH_RELEASE/PerlChemistry
export MATCH=${HOME}/data/scripts/MATCH_RELEASE/MATCH
export PATH=${PATH}:${MATCH}/scripts

# STaGE3 settings
stage3path=${HOME}/data/work/kyoto_pro-lig_md_2021q4/stage3/stage3

if [ -d ./${outbasename}_${ff_ligand} ]; then
    rm -rf ./${outbasename}_${ff_ligand}
fi

# Run STaGE3
. ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate ${conda_venv}

python ${stage3path}/stage.py -i ${ligand} -o ${outbasename} --forcefields ${ff_ligand} -q ${charge_model} -c ${protein} -b ${box_type} -d ${box_buff} -w ${water_model} --ffprotein ${ff_protein} --pname ${pname} --nname ${nname} -v

conda deactivate
conda deactivate

# Run gromacs
if [ ${MPI_lib,,} = 'intel' ]; then
    export I_MPI_PERHOST=`expr ${NUM_CORES} / ${NUM_THREADS}`
    export I_MPI_FABRICS=shm:ofi
    export I_MPI_PIN_DOMAIN=omp
    export I_MPI_PIN_CELL=core
    export I_MPI_HYDRA_BOOTSTRAP=rsh
    export I_MPI_HYDRA_BOOTSTRAP_EXEC=/bin/pjrsh
    export I_MPI_HYDRA_HOST_FILE=${PJM_O_NODEINF}
    MPIRUN=mpiexec.hydra
elif [ ${MPI_lib,,} = 'openmpi' ]; then
    MPIRUN=mpirun
elif [ ${MPI_lib,,} = 'tcs' ]; then
    MPIRUN=mpiexec
else
    MPIRUN=''
fi

#if [ ${quesystem,,} = 'pjm' ]; then
#    NUM_NODES=${PJM_VNODES}
#fi
NUM_PROCS=`expr ${NUM_NODES} "*" ${NUM_CORES} / ${NUM_THREADS}`
echo NUM_NODES=${NUM_NODES} NUM_CORES=${NUM_CORES} NUM_PROCS=${NUM_PROCS} NUM_THREADS=${NUM_THREADS}
export OMP_NUM_THREADS=${NUM_THREADS}
#export KMP_STACKSIZE=8m
#export KMP_AFFINITY=compact

if [ ${USE_GPU,,} = 'true' ]; then
    #export GMX_GPU_DD_COMMS=true
    #export GMX_GPU_PME_PP_COMMS=true
    export GMX_FORCE_UPDATE_DEFAULT_GPU=true
fi

gmx_run_md() {
    local mdp=$1
    local conf=$2
    local restraint=$3
    local tpr=$4
    if [ -n "${MPIRUN}" ]; then
	if [ -f ${ndx} ]; then
	    ${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} grompp -f ${mdp} -c ${conf} -r ${restraint} -p ${top} -n ${ndx} -o ${tpr} -maxwarn ${maxwarn}
	else
	    ${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} grompp -f ${mdp} -c ${conf} -r ${restraint} -p ${top} -o ${tpr} -maxwarn ${maxwarn}
	fi
	if [ ${USE_GPU,,} = 'true' ]; then
	    ${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} mdrun -ntomp ${NUM_THREADS} -rdd ${rdd} -nb gpu -gpu_id ${GPU_ID} -deffnm ${tpr%.*}
	elif [ ${NUM_PROCS} -gt 1 ]; then
	    ${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} mdrun -ntomp ${NUM_THREADS} -rdd ${rdd} -deffnm ${tpr%.*}
	fi
    else
	if [ -f ${ndx} ]; then
	    ${GMX_CMD} grompp -f ${mdp} -c ${conf} -r ${restraint} -p ${top} -n ${ndx} -o ${tpr} -maxwarn ${maxwarn}
	else
	    ${GMX_CMD} grompp -f ${mdp} -c ${conf} -r ${restraint} -p ${top} -o ${tpr} -maxwarn ${maxwarn}
	fi
	if [ ${USE_GPU,,} = 'true' ]; then
	    ${GMX_CMD} mdrun -ntmpi ${NUM_PROCS} -ntomp ${NUM_THREADS} -nb gpu -gpu_id ${GPU_ID} -deffnm ${tpr%.*}
	else
	    ${GMX_CMD} mdrun -ntmpi ${NUM_PROCS} -ntomp ${NUM_THREADS} -deffnm ${tpr%.*}
	fi
    fi
}

gmx_run_md_restart() {
    local mdp=$1
    local conf=$2
    local restraint=$3
    local cpt=$4
    local tpr=$5
    if [ -n "${MPIRUN}" ]; then
	if [ -f ${ndx} ]; then
	    ${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} grompp -f ${mdp} -c ${conf} -r ${restraint} -t ${cpt} -p ${top} -n ${ndx} -o ${tpr} -maxwarn ${maxwarn}
	else
	    ${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} grompp -f ${mdp} -c ${conf} -r ${restraint} -t ${cpt} -p ${top} -o ${tpr} -maxwarn ${maxwarn}
	fi
	if [ ${USE_GPU,,} = 'true' ]; then
	    ${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} mdrun -ntomp ${NUM_THREADS} -rdd ${rdd} -nb gpu -gpu_id ${GPU_ID} -deffnm ${tpr%.*}
	elif [ ${NUM_PROCS} -gt 1 ]; then
	    ${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} mdrun -ntomp ${NUM_THREADS} -rdd ${rdd} -deffnm ${tpr%.*}
	fi
    else
	if [ -f ${ndx} ]; then
	    ${GMX_CMD} grompp -f ${mdp} -c ${conf} -r ${restraint} -t ${cpt} -p ${top} -n ${ndx} -o ${tpr} -maxwarn ${maxwarn}
	else
	    ${GMX_CMD} grompp -f ${mdp} -c ${conf} -r ${restraint} -t ${cpt} -p ${top} -o ${tpr} -maxwarn ${maxwarn}
	fi
	if [ ${USE_GPU,,} = 'true' ]; then
	    ${GMX_CMD} mdrun -ntmpi ${NUM_PROCS} -ntomp ${NUM_THREADS} -nb gpu -gpu_id ${GPU_ID} -deffnm ${tpr%.*}
	else
	    ${GMX_CMD} mdrun -ntmpi ${NUM_PROCS} -ntomp ${NUM_THREADS} -deffnm ${tpr%.*}
	fi
    fi
}

gmx_rerun_md() {
    local mdp=$1
    local conf=$2
    local restraint=$3
    local cpt=$4
    local tpr=$5
    local xtc=$6
    if [ -n "${MPIRUN}" ]; then
	if [ -f ${ndx} ]; then
	    ${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} grompp -f ${mdp} -c ${conf} -r ${restraint} -t ${cpt} -p ${top} -n ${ndx} -o ${tpr} -maxwarn ${maxwarn}
	else
	    ${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} grompp -f ${mdp} -c ${conf} -r ${restraint} -t ${cpt} -p ${top} -o ${tpr} -maxwarn ${maxwarn}
	fi
	${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} mdrun -ntomp ${NUM_THREADS} -rdd ${rdd} -deffnm ${tpr%.*} -rerun ${xtc}
    else
	if [ -f ${ndx} ]; then
	    ${GMX_CMD} grompp -f ${mdp} -c ${conf} -r ${restraint} -t ${cpt} -p ${top} -n ${ndx} -o ${tpr} -maxwarn ${maxwarn}
	else
	    ${GMX_CMD} grompp -f ${mdp} -c ${conf} -r ${restraint} -t ${cpt} -p ${top} -o ${tpr} -maxwarn ${maxwarn}
	fi
	${GMX_CMD} mdrun -ntmpi ${NUM_PROCS} -ntomp ${NUM_THREADS} -deffnm ${tpr%.*} -rerun ${xtc}
    fi
}

gmx_energy_intr() {
    local edr=$1
    local xvg_comp=$2
    local xvg_sum=$3
    if [ -n "${MPIRUN}" ]; then
	${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} energy -f ${edr} -o ${xvg_comp} <<EOF
Coul-SR:Protein-LIG
LJ-SR:Protein-LIG
0
EOF
	${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} energy -f ${edr} -o ${xvg_sum} -sum <<EOF
Coul-SR:Protein-LIG
LJ-SR:Protein-LIG
0
EOF
    else
	${GMX_CMD} energy -f ${edr} -o ${xvg_comp} <<EOF
Coul-SR:Protein-LIG
LJ-SR:Protein-LIG
0
EOF
	${GMX_CMD} energy -f ${edr} -o ${xvg_sum} -sum <<EOF
Coul-SR:Protein-LIG
LJ-SR:Protein-LIG
0
EOF
    fi
}

gmx_trjconv() {
    local tpr=$1
    local xtc=$2
    local xtc_center=${xtc%.*}_center.xtc
    local xtc_fit=${xtc%.*}_fit.xtc
    local xtc_fit_nowat=${xtc%.*}_fit_nowat.xtc
    #local gro_center_start=${xtc%.*}_center_start_nowat.pdb
    local gro_fit_nowat_start=${xtc%.*}_fit_nowat_start.gro
    if [ -n "${MPIRUN}" ]; then
	${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} trjconv -s ${tpr} -f ${xtc} -o ${xtc_center} -center -pbc mol -ur compact <<EOF
Protein
System
EOF
        #${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} trjconv -s ${tpr} -f ${xtc_center} -o ${gro_center_start} -dump 0
        ${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} trjconv -s ${tpr} -f ${xtc_center} -o ${xtc_fit} -fit rot+trans <<EOF
Backbone
System
EOF
        ${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} trjconv -s ${tpr} -f ${xtc_center} -n ${ndx} -o ${gro_fit_nowat_start} -dump 0 <<EOF
TempCouplingA
EOF
        ${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} trjconv -s ${tpr} -f ${xtc_fit} -n ${ndx} -o ${xtc_fit_nowat} <<EOF
TempCouplingA
EOF
    else
	${GMX_CMD} trjconv -s ${tpr} -f ${xtc} -o ${xtc_center} -center -pbc mol -ur compact <<EOF
Protein
System
EOF
        #${GMX_CMD} trjconv -s ${tpr} -f ${xtc_center} -o ${gro_center_start} -dump 0
        ${GMX_CMD} trjconv -s ${tpr} -f ${xtc_center} -o ${xtc_fit} -fit rot+trans <<EOF
Backbone
System
EOF
        ${GMX_CMD} trjconv -s ${tpr} -f ${xtc_center} -n ${ndx} -o ${gro_fit_nowat_start} -dump 0 <<EOF
TempCouplingA
EOF
        ${GMX_CMD} trjconv -s ${tpr} -f ${xtc_fit} -n ${ndx} -o ${xtc_fit_nowat} <<EOF
TempCouplingA
EOF
    fi
}

gmx_rms() {
    local tpr=$1
    local xtc=$2
    local xvg=$3
    if [ -n "${MPIRUN}" ]; then
	${MPIRUN} -n ${NUM_PROCS} ${GMX_CMD} rms -s ${tpr} -f ${xtc} -n ${ndx} -o ${xvg} -tu ns <<EOF
Backbone
2
EOF
    else
	${GMX_CMD} rms -s ${tpr} -f ${xtc} -n ${ndx} -o ${xvg} -tu ns <<EOF
Backbone
2
EOF
    fi
}

top=${outbasename}.top
ndx=index.ndx

cd ${outbasename}_${ff_ligand}

if [ -e solvated_ionised.gro ]; then
    gro=solvated_ionised.gro
elif [ -e solvated.gro ]; then
    gro=solvated.gro
else
    gro=${outbasename}.gro
fi

gmx_run_md ${mdppath}/em.mdp ${gro} ${gro} em.tpr

gmx_run_md ${mdppath}/nvt.mdp em.gro em.gro nvt.tpr

gmx_run_md_restart ${mdppath}/npt.mdp nvt.gro nvt.gro nvt.cpt npt.tpr

gmx_run_md_restart ${mdppath}/md.mdp npt.gro npt.gro npt.cpt md_0_10.tpr

gmx_rerun_md ${mdppath}/ie.mdp npt.gro npt.gro npt.cpt ie.tpr md_0_10.xtc

gmx_energy_intr ie.edr md_0_10_interaction_energy_component.xvg md_0_10_interaction_energy_sum.xvg

gmx_trjconv md_0_10.tpr md_0_10.xtc

gmx_rms md_0_10.tpr md_0_10_center.xtc md_0_10_ligand_rmsd.xvg

cd ../

exit 0
