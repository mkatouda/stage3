#!/bin/bash
#!/bin/bash
# ITO-A
# PJM -L rscunit=ito-a
# PJM -L rscgrp=ito-ss
# PJM -L rscgrp=ito-ss-dbg
# PJM -L vnode=1
# PJM -L vnode-core=36
# ITO-B
# PJM -L rscunit=ito-b
# PJM -L rscgrp=ito-g-4
# PJM -L rscgrp=ito-g-4-dbg
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
#PJM -L elapse=0:15:00
#PJM -j
#PJM -S

# MD input settings
inpbasename=3HTB-jz4-wat
ff_ligand=cgenff
maxwarn=10
mdppath=${HOME}/data/work/kyoto_pro-lig_md_2021q4/stage3/templates
#mdppath=../mdp

# computing system settings
system=flow-cx
export OMP_NUM_THREADS=10
LANG=C

# software envieonment settings
# Groamcs settings
GMX_CMD=gmx_mpi
#export GMXLIB=${HOME}/data/param/gromacs/top

if [ ${system,,} = 'ito-a' ]; then
    module load gcc/6.3.0
    . /home/app/intel/intel2018_up3/bin/compilervars.sh intel64
    . /home/app/gromacs/2020.6/cpu/bin/GMXRC.bash
    export PATH=${PATH}:/home/app/a/GAMESS-20170420
elif [ ${system,,} = 'ito-b' ]; then
    module load gcc/6.3.0
    . /home/app/intel/intel2018_up3/bin/compilervars.sh intel64
    . /home/app/gromacs/2020.6/gpu/bin/GMXRC.bash
    export PATH=${PATH}:/home/app/a/GAMESS-20170420
elif [ ${system,,} = 'flow-fx' ]; then
    module load fftw/3.3.9-tune
    module load gromacs/2021.2
elif [ ${system,,} = 'flow-cx' ]; then
    module load gcc/8.4.0 cuda/11.2.1 openmpi_cuda/4.0.4 nccl/2.8.4
    GMX_BIN=${HOME}/data/bin/x86_64/gromacs/2021.4/gpu/gcc/bin
    . ${GMX_BIN}/GMXRC.bash
elif [ ${system,,} = 'flow-cloud' ]; then
    module load gcc/8.4.0
    . ${HOME}/data/bin/x86_64/gromacs/2021.1/gcc/bin/GMXRC.bash
else
    GMX_BIN=${HOME}/gromacs/2021.4/bin
    . ${GMX_BIN}/GMXRC.bash
fi

cd ${outbasename}_${ff_ligand}
if [ -e solvated_ionised.gro ]:
    gro=solvated_ionised.gro
elif [ -e solvated.gro ]:
    gro=solvated.gro
else
    gro=${outbasebame}.gro
fi
${GMX_CMD} grompp -f ${mdppath}/minim.mdp -c ${gro} -r ${gro} -p ${inpbasename}.top -o min.tpr -po mdout.mdp -maxwarn ${maxwarn}
${GMX_CMD} mdrun -deffnm min
cd ../

