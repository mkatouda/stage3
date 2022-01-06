#!/bin/bash
#PJM -L rscunit=ito-a
# PJM -L rscgrp=ito-ss
#PJM -L rscgrp=ito-ss-dbg
#PJM -L vnode=1
#PJM -L vnode-core=36
# PJM -L rscunit=cl
# PJM -L rscgrp=cl-share
# PJM -L node=1
# PJM --mpi proc=20
#PJM -L elapse=0:15:00
#PJM -j
#PJM -S

inpbasename=3HTB-jz4-wat
ff=cgenff

export OMP_NUM_THREADS=4

module load gcc/6.3.0
. /home/app/intel/intel2018_up3/bin/compilervars.sh intel64
. /home/app/gromacs/2020.6/cpu/bin/GMXRC.bash

cd ${inpbasename}_${ff}

gmx_mpi mdrun -deffnm min

cd ../

