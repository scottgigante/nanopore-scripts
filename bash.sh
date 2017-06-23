#!/bin/bash
#PBS -q small 

# Usage: qsub -F "echo 'hello'\; echo 'world'" bash.sh

CMD="$*"
set -x
cd $PBS_O_WORKDIR
eval "$CMD"
set +x
echo "COMPLETE"
