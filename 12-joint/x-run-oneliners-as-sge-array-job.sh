#!/bin/bash
#$ -S /bin/bash

if [ "x${SGE_TASK_STEPSIZE}" == "x" -o "$SGE_TASK_STEPSIZE" == "undefined" ]; then
    export SGE_TASK_STEPSIZE=1
fi

echo "# SGE_TASK_ID = ${SGE_TASK_ID}" 2>&1
exec sh -c "$(cat $1 | awk -v id=$SGE_TASK_ID -v size=$SGE_TASK_STEPSIZE 'id <= NR && NR < id + size { print $0 }')"

