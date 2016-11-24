#!/bin/bash -e

log=../slurm-pipeline.log

echo "06-stop sbatch.sh running at `date`" >> $log
echo "  dependencies are $SP_DEPENDENCY_ARG" >> $log
for task in "$@"
do
    echo "  task $task" >> $log
done

echo >> $log

jobid=`sbatch -n 1 $SP_DEPENDENCY_ARG stop.sh "$@" | cut -f4 -d' '`
echo "TASK: stop $jobid"
