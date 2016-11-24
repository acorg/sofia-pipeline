#!/bin/bash -e

#SBATCH -J trim
#SBATCH -A DSMITH-BIOCLOUD
#SBATCH -o slurm-%A.out
#SBATCH -p biocloud-normal
#SBATCH --time=05:00:00

. /home/tcj25/.virtualenvs/35/bin/activate

task=$1
fastq=../$task.fastq.gz
log=../$task.log
out=$task-trimmed.fastq.gz

echo "01-trim on task $task started at `date`" >> $log
echo "  fastq is $fastq" >> $log


function trim()
{
    echo "  No trimming needed." >> $log
}


function original_trim()
{
    # The adapter is either in field 7 or 8 (underscore separated) and must
    # have 6 characters from the set A, C, G, T, N (this is what
    # AdapterRemoval says when you give it something that it cannot
    # handle).
    #
    # We also temporarily allow commands to exit non-zero, as the greps
    # below may fail.
    set +e
    adapter=`echo $task | cut -f7 -d_ | egrep '^[ACGTN]{6}$'`

    if [ -z "$adapter" ]
    then
        adapter=`echo $task | cut -f8 -d_ | egrep '^[ACGTN]{6}$'`

        if [ -z "$adapter" ]
        then
            echo "  WARNING: Could not find adapter sequence in task $task" >> $log
            exit 1
        fi
    fi
    set -e

    echo "  Adapter is $adapter" >> $log

    echo "  AdapterRemoval started at `date`" >> $log
    srun -n 1 AdapterRemoval \
      --basename $task \
      --adapter1 $adapter \
      --file1 $fastq \
      --output1 $out \
      --gzip \
      --trimns \
      --trimqualities
    echo "  AdapterRemoval stopped at `date`" >> $log
}


if [ $SP_SIMULATE = "0" ]
then
    echo "  This is not a simulation." >> $log
    if [ -f $out ]
    then
        if [ $SP_FORCE = "1" ]
        then
            echo "  Pre-existing output file $out exists, but --force was used. Overwriting." >> $log
            trim
        else
            echo "  Will not overwrite pre-existing output file $out. Use --force to make me." >> $log
        fi
    else
        echo "  No pre-existing output file $out exists. Trimming." >> $log
        trim
    fi
else
    echo "  This is a simulation." >> $log
fi

echo "01-trim on task $task stopped at `date`" >> $log
echo >> $log
