#!/bin/bash -e

. /home/tcj25/.virtualenvs/35/bin/activate

task=$1
out=$task-unmapped.fastq.gz
bam=../02-map/$task.bam
log=../$task.log

function unmapped()
{
    echo "  print-unmapped-sam.py started at `date`" >> $log
    print-unmapped-sam.py $bam | gzip > $out
    echo "  print-unmapped-sam.py stopped at `date`" >> $log
}

echo "03-find-unmapped on task $task started at `date`" >> $log
echo "  bam file is $bam" >> $log

if [ $SP_SIMULATE = "0" ]
then
    echo "  This is not a simulation." >> $log
    if [ -f $out ]
    then
        if [ $SP_FORCE = "1" ]
        then
            echo "  Pre-existing output file $out exists, but --force was used. Overwriting." >> $log
            unmapped
        else
            echo "  Will not overwrite pre-existing output file $out. Use --force to make me." >> $log
        fi
    else
        echo "  No pre-existing output file $out exists. Finding unmapped." >> $log
        unmapped
    fi
else
    echo "  This is a simulation." >> $log
fi

echo "03-find-unmapped on task $task stopped at `date`" >> $log
echo >> $log
