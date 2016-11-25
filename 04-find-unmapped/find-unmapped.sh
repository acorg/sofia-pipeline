#!/bin/bash -e

. /home/tcj25/.virtualenvs/35/bin/activate

task=$1
out=$task-unmapped.fastq.gz

bamPairs=../03-map/$task-paired.bam
bamMerged=../03-map/$task-merged.bam

log=../$task.log

function unmapped()
{
    echo "  print-unmapped-sam.py started on $bamPairs at `date`" >> $log
    print-unmapped-sam.py $bamPairs | gzip > $out
    echo "  print-unmapped-sam.py stopped on $bamPairs at `date`" >> $log

    # Note we concatenate onto the end of the previous gzip output.
    echo "  print-unmapped-sam.py started on $bamMerged at `date`" >> $log
    print-unmapped-sam.py $bamMerged | gzip >> $out
    echo "  print-unmapped-sam.py stopped on $bamMerged at `date`" >> $log

    echo "  Sleeping 30 seconds to give the filesystem a chance to settle."
    sleep 30
}

echo "04-find-unmapped on task $task started at `date`" >> $log
echo "  bam pairs file is $bamPairs" >> $log
echo "  bam merged file is $bamMerged" >> $log

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

echo "04-find-unmapped on task $task stopped at `date`" >> $log
echo >> $log
