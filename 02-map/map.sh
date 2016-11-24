#!/bin/bash -e

. /home/tcj25/.virtualenvs/35/bin/activate

task=$1
out=$task.bam
bamtmp=/ramdisks/terry-$task.bam
sam=/ramdisks/terry-$task.sam
bwadb=$HOME/scratch/homo-sapiens/homo-sapiens
fastq=../01-trim/$task-trimmed.fastq.gz
log=../$task.log

function map()
{
    # Map to human genome, save to ramdisk SAM file.
    echo "  bwa mem started at `date`" >> $log
    bwa mem -t 24 $bwadb $fastq > $sam
    echo "  bwa mem stopped at `date`" >> $log

    # Convert SAM to BAM, on ramdisk.
    echo "  sam -> bam conversion started at `date`" >> $log
    samtools view --threads 24 -bS $sam > $bamtmp
    echo "  sam -> bam conversion stopped at `date`" >> $log

    # Sort ramdisk BAM, writing to persisted local storage.
    echo "  bam sort started at `date`" >> $log
    samtools sort --threads 24 -o $out $bamtmp
    echo "  bam sort stopped at `date`" >> $log
}


echo "02-map on task $task started at `date`" >> $log
echo "  fastq is $fastq" >> $log

if [ $SP_SIMULATE = "0" ]
then
    echo "  This is not a simulation." >> $log
    if [ -f $out ]
    then
        if [ $SP_FORCE = "1" ]
        then
            echo "  Pre-existing output file $out exists, but --force was used. Overwriting." >> $log
            map
        else
            echo "  Will not overwrite pre-existing output file $out. Use --force to make me." >> $log
        fi
    else
        echo "  No pre-existing output file $out exists. Mapping." >> $log
        map
    fi
else
    echo "  This is a simulation." >> $log
fi

echo "02-map on task $task stopped at `date`" >> $log
echo >> $log
