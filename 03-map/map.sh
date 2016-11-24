#!/bin/bash -e

. /home/tcj25/.virtualenvs/35/bin/activate

task=$1

fastq1=../02-flash/${task}_R1_001-flash.fastq.gz
fastq2=../02-flash/${task}_R2_001-flash.fastq.gz
fastqMerged=../02-flash/${task}-merged-flash.fastq.gz

outPairs=$task-paired.bam
outMerged=$task-merged.bam

bamtmp=/ramdisks/terry-$task.bam
sam=/ramdisks/terry-$task.sam

bwadb=$HOME/scratch/homo-sapiens/homo-sapiens
log=../$task.log

function map()
{
    # Map paired FASTQ to human genome, save to ramdisk SAM file.
    echo "  bwa mem on pairs started at `date`" >> $log
    bwa mem -t 24 $bwadb $fastq1 $fastq2 > $sam
    echo "  bwa mem on pairs stopped at `date`" >> $log

    # Convert SAM to BAM, on ramdisk.
    echo "  sam -> bam on pairs conversion started at `date`" >> $log
    samtools view --threads 24 -bS $sam > $bamtmp
    echo "  sam -> bam on pairs conversion stopped at `date`" >> $log

    # Sort ramdisk BAM, writing to persisted local storage.
    echo "  bam sort on pairs started at `date`" >> $log
    samtools sort --threads 24 -o $outPairs $bamtmp
    echo "  bam sort on pairs stopped at `date`" >> $log

    # Map (flash) merged to human genome, save to ramdisk SAM file.
    echo "  bwa mem on (flash) merged started at `date`" >> $log
    bwa mem -t 24 $bwadb $fastqMerged > $sam
    echo "  bwa mem on (flash) merged stopped at `date`" >> $log

    # Convert SAM to BAM, on ramdisk.
    echo "  sam -> bam on (flash) merged conversion started at `date`" >> $log
    samtools view --threads 24 -bS $sam > $bamtmp
    echo "  sam -> bam on (flash) merged conversion stopped at `date`" >> $log

    # Sort ramdisk BAM, writing to persisted local storage.
    echo "  bam sort on (flash) merged started at `date`" >> $log
    samtools sort --threads 24 -o $outMerged $bamtmp
    echo "  bam sort on (flash) merged stopped at `date`" >> $log
}


echo "03-map on task $task started at `date`" >> $log
echo "  fastq is $fastq" >> $log

if [ $SP_SIMULATE = "0" ]
then
    echo "  This is not a simulation." >> $log
    if [ -f $outPairs -a -f $outMerged ]
    then
        if [ $SP_FORCE = "1" ]
        then
            echo "  Pre-existing output files $outPairs and $outMerged exist, but --force was used. Overwriting." >> $log
            map
        else
            echo "  Will not overwrite pre-existing output files $outPairs and $outMerged. Use --force to make me." >> $log
        fi
    else
        echo "  Pre-existing output files $outPairs and $outMerged do not both exist. Mapping." >> $log
        map
    fi
else
    echo "  This is a simulation." >> $log
fi

echo "03-map on task $task stopped at `date`" >> $log
echo >> $log
