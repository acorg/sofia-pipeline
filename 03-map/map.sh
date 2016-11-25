#!/bin/bash -e

. /home/tcj25/.virtualenvs/35/bin/activate

task=$1
log=../slurm-pipeline.log

fastq1=../02-flash/${task}_R1_001-flash.fastq.gz
fastq2=../02-flash/${task}_R2_001-flash.fastq.gz
fastqMerged=../02-flash/${task}-merged-flash.fastq.gz

outPairs=$task-paired.bam
outMerged=$task-merged.bam

if [ -w /ramdisks ]
then
    tmp=/ramdisks
else
    tmp=/tmp
fi

bamtmp=$tmp/terry-$task.bam
sam=$tmp/terry-$task.sam

function get_genome()
{
    sample=$(basename $(dirname $(/bin/pwd)))
    sample_to_genome='../../../../sample-to-genome'

    if [ ! -f $sample_to_genome ]
    then
        echo "  Sample to genome mapping file '$sample_to_genome' does not exist." >> $log
        exit 1
    fi

    genome=$(egrep "^$sample " $sample_to_genome | awk '{print $2}')

    if [ -z "$genome" ]
    then
        echo "  Could not find sample '$sample' in genome mapping file '$sample_to_genome'." >> $log
        exit 1
    fi
    
    genomeDir=$HOME/scratch/genomes/$genome

    if [ ! -d $genomeDir ]
    then
        echo "  Genome directory '$genomeDir' does not exist." >> $log
        exit 1
    fi

    # BWA takes the basename of the index.
    genomeForBWA=$genomeDir/$genome

    # Make sure at least one of the actual BWA index files exists.
    if [ ! -f $genomeForBWA.amb ]
    then
        echo "  Genome index file '$genomeForBWA.amb' does not exist." >> $log
        exit 1
    fi

    echo $genomeForBWA
}

function map()
{
    genome=`get_genome`
    echo "  Genome base name is '$genome'." >> $log

    # Map paired FASTQ to human genome, save to ramdisk SAM file.
    echo "  bwa mem on pairs started at `date`" >> $log
    bwa mem -t 24 $genome $fastq1 $fastq2 > $sam
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
    bwa mem -t 24 $genome $fastqMerged > $sam
    echo "  bwa mem on (flash) merged stopped at `date`" >> $log

    # Convert SAM to BAM, on ramdisk.
    echo "  sam -> bam on (flash) merged conversion started at `date`" >> $log
    samtools view --threads 24 -bS $sam > $bamtmp
    echo "  sam -> bam on (flash) merged conversion stopped at `date`" >> $log

    # Sort ramdisk BAM, writing to persisted local storage.
    echo "  bam sort on (flash) merged started at `date`" >> $log
    samtools sort --threads 24 -o $outMerged $bamtmp
    echo "  bam sort on (flash) merged stopped at `date`" >> $log

    echo "  Sleeping 5 seconds to give the filesystem a chance to settle." >> $log
    sleep 5
}


echo "03-map on task $task started at `date`" >> $log
echo "  fastq1 is $fastq1" >> $log
echo "  fastq2 is $fastq2" >> $log
echo "  merged is $fastqMerged" >> $log

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
