#!/bin/bash -e

. /home/tcj25/.virtualenvs/35/bin/activate

task=$1
log=../slurm-pipeline.log

fastq1=../02-flash/${task}_R1_001-flash.fastq.gz
fastq2=../02-flash/${task}_R2_001-flash.fastq.gz
fastqMerged=../02-flash/${task}-merged-flash.fastq.gz

out=$task-unmapped.fastq.gz
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

function skip()
{
    # Copy our input FASTQ to our output (almost) unchanged. We need put /1
    # at the end of the ids in $fastq1 and /2 at the end of $fastq2,
    # otherwise we'll end up with duplicate read ids from the paired reads.

    # In the following, filter-fasta.py is used to make sure sequences and
    # quality strings occupy just one line, so that awk can append /1 or /2
    # to the FASTQ id lines (matching those lines via NR % 4 == 1).
    (
        zcat $fastq1 | filter-fasta.py --readClass fastq | awk '{if (NR % 2 == 1){printf "%s/1\n", $0} else {print}}'
        zcat $fastq2 | filter-fasta.py --readClass fastq | awk '{if (NR % 2 == 1){printf "%s/1\n", $0} else {print}}'
    ) | gzip > $out

    cat $fastqMerged >> $out
}

function get_genome()
{
    # Look for the sample id in ../../../../sample-to-genome and see if it
    # has a genome that we should map to.

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

function find_unmapped()
{
    echo "  print-unmapped-sam.py started on $bamPairs at `date`" >> $log
    print-unmapped-sam.py $outPairs | gzip > $out
    echo "  print-unmapped-sam.py stopped on $bamPairs at `date`" >> $log

    # Note we concatenate onto the end of the previous gzip output.
    echo "  print-unmapped-sam.py started on $bamMerged at `date`" >> $log
    print-unmapped-sam.py $outMerged | gzip >> $out
    echo "  print-unmapped-sam.py stopped on $bamMerged at `date`" >> $log
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

if [ $SP_SIMULATE = "1" ]
then
    echo "  This is a simulation." >> $log
else
    echo "  This is not a simulation." >> $log
    if [ $SP_SKIP = "1" ]
    then
        echo "  Mapping is being skipped on this run." >> $log
        skip
    elif [ -f $out ]
    then
        if [ $SP_FORCE = "1" ]
        then
            echo "  Pre-existing output file $out exists, but --force was used. Overwriting." >> $log
            map
            find_unmapped
        else
            echo "  Will not overwrite pre-existing output file $out. Use --force to make me." >> $log
        fi
    else
        echo "  Pre-existing output file $out does not exist. Mapping." >> $log
        map
        find_unmapped
    fi
fi

echo "03-map on task $task stopped at `date`" >> $log
echo >> $log
