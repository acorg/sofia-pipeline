#!/bin/bash -e

. /home/tcj25/.virtualenvs/35/bin/activate

task=$1
log=../slurm-pipeline.log
fastq=../03-map/$task-unmapped.fastq.gz
out=$task.json.bz2

echo "04-diamond on task $task started at `date`" >> $log
echo "  fastq file is $fastq" >> $log


dbfile=$HOME/scratch/root/share/ncbi/diamond-dbs/viral-protein-OKIAV-ECH.dmnd

if [ ! -f $dbfile ]
then
    echo "  DIAMOND database file $dbfile does not exist!" >> $log
    exit 1
fi

function run_diamond()
{
    echo "  DIAMOND blastx started at `date`" >> $log
    diamond blastx \
        --tmpdir /ramdisks \
        --threads 24 \
        --query $fastq \
        --db $dbfile \
        --outfmt 6 qtitle stitle bitscore evalue qframe qseq qstart qend sseq sstart send slen btop |
    convert-diamond-to-json.py | bzip2 > $out

    echo "  Sleeping 5 seconds to give the filesystem a chance to settle." >> $log
    sleep 5

    echo "  DIAMOND blastx stopped at `date`" >> $log
}


if [ $SP_SIMULATE = "0" ]
then
    echo "  This is not a simulation." >> $log
    if [ -f $out ]
    then
        if [ $SP_FORCE = "1" ]
        then
            echo "  Pre-existing output file $out exists, but --force was used. Overwriting." >> $log
            run_diamond
        else
            echo "  Will not overwrite pre-existing output file $out. Use --force to make me." >> $log
        fi
    else
        echo "  No pre-existing output file $out exists. Running DIAMOND." >> $log
        run_diamond
    fi
else
    echo "  This is a simulation." >> $log
fi

echo "04-diamond on task $task stopped at `date`" >> $log
echo >> $log
