#!/bin/bash -e

#SBATCH -J flash
#SBATCH -A DSMITH-BIOCLOUD
#SBATCH -o slurm-%A.out
#SBATCH -p biocloud-normal
#SBATCH --time=05:00:00

. /home/tcj25/.virtualenvs/35/bin/activate

task=$1
log=../slurm-pipeline.log
fastq1=../01-trim/${task}_R1_001-trimmed.fastq.gz
fastq2=../01-trim/${task}_R2_001-trimmed.fastq.gz
out1=${task}_R1_001-flash.fastq.gz
out2=${task}_R2_001-flash.fastq.gz
outMerged=${task}-merged-flash.fastq.gz

function run_flash()
{
    # FLASH will print a warning such as:
    #
    #   [FLASH] WARNING: An unexpectedly high proportion of combined pairs
    #   (30.25%) overlapped by more than 65 bp, the --max-overlap (-M)
    #   parameter.  Consider increasing this parameter.  (As-is, FLASH is
    #   penalizing overlaps longer than 65 bp when considering them for
    #   possible combining!)
    #
    # It's not clear what to do about this (if anything) because our reads
    # are not of a (nearly) constant length. So passing a fixed value to
    # --max-overlap doesn't help. E.g., we can pass 99 and still get the
    # error.

    echo "  FLASH started at `date`" >> $log
    srun -n 1 flash --compress --threads 1 $fastq1 $fastq2
    mv out.notCombined_1.fastq.gz $out1
    mv out.notCombined_2.fastq.gz $out2
    mv out.extendedFrags.fastq.gz $outMerged
    echo "  Sleeping 5 seconds to give the filesystem a chance to settle." >> $log
    sleep 5
    echo "  FLASH stopped at `date`" >> $log
}


echo "02-flash on task $task started at `date`" >> $log
echo "  fastq1 is $fastq1" >> $log
echo "  fastq2 is $fastq2" >> $log

if [ $SP_SIMULATE = "0" ]
then
    echo "  This is not a simulation." >> $log
    if [ -f $out1 -a -f $out2 ]
    then
        if [ $SP_FORCE = "1" ]
        then
            echo "  Pre-existing output files $out1 and $out2 exist, but --force was used. Overwriting." >> $log
            run_flash
        else
            echo "  Will not overwrite pre-existing output files $out1 and $out2. Use --force to make me." >> $log
        fi
    else
        echo "  Pre-existing output files $out1 and $out2 do not both exist. Running flash." >> $log
        run_flash
    fi
else
    echo "  This is a simulation." >> $log
fi

echo "02-flash on task $task stopped at `date`" >> $log
echo >> $log
