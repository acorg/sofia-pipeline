#!/bin/bash -e

#SBATCH -J panel
#SBATCH -A DSMITH-BIOCLOUD
#SBATCH -o slurm-%A.out
#SBATCH -p biocloud-normal
#SBATCH --time=10:00:00

. /home/tcj25/.virtualenvs/35/bin/activate

# The log file is the overall top-level job log file, seeing as this step
# is a 'collect' step that is only run once.
log=../slurm-pipeline.log
outputDir=out
out=$outputDir/index.html

echo "06-panel started at `date`" >> $log

json=
fastq=
for task in "$@"
do
    echo "  task $task" >> $log
    json="$json ../05-diamond/$task.json.bz2"
    fastq="$fastq ../04-find-unmapped/$task-unmapped.fastq.gz"
done

dbfile=$HOME/scratch/root/share/ncbi/diamond-dbs/viral-protein-OKIAV-ECH.dmnd

if [ ! -f $dbfile ]
then
    echo "  DIAMOND database FASTA file $dbfile does not exist!" >> $log
    exit 1
fi


function panel()
{
    echo "  noninteractive-alignment-panel.py started at `date`" >> $log
    srun -n 1 noninteractive-alignment-panel.py \
      --json $json \
      --fastq $fastq \
      --matcher diamond \
      --outputDir $outputDir \
      --withScoreBetterThan 40 \
      --maxTitles 200 \
      --minMatchingReads 3 \
      --minCoverage 0.1 \
      --negativeTitleRegex phage \
      --diamondDatabaseFastaFilename $dbfile > summary-proteins
    echo "  noninteractive-alignment-panel.py stopped at `date`" >> $log

    echo "  proteins-to-viruses.py started at `date`" >> $log
    echo summary-proteins | proteins-to-viruses.py > summary-virus
    echo "  proteins-to-viruses.py stopped at `date`" >> $log
}

if [ $SP_SIMULATE = "0" ]
then
    echo "  This is not a simulation." >> $log
    if [ -f $out ]
    then
        if [ $SP_FORCE = "1" ]
        then
            echo "  Pre-existing output file $out exists, but --force was used. Overwriting." >> $log
            panel
        else
            echo "  Will not overwrite pre-existing output file $out. Use --force to make me." >> $log
        fi
    else
        echo "  No pre-existing output file $out exists. Making panel." >> $log
        panel
    fi
else
    echo "  This is a simulation." >> $log
fi

echo "06-panel stopped at `date`" >> $log
echo >> $log
