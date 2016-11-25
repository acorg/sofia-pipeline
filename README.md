# Sofia pipeline spec

This repo contains a
[slurm-pipeline](https://github.com/acorg/slurm-pipeline) specification
file (`specification.json`) and associated scripts for processing
sequencing reads from Sofia.

## Usage

So you should have a directory structure that looks like:

    drwxr-xr-x  2 tcj25 tcj25 4096 Nov 24 22:20 140715-12
    drwxr-xr-x  2 tcj25 tcj25 4096 Nov 24 22:20 140715-24
    drwxr-xr-x 10 tcj25 tcj25 4096 Nov 24 22:44 pipelines/initial/140715-12
    drwxr-xr-x 10 tcj25 tcj25 4096 Nov 24 22:23 pipelines/initial/140715-24

### Adding and processing new data

Working from the top level of this structure, to add and process a new
sample, `141110-79` you'd do this:

```sh
$ mkdir 141110-79
$ cp 1.fastq.gz 2.fastq.gz 141110-79  # Or whatever your two paired-end reads are.
$ mkdir -p pipelines/initial/141110-79
$ cd pipelines/initial/141110-79
$ git clone https://github.com/acorg/sofia-pipeline .
$ make run
```

### Adding a second pipeline

Assuming you already had another SLURM pipeline on GitHub at
`git clone https://github.com/acorg/second-stage-pipeline`,
you could add a second stage of processing via

```sh
$ mkdir pipelines/second-stage
```

And run the pre-existing `140715-12` data set through it via

```sh
$ mkdir pipelines/second-stage/140715-12
$ cd pipelines/second-stage/140715-12
$ git clone https://github.com/acorg/second-stage-pipeline .
$ make run
```

## Pipeline steps

    00-start
    01-trim
    02-flash
    03-map
    04-find-unmapped
    05-diamond
    06-panel
    07-stop

The scripts in `00-start` and `01-trim` assume they can find their
input FASTQ files in `../../../../SAMPLE-NAME` where `SAMPLE-NAME` is
something like `141110-79` (as in the above directory hierarchy).

## Output

The final step, `06-panel` leaves its output in `06-panel/out`.

### Logging

Logging goes to two places when you run `make run` in a directory such as
`pipelines/initial/140715-12`. 
* Top-level output for the initial and final scripts and for the `collect`
script in `06-panel` appears in `slurm-pipeline.log`.
* Per-FASTA output for the scripts that treat just one FASTA file
(`01-trim` to `05-diamond`) appears in a `.log` file whose name is
based on the FASTA files being processed.

### Cleaning up

There are 3 `Makefile` targets for cleaning up, `clean-1`, `clean-2`, and
`clean-3`. These increase in severity. See the `Makefile` to see what
exactly is removed by each target. You can run these by, e.g.,

```sh
$ make clean-1
```

`clean-1` throws away the large intermediate files created by the pipeline.
E.g., information about how many reads mapped to the human genome can be
found in the `03-find-unmapped/slurm-*.out`. Make sure you know what you're
doing if you clean up, else you may need to recreate these files at some
point.
