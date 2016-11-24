## Sofia P pipeline spec

This repo contains a
[slurm-pipeline](https://github.com/acorg/slurm-pipeline) specification
file (`specification.json`) and associated scripts for processing
sequencing reads from Sofia P.

### Usage

```sh
$ git clone https://github.com/acorg/sofia-pipeline sample-dir
$ cd sample-dir
$ cp /your/data/*.fastq.gz .
$ make run
```

### Output

The scripts in `01-trim`, `02-map`, etc. are all submitted by `sbatch` for
execution under [SLURM](http://slurm.schedmd.com/). The final step,
`05-panel` leaves its output in `05-panel/out`.

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
