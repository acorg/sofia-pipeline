.PHONY: x, run, clean-1, clean-2, clean-3

x:
	@echo "There is no default make target. Use 'make run' to run the SLURM pipeline."

run:
	slurm-pipeline.py --sleep 2 --specification specification.json > status.json

# Remove all large intermediate files. Only run this if you're sure you
# want to throw away all that work!
clean-1:
	rm -f \
              01-trim/*.discarded.gz \
              01-trim/*.settings \
              01-trim/*.fastq.gz \
              01-trim/*.truncated.gz \
              02-flash/*.fastq.gz \
              03-map/*.bam \
              04-find-unmapped/*.fastq.gz \
              04-diamond/*.json.bz2

# Remove even more intermediate files.
clean-2: clean-1
	rm -f \
              01-trim/*.out \
              02-flash/out.* \
              02-flash/*.out \
              03-map/*.out \
              04-find-unmapped/*.out \
              04-diamond/*.out \
              05-panel/*.out \
              06-stop/*.out

# Remove all products, including the final panel output.
clean-3: clean-2
	rm -fr \
              05-panel/out \
              05-panel/summary-proteins \
              05-panel/summary-virus \
              slurm-pipeline.log \
              sbatch.log \
              slurm-pipeline.done \
              status.json
