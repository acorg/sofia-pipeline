#!/usr/bin/env python

"""
Read a SAM/BAM file (given as the first argument on the command line
and print FASTQ of the unmapped reads.
"""

from __future__ import division, print_function

import pysam
import sys

filename = sys.argv[1]
samfile = pysam.AlignmentFile(filename, 'rb')

# print('According to pysam file %r has %d mapped reads and %d unmapped.' %
# (filename, samfile.mapped, samfile.unmapped), file=sys.stderr)

mappedCount = unmappedCount = 0

for read in samfile:
    if read.is_unmapped:
        if read.is_paired:
            readId = read.query_name + ('/1' if read.is_read1 else '/2')
        else:
            readId = read.query_name
        unmappedCount += 1
        quality = bytes(map(lambda x: x + 33,
                            read.query_qualities)).decode('utf-8')
        # Note that printing the query name a second time (after the '+')
        # is optional in FASTQ. Unfortunately, some tools require it.
        print('@%s\n%s\n+%s\n%s' % (readId, read.query_sequence,
                                    readId, quality))
    else:
        mappedCount += 1

# Check that what we counted is what the SAM file says it contains.
# assert samfile.mapped == mappedCount, (
#     'samfile.mapped (%d) != mappedCount (%d)' %
#     (samfile.mapped == mappedCount))

# assert samfile.unmapped == unmappedCount, (
#     'samfile.unmapped (%d) != unmappedCount (%d)' %
#     (samfile.unmapped == unmappedCount))

samfile.close()

total = mappedCount + unmappedCount

print('Read %d reads in total. %d (%.2f%%) mapped, %d unmapped.' % (
    total, mappedCount,
    (mappedCount / total * 100.0) if total else 0.0,
    unmappedCount), file=sys.stderr)
