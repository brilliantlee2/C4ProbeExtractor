# C4ProbeExtractor
Extracting reads containing the C4 probe sequence using BWA

```
usage: GetProbedReads.py [-h] [--fastq FASTQ] [--probe PROBE] [--output_fastq OUTPUT_FASTQ] [-t THREADS]

options:
  -h, --help            show this help message and exit
  --fastq FASTQ         FASTQ of Cyclone reads
  --probe PROBE         FASTA of Cyclone probe sequence(s)
  --output_fastq OUTPUT_FASTQ
                        Output file name for stranded FASTQ entries
  -t THREADS, --threads THREADS
                        Threads to use [8]

```