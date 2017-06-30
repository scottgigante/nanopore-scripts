# nanopore-scripts
Various scripts and utilities for manipulation of nanopore data and general use of Torque (Milton).

## albacore_run.sh (and albacore.sh)
Runs albacore as an array job on Torque. Specifically for WEHI Milton, but should work on any Torque instance. 

Usage: `./albacore_run.sh [-h] [-i INPUT_DIR] [-o OUTPUT_DIR] [-c CONFIG] /path/to/fast5`
Default arguments: INPUT_DIR=raw, OUTPUT_DIR=albacore, CONFIG=r94_450bps_linear.cfg

Directory structure:
```
/path/to/fast5
/path/to/fast5/INPUT_DIR
/path/to/fast5/INPUT_DIR/0001_12345/[...].fast5
[...]
```

The basecalled reads will be output to `/path/to/fast5/OUTPUT_DIR`.

## bash.sh
Runs any bash command o Torque. Commands should be separated by escaped semicolons.

Usage: `qsub -F "echo 'hello'\; echo 'world'" bash.sh`

## check_alignment_identity.py
Outputs insertion, deletion, mismatch, total error and yield for any number of BAM files.

Usage: `./check_alignment_identity.py in1.bam in2.bam [...] out.tsv`

## gmap_run.sh (and gmap_check.sh, gmap_clean.sh, gmap.sh)
Runs GMAP as an array job on Torque. Specifically for WEHI Milton, but should work on any Torque instance.

Usage: `qsub -F "wehi_home/grpu_mritchie_0/AGRF_Data/2017_05_09/reads/fasta/albacore_no_phage.fasta" gmap_run.sh`
After completion of the above: `qsub -F "wehi_home/grpu_mritchie_0/AGRF_Data/2017_05_09/reads/fasta/albacore_no_phage.fasta" gmap_check.sh`
After completion with no new jobs launched: qsub -F "wehi_home/grpu_mritchie_0/AGRF_Data/2017_05_09/reads/fasta/albacore_no_phage.fasta" gmap_clean.sh`

TODO: incorporate this as a workflow with dependencies. Currently not viable due to Milton limitations.

## monitor_cpu.sh
Monitors ongoing CPU and memory usage via `ps` for jobs running on Milton.

Usage: `./monitor_cpu.sh "gmap.avx"`

## split_fasta.py
Splits a large FASTA file into many smaller files of equal size.

Usage: `./split_fasta.py /path/to/fasta.fa 2000`
`./split_fasta.py /path/to/reads.fastq 100 fastq`

## split_large_folder.sh
Splits a folder with far too many files inside into subfolders of a fixed size. Will also reduce size of existing subfolders, but will not enlarge existing subfolders at present.

Usage: `./split_large_folder.sh /path/to/folder 4000`
