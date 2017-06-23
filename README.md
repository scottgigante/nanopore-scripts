# nanopore-scripts
Various scripts and utilities for manipulation of nanopore data

# albacore_run.sh (and albacore.sh)
Runs albacore as an array job on Torque. Specifically for WEHI Milton, but should work on any Torque instance. 

Usage: ./albacore_run.sh [-h] [-i INPUT_DIR] [-o OUTPUT_DIR] [-c CONFIG] /path/to/fast5

Directory structure:
/path/to/fast5
/path/to/fast5/INPUT_DIR
/path/to/fast5/INPUT_DIR/0001_12345/[...].fast5
[...]
