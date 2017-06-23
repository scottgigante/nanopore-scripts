#!/bin/bash

# Usage: ./split_large_folder.sh /path/to/folder 4000

DIR=$1
N=$2
BLOCK=200
M=$(($N / $BLOCK))
i=0
while [ -d "$DIR/$i" ]; do
  i=$(($i+1))
done

export DIR
export i
export M
export N
find $DIR -mindepth 1 -maxdepth 1 -type d | parallel -n 1 -X 'find {} -maxdepth 1 -mindepth 1 -type f | tail -n +$(($N+1))' | cat - <(find $DIR -mindepth 1 -maxdepth 1 -type f) | parallel -j 40 -n $BLOCK 'k=$(((({#}-1)/$M)+$i)); mkdir -p $DIR/$k; mv {} $DIR/$k'
