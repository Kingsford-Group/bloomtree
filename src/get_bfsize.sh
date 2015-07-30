#!/bin/bash

# Directory storing fasta.gz files to be counted
fastaList="$1"
kmer=20
jfsize="200M"
threads=10

tempfile="mer_counts.jf"

#Set up jellyfish generators file
rm generators

while read -r file ; do
	echo gunzip -c $file >> generators
done < $fastaList

jellyfish count -m $kmer -s $jfsize -t $threads -C -o $tempfile -g generators -G 1

#Output counts for 1, 2, 3+ occurence kmers
echo "Unique canonical kmer counts by number of occurences:" 
jellyfish histo $tempfile -h 2 
