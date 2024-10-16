#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1

#find all files, return their md5sums to std out
#find . -xtype f -iname "*.fastq.gz" -exec zcat {} | md5sum
for f in $(find . -xtype f -name "*.fastq.gz" | sort -V);do echo $f;zcat $f | md5sum;done | paste - -
