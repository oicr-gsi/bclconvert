# bclconvert

Workflow to produce FASTQ files from an Illumina instrument's run directory using DRAGEN bclconvert

## Overview

## Dependencies

* [bclconvert](https://emea.support.illumina.com/sequencing/sequencing_software/bcl-convert.html)
* [bclconvert-scripts](https://gitlab.oicr.on.ca/ResearchIT/modulator/-/blob/master/code/gsi/70_bclconvert_scripts.yaml)


## Usage

### Cromwell
```
java -jar cromwell.jar run bclconvert.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`runDirectory`|String|{'description': 'Illumina run directory (e.g. /path/to/191219_M00000_0001_000000000-ABCDE).', 'vidarr_type': 'directory'}
`mode`|String|Either dragen or hpc, selected mode will determined the backend used to run all this
`sampleName`|String|Sample (we accept only one)
`barcodes`|Array[String]+|Array of barcode(s)


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`lanes`|Array[Int]|[]|Extract reads only for specified lanes
`basesMask`|String?|None|An Illumina bases mask string to use. If absent, the one written by the instrument will be used.
`mismatches`|Int|1|Number of mismatches to allow in the barcodes (usually, 1)


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`buildSamplesheet.modules`|String|"bclconvert-scripts/1.1"|Modules for running bclconvert task
`buildSamplesheet.samplesheetScript`|String|"$BCLCONVERT_SCRIPTS_ROOT/bin/buildSamplesheet.py"|Script for generating sample sheet
`buildSamplesheet.memory`|Int|4|Memory allocated for running this task
`buildSamplesheet.timeout`|Int|2|Timeout for building a samplesheet
`runBclconvertHpc.runName`|String|basename(runFolder)|Run name
`runBclconvertHpc.firstTileOnly`|Boolean|false|Flag for processing first tile only
`runBclconvertHpc.noLaneSplitting`|Boolean|false|Flag to disable lane splitting
`runBclconvertHpc.fastqCompression`|String|'gzip'|Compression type of fastq files
`runBclconvertHpc.fastqCompressionLevel`|Int|1|Fastq compression level
`runBclconvertHpc.timeout`|Int|40|Timeout for this task
`runBclconvertHpc.memory`|Int|32|Memory allocated for running this task
`runBclconvertHpc.modules`|String|"bclconvert/4.2.7-2"|Modules for running bclconvert task
`runBclconvertHpc.additionalParameters`|String?|None|Pass parameters which were not exposed
`runBclconvertDragen.runName`|String|basename(runFolder)|Run name
`runBclconvertDragen.firstTileOnly`|Boolean|false|Flag for processing first tile only
`runBclconvertDragen.noLaneSplitting`|Boolean|false|Flag to disable lane splitting
`runBclconvertDragen.fastqCompression`|String|'gzip'|Compression type of fastq files
`runBclconvertDragen.fastqCompressionLevel`|Int|1|Fastq compression level
`runBclconvertDragen.timeout`|Int|40|Timeout for this task
`runBclconvertDragen.additionalParameters`|String?|None|Pass parameters which were not exposed
`postprocessResults.runName`|String|basename(runFolder)|Run name
`postprocessResults.modules`|String|"bclconvert-scripts/1.1"|Module with python bclconvert scripts
`postprocessResults.bclconvertScript`|String|"$BCLCONVERT_SCRIPTS_ROOT/bin/runBclconvert.py"|Script for generating sample sheet
`postprocessResults.timeout`|Int|12|Timeout for this task
`postprocessResults.memory`|Int|8|Memory allocated for running this task


### Outputs

Output | Type | Description | Labels
---|---|---|---
`fastq_read1`|Pair[File,Map[String,String]]|FASTQ reads 1st in pair.|
`fastq_read2`|Pair[File,Map[String,String]]|FASTQ reads 2nd in pair.|


## Commands
This section lists commands run by bclconvert workflow
 
* Running bclconvert
 
### Create a samplesheet
 
This step creates a csv file with sample and barcode information. Required for bclconvert task
 
```
   python3 ~{samplesheetScript} -i "~{write_json(sample)}" -m "~{basesMask}" -l "~{sep=',' lanes}"
```
 
### Run bclconvert in HPC mode
 
Run Illumina's bclconvert to produce fastq files
 
```
   set -euo pipefail
   bcl-convert -f \
   --bcl-input-directory ~{runFolder} \
   --output-directory . \
   --sample-sheet ~{sampleSheet} \
   --no-lane-splitting ~{noLaneSplitting} \
   --first-tile-only ~{firstTileOnly} \
   --bcl-only-matched-reads true \
   --fastq-gzip-compression-level ~{fastqCompressionLevel} ~{additionalParameters}
   
   zip ~{runName}.reports.gz Reports/*
 ```
 ### Run bclconvert in DRAGEN mode
 
### Run bclconvert in DRAGEN mode
 
```
   set -euo pipefail
   dragen -f --bcl-conversion-only true \
   --bcl-input-directory ~{runFolder} \
   --output-directory . \
   --sample-sheet ~{sampleSheet} \
   --no-lane-splitting ~{noLaneSplitting} \
   --first-tile-only ~{firstTileOnly} \
   --bcl-only-matched-reads true \
   --fastq-compression-format ~{fastqCompression} \
   --fastq-gzip-compression-level ~{fastqCompressionLevel} ~{additionalParameters}
 
   zip ~{runName}.reports.gz Reports/*
 ```
 ### Post-process files
 
 Rename files using library id, sequencing run, lane and barcode. Assemble additional metadata into an object for the final output
 
 ```
   python3 ~{bclconvertScript} -r ~{runName} -d ~{demultiplexStats} -l ~{fastqList} -f ~{sep="," fastqs}
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
