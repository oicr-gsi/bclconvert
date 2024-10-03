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
`buildSamplesheet.modules`|String|"bclconvert-scripts/1.0"|Modules for running bclconvert task
`buildSamplesheet.samplesheetScript`|String|"$BCLCONVERT_SCRIPTS_ROOT/bin/buildSamplesheet.py"|Script for generating sample sheet
`buildSamplesheet.memory`|Int|4|Memory allocated for running this task. Default 4
`buildSamplesheet.timeout`|Int|2|Timeout for building a samplesheet. Default 2
`runBclconvertHpc.runName`|String|basename(runFolder)|Run name
`runBclconvertHpc.firstTileOnly`|Boolean|false|Flag for processing first tile only. Default false
`runBclconvertHpc.noLaneSplitting`|Boolean|false|Flag to disable lane splitting. Default false
`runBclconvertHpc.onlyMatchedReads`|Boolean|true|Process only matched reads. Default true
`runBclconvertHpc.fastqCompression`|String|'gzip'|Compression type of fastq files. Default gzip
`runBclconvertHpc.fastqCompressionLevel`|Int|1|Fastq compression level. Default 1
`runBclconvertHpc.timeout`|Int|40|Timeout for this task. Default 40
`runBclconvertHpc.memory`|Int|32|Memory allocated for running this task. Default 32
`runBclconvertHpc.modules`|String|"bclconvert-scripts/1.0 bclconvert/4.2.7-2"|Modules for running bclconvert task
`runBclconvertHpc.bclconvertScript`|String|"$BCLCONVERT_SCRIPTS_ROOT/bin/runBclconvert.py"|Script for generating sample sheet
`runBclconvertHpc.additionalParameters`|String?|None|Pass parameters which were not exposed
`runBclconvertDragen.runName`|String|basename(runFolder)|Run name
`runBclconvertDragen.firstTileOnly`|Boolean|false|Flag for processing first tile only. Default false
`runBclconvertDragen.noLaneSplitting`|Boolean|false|Flag to disable lane splitting. Default false
`runBclconvertDragen.onlyMatchedReads`|Boolean|true|Process only matched reads. Default true
`runBclconvertDragen.fastqCompression`|String|'gzip'|Compression type of fastq files. Default gzip
`runBclconvertDragen.modules`|String|"bclconvert-scripts/1.0"|Modules for running bclconvert task
`runBclconvertDragen.bclconvertScript`|String|"$BCLCONVERT_SCRIPTS_ROOT/bin/runBclconvert.py"|Script for generating sample sheet
`runBclconvertDragen.fastqCompressionLevel`|Int|1|Fastq compression level. Default 1
`runBclconvertDragen.timeout`|Int|40|Timeout for this task. Default 40
`runBclconvertDragen.additionalParameters`|String?|None|Pass parameters which were not exposed


### Outputs

Output | Type | Description | Labels
---|---|---|---
`fastqs`|Array[FastqFile]+|A list of FASTQs generated and annotations that should be applied to them.|


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
   --bcl-only-matched-reads ~{onlyMatchedReads} \
   --fastq-gzip-compression-level ~{fastqCompressionLevel} ~{additionalParameters}
   
   zip ~{runName}.reports.gz Reports/*
  
   python3 ~{bclconvertScript}
```
 
### Run bclconvert in DRAGEN mode
 
```
   set -euo pipefail
   dragen -f --bcl-conversion-only true \
   --bcl-input-directory ~{runFolder} \
   --output-directory . \
   --sample-sheet ~{sampleSheet} \
   --no-lane-splitting ~{noLaneSplitting} \
   --first-tile-only ~{firstTileOnly} \
   --bcl-only-matched-reads ~{onlyMatchedReads} \
   --fastq-compression-format ~{fastqCompression} \
   --fastq-gzip-compression-level ~{fastqCompressionLevel} ~{additionalParameters}
 
   zip ~{runName}.reports.gz Reports/*
 
   python3 ~{bclconvertScript}
```
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
