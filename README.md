# bclconvert

Workflow to produce FASTQ files from an Illumina instrument's run directory using DRAGEN bclconvert

## Overview

## Dependencies

* [bclconvert](https://emea.support.illumina.com/sequencing/sequencing_software/bcl-convert.html)


## Usage

### Cromwell
```
java -jar cromwell.jar run bclconvert.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`runDirectory`|String|The path to the instrument's output directory.
`runName`|String|The name of the run, this will be used for the output folder and as a file prefix
`samples`|Array[Sample]|array of Samples, that will includes names and barcodes
`modules`|String|Modules to run on hpc


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`basesMask`|String?|None|An Illumina bases mask string to use. If absent, the one written by the instrument will be used.
`mismatches`|Int|1|Number of mismatches to allow in the barcodes (usually, 1)
`timeout`|Int|40|The maximum number of hours this workflow can run for.


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`runBclconvert.firstTileOnly`|Boolean|false|Flag for processing first tile only. Default false
`runBclconvert.noLaneSplitting`|Boolean|false|Flag to disable lane splitting. Default false
`runBclconvert.onlyMatchedReads`|Boolean|true|Process only matched reads. Default true
`runBclconvert.fastqCompression`|String|'gzip'|Compression type of fastq files. Default gzip
`runBclconvert.fastqCompressionLevel`|Int|1|Fastq compression level. Default 1
`runBclconvert.memory`|Int|32|Memory allocated for running this task. Default 32


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
     python3 <<CODE
     import json
     with open("samplesheet.csv", "w") as ss:
       ss.write("[Data]\n")
       ss.write("Sample_ID,index,index2\n")
       with open("~{write_json(samples)}") as js:
         d=json.load(js)
         print(d)
         for sample in d['samples']:
            name=sample['name']
            barcodes=sample['barcodes']
            for barcode in barcodes:
              (bc1,bc2)=barcode.split("-")
              ss.write("%s,%s,%s\n" %(name,bc1,bc2))
     ss.close()
     CODE
 ```
 
 ### Run bclconvert
 
 Run Illumina's bclconvert to produce fastq files
 
 ```
   bcl-convert -f \
   --bcl-input-directory ~{runFolder} \
   --output-directory . \
   --sample-sheet ~{sampleSheet} \
   --no-lane-splitting ~{noLaneSplitting} \
   --first-tile-only ~{firstTileOnly} \
   --bcl-only-matched-reads ~{onlyMatchedReads} \
   --fastq-gzip-compression-level ~{fastqCompressionLevel}
   
   zip ~{runName}.reports.gz Reports/*
   
   python3 <<CODE
   # code to create an output json, cataloging all the fastq files
   import json
   import csv
   import os
   import glob
   import re  
   
   report_fastq="Reports/fastq_list.csv"
   report_demultiplex="Reports/Demultiplex_Stats.csv"
   samples={}
   with open(report_fastq,'r') as f:
       reader = csv.DictReader(f,delimiter=",")
       for row in reader:
           sid=row['RGSM']
           R1=row.get("Read1File","")
           R2=row.get("Read2File","")
           lane=row['Lane']
           if sid not in samples:
               samples[sid]={}
           if lane not in samples[sid]:
               samples[sid][lane]={"R1":"","R2":"","read_count":0}
 
           samples[sid][lane]["R1"]=os.path.basename(R1)
           samples[sid][lane]["R2"]=os.path.basename(R2)
         
 
   f.close()
   ufastqs=glob.glob('Undetermined*.fastq.gz')
   samples['Undetermined']={}
   for ufastq in ufastqs:
       fname=os.path.basename(ufastq)
       cap=re.search('_(R\d)_',fname)
       read=cap[1] if cap is not None else "R0"
       cap=re.search('_L00(\d)_',fname)
       lane=cap[1] if cap is not None else '1'
       if lane not in samples['Undetermined']:
           samples['Undetermined'][lane]={"R1":"","R2":"","read_count":0}
       samples['Undetermined'][lane][read]=fname
 
   with open(report_demultiplex) as f:
         reader = csv.DictReader(f,delimiter=",")
         for row in reader:
             sid=row["SampleID"]
             lane=row["Lane"]
             barcodes=row["Index"]
             read_count=row["# Reads"]
             read_count_mm0=row["# Perfect Index Reads"]
             read_count_mm1=row["# One Mismatch Index Reads"]
             read_count_mm2=row["# Two Mismatch Index Reads"]
           
             if sid not in samples:
                 samples[sid]={}
             if lane not in samples[sid]:
                 samples[sid][lane]={"R1":"","R2":"","read_count":0}
           
             samples[sid][lane]["read_count"]=read_count
   output={"fastqCollection":[]}
   for sid in samples:
     for lane in samples[sid]:
       read_count=samples[sid][lane]['read_count']
       if samples[sid][lane]['R1']:
         output["fastqCollection"].append({"name":sid,"fastqFile":{"left":samples[sid][lane]['R1'],"right":{"read_number":1,"read_count":read_count}}})
       if samples[sid][lane]['R2']:
         output["fastqCollection"].append({"name":sid,"fastqFile":{"left":samples[sid][lane]['R2'],"right":{"read_number":2,"read_count":read_count}}})
         
   with open('outputs.json',"w") as f:
       json.dump(output,f)
   f.close()
   CODE
 
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
