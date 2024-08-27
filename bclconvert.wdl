version 1.0

struct Sample {
    Array[String]+ barcodes
    String name
}

struct SampleList {
    Array[Sample]+ samples
}

struct FastqFile {
    String name
    Pair[File,Map[String,String]] fastqFile
}

struct FastqCollection {
    Array[FastqFile]+ fastqCollection
}




workflow bclconvert {

  input {
    String runDirectory
    String runName
    Array[Sample] samples 
    String? basesMask
    Int mismatches = 1
    String modules
    Int timeout = 40
  }

  parameter_meta {
    runDirectory: {
      description: "The path to the instrument's output directory.",
      vidarr_type: "directory"
    }
    runName: "The name of the run, this will be used for the output folder and as a file prefix"
    samples: "array of Samples, that will includes names and barcodes"
    basesMask: "An Illumina bases mask string to use. If absent, the one written by the instrument will be used."
    mismatches: "Number of mismatches to allow in the barcodes (usually, 1)"
    timeout: "The maximum number of hours this workflow can run for."
    modules: "Modules to run on hpc"
  }

  output {
    Array[FastqFile]+ fastqs = runBclconvert.fastqs.fastqCollection
  }  

  meta {
    author: "Lawrence Heisler"
    email: "lheisler@oicr.on.ca"
    description: "Workflow to produce FASTQ files from an Illumina instrument's run directory using DRAGEN bclconvert"
    dependencies: [{
      name: "bclconvert",
      url: "https://emea.support.illumina.com/sequencing/sequencing_software/bcl-convert.html"
    }]
    output_meta: {
      fastqs: "A list of FASTQs generated and annotations that should be applied to them."
    }
  }

  call buildSamplesheet {
     input:
      samples = object { samples: samples }
  }

  call runBclconvert {
    input:
      runFolder = runDirectory,
      runName = runName,
      sampleSheet = buildSamplesheet.samplesheet,
      timeout = timeout,
      modules = modules
  }
}

task buildSamplesheet{
  input {
    SampleList samples
  }

  parameter_meta {
    samples: "Object which holds an array of samples"
  }

  command <<<
    python3 <<CODE
    import json
    with open("samplesheet.csv", "w") as ss:
      ss.write("[Data]\n")
      ss_lines = []
      dualBarcodes = False
      with open("~{write_json(samples)}") as js:
        d = json.load(js)
        print(d)
        for sample in d['samples']:
          name = sample['name']
          barcodes = sample['barcodes']
          for barcode in barcodes:
            try:
              (bc1, bc2) = barcode.split("-")
              ss_lines.append(f'{name},{bc1},{bc2}\n')
              dualBarcodes = True
            except:
              ss_lines.append(f'{name},{barcode}\n')
      if dualBarcodes:
        ss.write("Sample_ID,index,index2\n")
      else:
        ss.write("Sample_ID,index\n")
      for line in ss_lines:
        ss.write(line)
      ss.close()
    CODE
  >>>

  output {
    File samplesheet = "samplesheet.csv"
  }
}


task runBclconvert {
  input {
    String runFolder
    String runName
    File sampleSheet
    Boolean firstTileOnly = false
    Boolean noLaneSplitting = false
    Boolean onlyMatchedReads = true
    String fastqCompression = 'gzip'
    Int fastqCompressionLevel = 1
    Int timeout = 40
    Int memory = 32
    String modules
  }

  parameter_meta {
    runFolder: "Run folder"
    runName: "Run name"
    sampleSheet: "File with sample sheet"
    firstTileOnly: "Flag for processing first tile only. Default false"
    noLaneSplitting: "Flag to disable lane splitting. Default false"
    onlyMatchedReads: "Process only matched reads. Default true"
    fastqCompression: "Compression type of fastq files. Default gzip"
    fastqCompressionLevel: "Fastq compression level. Default 1"
    timeout: "Timeout for this task. Default 40"
    memory: "Memory allocated for running this task. Default 32"
    modules: "Modules for running bclconvert task"
  }

  command <<<
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
  >>>

  runtime {
      memory: "~{memory}G"
      modules: "~{modules}"
      timeout: "~{timeout}"
  }
  
  output { 
     File reports = "~{runName}.reports.gz"
     FastqCollection fastqs = read_json("outputs.json")
  }
  

}

 
  
  
  
  
  
