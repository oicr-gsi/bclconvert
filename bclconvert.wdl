version 1.0

struct Sample {
    Array[String] barcodes
    String name
}

struct SampleList {
    Array[Sample]+ samples
}

workflow bclconvert {
  input {
    String runDirectory
    String runName
    Array[Sample] samples
    String? basesMask
    Int mismatches = 1
    Int timeout = 40
  }
  parameter_meta {
    runDirectory: "The path to the instrument's output directory."
    runName: "The name of the run, this will be used for the output folder and as a file prefix"
    samples: "array of Samples, that will includes names and barcodes"
    basesMask: "An Illumina bases mask string to use. If absent, the one written by the instrument will be used."
    mismatches: "Number of mismatches to allow in the barcodes (usually, 1)"
    timeout: "The maximum number of hours this workflow can run for."
  }

  output {
    File samplesheet = buildsamplesheet.samplesheet
    File metrics = bclconvert.metrics
    File reports = bclconvert.reports
    Array[File] fastqFiles = bclconvert.fastqFiles
    Array[File] undeterminedFastqFiles =bclconvert.undeterminedFastqFiles
  }  

  meta {
    author: "Lawrence Heisler"
    description: "Workflow to produce FASTQ files from an Illumina instrument's run directory using DRAGEN bclconvert"
    dependencies: [{
      name: "bclconvert",
      url: "https://emea.support.illumina.com/sequencing/sequencing_software/bcl-convert.html"
    }]
    output_meta: {
      fastqs: "A list of FASTQs generated and annotations that should be applied to them."
    }
  }

  call buildsamplesheet {
     input:
      samples = object { samples: samples }
	  #samples = samples
  }

  call bclconvert {
    input:
      runFolder = runDirectory,
      runName = runName,
      sampleSheet = buildsamplesheet.samplesheet,
      timeout = timeout
  }
}

task buildsamplesheet{
  input {
	SampleList samples
  }	

  command <<<
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
  >>>

  output {
    File samplesheet = "samplesheet.csv"
  }
}


task bclconvert {
  input {
    String runFolder
    String runName
    File sampleSheet
    Boolean firstTileOnly = false
    Boolean noLaneSplitting = false
    Boolean onlyMatchedReads = false
    String fastqCompression = 'gzip'
    Int fastqCompressionLevel = 1
    Int timeout = 40
  }

  command <<<
  dragen -f --bcl-conversion-only true \
  --bcl-input-directory ~{runFolder} \
  --output-directory . \
  --sample-sheet ~{sampleSheet} \
  --no-lane-splitting ~{noLaneSplitting} \
  --first-tile-only ~{firstTileOnly} \
  --bcl-only-matched-reads ~{onlyMatchedReads} \
  --fastq-compression-format ~{fastqCompression} \
  --fastq-gzip-compression-level ~{fastqCompressionLevel}
  
  zip ~{runName}.reports.gz Reports/*
  
  ## https://bioinformatics.stackexchange.com/questions/19068/how-can-you-put-a-bash-array-into-a-wdl-variable
  #ls *fastq.gz | jq --raw-input . | jq --slurp > fastq_list.json
  
  
  ls *fastq.gz | grep -v Undetermined > fastq_list.txt
  ls Undetermined*fastq.gz > undetermined_fastq_list.txt 
  >>>

  runtime {
    backend: "DRAGEN"
    timeout: "~{timeout}"
  }
  
  output { 
     File metrics = "dragen.time_metrics.csv"
     File reports = "~{runName}.reports.gz"
     # Array[Fastq] fastqFiles = read_json("fastq_list.json")
     Array[File] fastqFiles = read_lines("fastq_list.txt")
     Array[File] undeterminedFastqFiles = read_lines("undetermined_fastq_list.txt")
   }
  

}

 
  
  
  
  
  