version 1.0

struct Sample {
    Array[String]+ barcodes
    String name
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
    String mode
    Sample sample 
    Array[Int] lanes = []
    String? basesMask
    Int mismatches = 1
  }

  parameter_meta {
    runDirectory: {
      description: "Illumina run directory (e.g. /path/to/191219_M00000_0001_000000000-ABCDE).",
      vidarr_type: "directory"
    }
    sample: "Sample (we accept only one) that will includes name and barcode(s)"
    mode: "Either dragen or hpc, selected mode will determined the backend used to run all this"
    lanes: "Extract reads only for specified lanes"
    basesMask: "An Illumina bases mask string to use. If absent, the one written by the instrument will be used."
    mismatches: "Number of mismatches to allow in the barcodes (usually, 1)"
  }


  meta {
    author: "Lawrence Heisler"
    email: "lheisler@oicr.on.ca"
    description: "Workflow to produce FASTQ files from an Illumina instrument's run directory using DRAGEN bclconvert"
    dependencies: [{
      name: "bclconvert",
      url: "https://emea.support.illumina.com/sequencing/sequencing_software/bcl-convert.html"
      },
      {
      name: "bclconvert-scripts",
      url: "https://gitlab.oicr.on.ca/ResearchIT/modulator/-/blob/master/code/gsi/70_bclconvert_scripts.yaml"
    }]
    output_meta: {
      fastqs: "A list of FASTQs generated and annotations that should be applied to them."
    }
  }

  call buildSamplesheet {
     input:
      sample = sample,
      lanes = lanes,
      basesMask = basesMask
  }

  if (mode == "hpc") {
   call runBclconvertHpc {
     input:
       runFolder = runDirectory,
       sampleSheet = buildSamplesheet.samplesheet
   }
  }

  if (mode == "dragen") {
   call runBclconvertDragen {
     input:
       runFolder = runDirectory,
       sampleSheet = buildSamplesheet.samplesheet
   }
  }

  FastqCollection finalResults = select_first([runBclconvertHpc.fastqs, runBclconvertDragen.fastqs])

  output {
    Array[FastqFile]+ fastqs = finalResults.fastqCollection
  }  

}

# ====================
# Build a sample sheet
# ====================
task buildSamplesheet{
  input {
    Sample sample
    Array[Int]? lanes
    String? basesMask
    String modules = "bclconvert-scripts/1.0"
    String samplesheetScript = "$BCLCONVERT_SCRIPTS_ROOT/bin/buildSamplesheet.py"
    Int memory = 4
    Int timeout = 2
  }

  parameter_meta {
    sample: "Object which holds a metadata for a sample"
    samplesheetScript: "Script for generating sample sheet"
    lanes: "Lanes to extract"
    basesMask: "An Illumina bases mask string to use. If absent, the one written by the instrument will be used."
    memory: "Memory allocated for running this task. Default 4"
    modules: "Modules for running bclconvert task"
    timeout: "Timeout for building a samplesheet. Default 2"
  }

  command <<<
    python3 ~{samplesheetScript} -i "~{write_json(sample)}" -m "~{basesMask}" -l "~{sep=',' lanes}"
  >>>

  runtime {
      memory: "~{memory}G"
      modules: "~{modules}"
      timeout: "~{timeout}"
  }

  output {
    File samplesheet = "samplesheet.csv"
  }
}

# =====================
# Run bclconvert on HPC
# =====================
task runBclconvertHpc {
  input {
    String runFolder
    String runName = basename(runFolder)
    File sampleSheet
    Boolean firstTileOnly = false
    Boolean noLaneSplitting = false
    Boolean onlyMatchedReads = true
    String fastqCompression = 'gzip'
    Int fastqCompressionLevel = 1
    Int timeout = 40
    Int memory = 32
    String modules = "bclconvert-scripts/1.0 bclconvert/4.2.7-2"
    String bclconvertScript = "$BCLCONVERT_SCRIPTS_ROOT/bin/runBclconvert.py"
    String? additionalParameters
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
    bclconvertScript: "Script for generating sample sheet"
    timeout: "Timeout for this task. Default 40"
    memory: "Memory allocated for running this task. Default 32"
    modules: "Modules for running bclconvert task"
    additionalParameters: "Pass parameters which were not exposed"
  }

  command <<<
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
  >>>

  runtime {
      memory: "~{memory}G"
      modules: "~{modules}"
      timeout: "~{timeout}"
  }
  
  output { 
     FastqCollection fastqs = read_json("outputs.json")
  }
}

# =================================
# Run bclconvert on DRAGEN hardware
# =================================
task runBclconvertDragen {
  input {
    String runFolder
    String runName = basename(runFolder)
    File sampleSheet
    Boolean firstTileOnly = false
    Boolean noLaneSplitting = false
    Boolean onlyMatchedReads = true
    String fastqCompression = 'gzip'
    String modules = "bclconvert-scripts/1.0"
    String bclconvertScript = "$BCLCONVERT_SCRIPTS_ROOT/bin/runBclconvert.py"
    Int fastqCompressionLevel = 1
    Int timeout = 40
    String? additionalParameters
  }

  parameter_meta {
    runFolder: "Run folder"
    runName: "Run name"
    sampleSheet: "File with sample sheet"
    bclconvertScript: "Script for generating sample sheet"
    firstTileOnly: "Flag for processing first tile only. Default false"
    noLaneSplitting: "Flag to disable lane splitting. Default false"
    onlyMatchedReads: "Process only matched reads. Default true"
    fastqCompression: "Compression type of fastq files. Default gzip"
    fastqCompressionLevel: "Fastq compression level. Default 1"
    modules: "Modules for running bclconvert task"
    timeout: "Timeout for this task. Default 40"
    additionalParameters: "Pass parameters which were not exposed"
  }

  command <<<
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
  >>>

  runtime {
      backend: "DRAGEN"
      timeout: "~{timeout}"
  }

  output {
     FastqCollection fastqs = read_json("outputs.json")
  }
}
  
  
  
  
