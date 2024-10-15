version 1.0

struct FastqCollection {
    Array[Pair[File,Map[String,String]]] fastqs
}

struct Sample {
    Array[String]+ barcodes
    String name
}

workflow bclconvert {

  input {
    String runDirectory
    String mode
    String sampleName 
    Array[String]+ barcodes
    Array[Int] lanes = []
    String? basesMask
    Int mismatches = 1
  }

  parameter_meta {
    runDirectory: {
      description: "Illumina run directory (e.g. /path/to/191219_M00000_0001_000000000-ABCDE).",
      vidarr_type: "directory"
    }
    sampleName: "Sample (we accept only one)"
    barcodes: "Array of barcode(s)"
    mode: "Either dragen or hpc, selected mode will determined the backend used to run all this"
    lanes: "Extract reads only for specified lanes"
    basesMask: "An Illumina bases mask string to use. If absent, the one written by the instrument will be used."
    mismatches: "Number of mismatches to allow in the barcodes (usually, 1)"
  }


  meta {
    author: "Peter Ruzanov, Lawrence Heisler"
    email: "pruzanov@oicr.on.ca, lheisler@oicr.on.ca"
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
      fastq_read1: "FASTQ reads 1st in pair.",
      fastq_read2: "FASTQ reads 2nd in pair."
    }
  }

  call buildSamplesheet {
     input:
      sample = object{name: sampleName, barcodes: barcodes},
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

  call postprocessResults {
    input:
      runFolder = runDirectory,
      fastqList = select_first([runBclconvertHpc.fastqList, runBclconvertDragen.fastqList]),
      demultiplexStats = select_first([runBclconvertHpc.demultiplexStats, runBclconvertDragen.demultiplexStats]),
      fastqs = select_first([runBclconvertHpc.fastqs, runBclconvertDragen.fastqs])
  }

  output {
    Pair[File,Map[String,String]] fastq_read1 = postprocessResults.out.fastqs[0]
    Pair[File,Map[String,String]] fastq_read2 = postprocessResults.out.fastqs[1]
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
    String modules = "bclconvert-scripts/1.1"
    String samplesheetScript = "$BCLCONVERT_SCRIPTS_ROOT/bin/buildSamplesheet.py"
    Int memory = 4
    Int timeout = 2
  }

  parameter_meta {
    sample: "Object which holds a metadata for a sample"
    samplesheetScript: "Script for generating sample sheet"
    lanes: "Lanes to extract"
    basesMask: "An Illumina bases mask string to use. If absent, the one written by the instrument will be used"
    memory: "Memory allocated for running this task"
    modules: "Modules for running bclconvert task"
    timeout: "Timeout for building a samplesheet"
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
    String fastqCompression = 'gzip'
    Int fastqCompressionLevel = 1
    Int timeout = 40
    Int memory = 32
    String modules = "bclconvert/4.2.7-2"
    String? additionalParameters
  }

  parameter_meta {
    runFolder: "Run folder"
    runName: "Run name"
    sampleSheet: "File with sample sheet"
    firstTileOnly: "Flag for processing first tile only"
    noLaneSplitting: "Flag to disable lane splitting"
    fastqCompression: "Compression type of fastq files"
    fastqCompressionLevel: "Fastq compression level"
    timeout: "Timeout for this task"
    memory: "Memory allocated for running this task"
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
  --bcl-only-matched-reads true \
  --fastq-gzip-compression-level ~{fastqCompressionLevel} ~{additionalParameters}
  
  zip ~{runName}.reports.gz Reports/*
  >>>

  runtime {
      memory: "~{memory}G"
      modules: "~{modules}"
      timeout: "~{timeout}"
  }
  
  output { 
     File fastqList = "Reports/fastq_list.csv"
     File demultiplexStats = "Reports/Demultiplex_Stats.csv"
     Array[File]+ fastqs = glob("*fastq.gz")
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
    String fastqCompression = 'gzip'
    Int fastqCompressionLevel = 1
    Int timeout = 40
    String? additionalParameters
  }

  parameter_meta {
    runFolder: "Run folder"
    runName: "Run name"
    sampleSheet: "File with sample sheet"
    firstTileOnly: "Flag for processing first tile only"
    noLaneSplitting: "Flag to disable lane splitting"
    fastqCompression: "Compression type of fastq files"
    fastqCompressionLevel: "Fastq compression level"
    timeout: "Timeout for this task"
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
  --bcl-only-matched-reads true \
  --fastq-compression-format ~{fastqCompression} \
  --fastq-gzip-compression-level ~{fastqCompressionLevel} ~{additionalParameters}

  zip ~{runName}.reports.gz Reports/*
  >>>

  runtime {
      backend: "DRAGEN"
      timeout: "~{timeout}"
  }

  output {
     File fastqList = "Reports/fastq_list.csv"
     File demultiplexStats = "Reports/Demultiplex_Stats.csv"
     Array[File]+ fastqs = glob("*fastq.gz")
  }
}
  
# ========================================================
#   A task to postprocess (rename) files appropriately
# ========================================================
task postprocessResults {
  input {
    String runFolder
    String runName = basename(runFolder)
    File fastqList
    File demultiplexStats
    Array[File]+ fastqs
    String modules = "bclconvert-scripts/1.1"
    String bclconvertScript = "$BCLCONVERT_SCRIPTS_ROOT/bin/runBclconvert.py"
    Int timeout = 12
    Int memory = 8
  }

  parameter_meta {
    runFolder: "Run folder"
    runName: "Run name"
    fastqList: "File produced by bclconvert, list of all fasq files"
    demultiplexStats: "File produced by bclconvert, demultiplexing info"
    fastqs: "Fastq files produced by bclconvert"
    modules: "Module with python bclconvert scripts"
    bclconvertScript: "Script for generating sample sheet"
    timeout: "Timeout for this task"
    memory: "Memory allocated for running this task"
  }

  command <<<
  python3 ~{bclconvertScript} -r ~{runName} -d ~{demultiplexStats} -l ~{fastqList} -f ~{sep="," fastqs}
  >>>

  runtime {
      modules: "~{modules}"
      memory:  "~{memory}G"
      timeout: "~{timeout}"
  }

  output {
     FastqCollection out = read_json("outputs.json")
  }
}
  
