import json
import csv
import os
import glob
import re

report_fastq = "Reports/fastq_list.csv"
report_demultiplex = "Reports/Demultiplex_Stats.csv"
samples = {}
with open(report_fastq, 'r') as f:
    reader = csv.DictReader(f, delimiter=",")
    for row in reader:
        sid = row['RGSM']
        R1 = row.get("Read1File", "")
        R2 = row.get("Read2File", "")
        lane = row['Lane']
        if sid not in samples:
            samples[sid] = {}
        if lane not in samples[sid]:
            samples[sid][lane] = {"R1": "", "R2": "", "read_count": 0}

        samples[sid][lane]["R1"] = os.path.basename(R1)
        samples[sid][lane]["R2"] = os.path.basename(R2)

f.close()
ufastqs = glob.glob('Undetermined*.fastq.gz')
samples['Undetermined'] = {}
for ufastq in ufastqs:
    fname = os.path.basename(ufastq)
    cap = re.search('_(R\d)_', fname)
    read = cap[1] if cap is not None else "R0"
    cap = re.search('_L00(\d)_', fname)
    lane = cap[1] if cap is not None else '1'
    if lane not in samples['Undetermined']:
        samples['Undetermined'][lane] = {"R1": "", "R2": "", "read_count": 0}
    samples['Undetermined'][lane][read] = fname

with open(report_demultiplex) as f:
    reader = csv.DictReader(f, delimiter=",")
    for row in reader:
        sid = row["SampleID"]
        lane = row["Lane"]
        barcodes = row["Index"]
        read_count = row["# Reads"]
        read_count_mm0 = row["# Perfect Index Reads"]
        read_count_mm1 = row["# One Mismatch Index Reads"]
        read_count_mm2 = row["# Two Mismatch Index Reads"]

        if sid not in samples:
            samples[sid] = {}
        if lane not in samples[sid]:
            samples[sid][lane] = {"R1": "", "R2": "", "read_count": 0}

        samples[sid][lane]["read_count"] = read_count
output = {"fastqCollection": []}
for sid in samples:
    for lane in samples[sid]:
        read_count = samples[sid][lane]['read_count']
        if samples[sid][lane]['R1']:
            output["fastqCollection"].append({"name": sid, "fastqFile": {"left": samples[sid][lane]['R1'],
                                                                         "right": {"read_number": 1,
                                                                                   "read_count": read_count}}})
        if samples[sid][lane]['R2']:
            output["fastqCollection"].append({"name": sid, "fastqFile": {"left": samples[sid][lane]['R2'],
                                                                         "right": {"read_number": 2,
                                                                                   "read_count": read_count}}})

with open('outputs.json', "w") as f:
    json.dump(output, f)
f.close()
