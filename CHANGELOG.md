# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2024-10-12
### Changed
- Regression testing for multiple scenarios, validates in Jenkins
- runDir type is changed to the same type bcl2fastq uses
- README re-generated using gsi wdl tools
- default parameters adjusted according to [GBS-5046](https://jira.oicr.on.ca/browse/GBS-5046)
- workflows bclconvert and dragen-bclconvert fused into one
- output type changed to Pair[File,Map{String,String]], one for each of the fastqs
### Added
- added commands.txt file
- support lanes and basesMasks parameters
- added two modes - hps and dragen
- simplified inputs/outputs
- post-processing step

## [Unreleased] - 2024-05-29
Initial Implementation
- [GBS-5046](https://jira.oicr.on.ca/browse/GBS-5046)


