# CliqueSNV-validation

The project provides scripts to validate CliqueSNV haplotyping tool and to compare it with other haplotyping tools such as aBayesQR, 2SNV, and PredictHaplo.
The details of validation can be found in the article https://www.biorxiv.org/content/10.1101/264242v1

## Project description and prerequisites
All the validation steps are stored in the `protocol.ipynb` file. It requires the haplotyping tools to be installed and sequencing reads to be downloaded.
The reads should be uploaded into `reads` folder at the root of the project. The project also requires alignment tool `bwa`.
Before running the tools the configuration files for aBayesQR and PredictHaplo may require to be modified.
The configuration files are located at `tool_configs`.
Begining of `protocol.ipynb` contains code chunks that set variables that are responsible for setting the paths. These chunks may also require modification.

## cite us
https://www.biorxiv.org/content/10.1101/264242v1
