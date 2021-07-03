# CliqueSNV-validation

The project provides scripts to validate CliqueSNV haplotyping tool and to compare it with other haplotyping tools such as aBayesQR, 2SNV, and PredictHaplo.
The details of validation can be found in the article https://pubmed.ncbi.nlm.nih.gov/34214168/

## Project description and prerequisites
All the validation steps are stored in the `protocol.ipynb` file. It requires the haplotyping tools to be installed and sequencing reads to be downloaded.
The reads should be uploaded into `reads` folder at the root of the project. The project also requires alignment tool `bwa`.
Before running the tools the configuration files for aBayesQR and PredictHaplo may require to be modified.
The configuration files are located at `tool_configs`.
Begining of `protocol.ipynb` contains code chunks that set variables that are responsible for setting the paths. These chunks may also require modification.

## cite us
Knyazev S, Tsyvina V, Shankar A, Melnyk A, Artyomenko A, Malygina T, Porozov YB, Campbell EM, Switzer WM, Skums P, Mangul S, Zelikovsky A. Accurate assembly of minority viral haplotypes from next-generation sequencing through efficient noise reduction. Nucleic Acids Res. 2021 Jul 2:gkab576. doi: 10.1093/nar/gkab576. Epub ahead of print. PMID: 34214168.

https://pubmed.ncbi.nlm.nih.gov/34214168/
