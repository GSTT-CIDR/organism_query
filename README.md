# CIDR Organism Query 
### Overview
This tool automates querying a subset of NGS reads based on Centrifuge taxanomic classification outputs using BLASTn. The output produces a report helping the reader understand BLAST results and find the most likely origin of the sequence within the confines of the database applied.

This tool has a streamlined and containerised implementation as part of the [CIDR metagenomics library](https://gstt-cidr.github.io/network_hub/).

**This tool is for reasearch use only - it has not been validated for diagnostic use**

### Database dependencies
An explicit list to build a conda env can be found in the ```./conda/``` directory. In short:
* NCBI BLAST DB
* NCBI taxa dump

### Inputs
The script parses the outputs of EPI2ME or the CIDR metagenomics workflows and raw sequencing reads. 

```
  -h, --help            show this help message and exit
  -e EPI2ME_REPORT, --epi2me_report EPI2ME_REPORT
                        Path to WIMP CSV file downloaded from EPI2ME
  -c CENTRIFUGE_REPORT_DIR, --centrifuge_report_dir CENTRIFUGE_REPORT_DIR
                        Path to the two Centrifuge outputs: report (TSV format) containing columns like Organism, Tax_ID, etc and the raw Centrifuge report (TSV format) containing columns 'readID, seqID, taxID' etc.
  -o 'ORGANISM NAME', --organism 'ORGANISM NAME'
                        Organism name to search for in the Centrifuge report. The script extracts corresponding taxonomic IDs.
  -f FASTQ_DIR, --fastq_dir FASTQ_DIR
                        Directory containing .fastq.gz files. The script processes these files to extract relevant reads.
  -d OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory for storing results including extracted reads and BLAST output.
  -b BLASTDB, --blastdb BLASTDB
                        Prefix for the blastDB. Needs to be in the $BLASTDB path
```

### Example usage 
Activate appropriate conda env:
```conda activate organism_query```
Run organism query:
```
python organism_report.py 
-c /media/grid/metagenomics/results/cidr_control/16_hours/centrifuge 
-o 'Human adenovirus 3' 
-f /data/GSTT_control_sample_01/GSTT_control_sample_01/20240424_1408_X4_FAY77387_d3868a4f/fastq_pass/barcode11/ 
-d ~/temp/organism_query/ 
-b /media/grid/metagenomics/db/blastdb/ref_viruses_rep_genomes
```
