# Metagenomics organism query 
### Overview
This Python script automates querying a subset of NGS reads, based on previous taxonomic classification, using BLASTn. The output produces a report helping the reader understand BLAST results and find the most likely origin of the sequence.

### Dependancies
An explicit list to build a conda env can be found in the ```./conda/``` directory. In short:
* NCBI BLAST
### Inputs
The script parses the outputs of EPI2ME or the CIDR metagenomics workflows and raw sequencing reads. 

```
  -h, --help            show this help message and exit
  -e EPI2ME_REPORT, --epi2me_report EPI2ME_REPORT
                        Path to WIMP CSV file downloaded from EPI2ME
  -c CENTRIFUGE_REPORT, --centrifuge_report CENTRIFUGE_REPORT
                        Path to the human-readable Centrifuge report (TSV format) containing columns like Organism, Tax_ID, etc.
  -r RAW_REPORT, --raw_report RAW_REPORT
                        Path to the raw Centrifuge report (TSV format) containing columns like readID, seqID, taxID, etc.
  -o 'ORGANISM NAME', --organism 'ORGANISM NAME'
                        Organism name to search for in the Centrifuge report. The script extracts corresponding taxonomic IDs.
  -f FASTQ_DIR, --fastq_dir FASTQ_DIR
                        Directory containing .fastq.gz files. The script processes these files to extract relevant reads.
  -d OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory for storing results including extracted reads and BLAST output.
  -b BLASTDB, --blastdb BLASTDB
                        Prefix for the blastDB. Needs to be in the $BLASTDB path
```

Prerequisites
Python 3
Pandas
Plotly
Dash
Installation
To use this script, ensure you have Python 3 installed. Then, install the required Python packages using pip:

Copy code
pip install pandas plotly dash
Usage
To run the script, navigate to the directory containing the script and execute it with Python:

css
Copy code
python organism_report.py [arguments]
The script accepts various arguments for customization. Use the -h flag to view all available options.

Features
Data Processing: The script uses Pandas to process data related to organisms.
Visualization: It utilizes Plotly for generating interactive visualizations.
Web Interface: Dash is employed to create a web-based user interface for presenting the results.
Customization: Arguments can be passed to customize the data processing and visualization aspects.
Example
Here's an example command to run the script:

css
Copy code
python organism_report.py --input data.csv --output report.html
Replace data.csv with your input data file and report.html with the desired output file name.

Contributing
Contributions to this project are welcome. Please fork the repository, make your changes, and submit a pull request.
```