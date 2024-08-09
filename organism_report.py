#!/usr/bin/env python3

# Import necessary libraries
import argparse
import os
import sys
import gzip
import subprocess
import pandas as pd
import csv
from collections import Counter
from colorama import init, Fore, Style
import plotly.express as px
import dash
from dash import html, dcc, Input, Output
from jinja2 import Environment, FileSystemLoader
from datetime import datetime
import pytaxonkit
from tabulate import tabulate
import shutil

# Argument parsing
def parse_args():
    parser = argparse.ArgumentParser(description=ascii_art + '''
    This script processes Centrifuge reports and FASTQ files for BLAST analysis and interpretation.
    For information contact Daniel Ward - GSTT
    ''', formatter_class=argparse.RawTextHelpFormatter)
    parser = argparse.ArgumentParser(description="For EPI2ME analysis use --epi2me_report. For CIDR metagenomics workflow use --centrifuge_report_dir.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-e', '--epi2me_report', required=False, help='Path to WIMP CSV file downloaded from EPI2ME')
    parser.add_argument('-c', '--centrifuge_report_dir', required=False, help='Path to the CIDR metagenomics Centrifuge report directory.')
    parser.add_argument('-o', '--organism', required=True, help='Organism name to search for in the Centrifuge report. The script extracts corresponding taxonomic IDs.')
    parser.add_argument('-f', '--fastq_dir', required=True, help='Directory containing .fastq.gz files. The script processes these files to extract relevant reads.')
    parser.add_argument('-d', '--output_dir', required=True, help='Output directory for storing results including extracted reads and BLAST output.')
    parser.add_argument('-b', '--blastdb', required=True, help='Prefix for the blastDB. Needs to be in the $BLASTDB path')
    
    return parser.parse_args()

# Parsing report files and extracting reads from FASTQ file
def ensure_trailing_slash(path):
    return os.path.join(path, '')

def extract_tax_ids(centrifuge_report, organism_name):
    print(ascii_art, flush=True)
    tax_ids = []
    with open(centrifuge_report, 'r') as file:
        for line in file:
            if organism_name in line:
                parts = line.strip().split('\t')
                tax_ids.append(parts[1])
    if not tax_ids:  # Check if tax_ids list is empty
        sys.stderr.write(f"{Fore.RED}Error: No tax IDs found for {organism_name}\n")
        sys.exit(1)  # Exit the script with an error status code, e.g., 1
    sys.stderr.write(f"{Fore.GREEN}Identified species: {organism_name}, Tax IDs: {', '.join(tax_ids)}\n")
    return tax_ids

def get_unique_children_taxids_and_names(taxids):
    """
    Calls pytaxonkit to get children taxonomic IDs and names for each taxonomic ID in the input list,
    including the query level itself. Ensures that each taxonomic ID and name is unique.
    
    Parameters:
    - taxids (List[int]): A list of taxonomic IDs.
    
    Returns:
    - List[int]: A list of unique children taxonomic IDs including the query level.
    - List[str]: A list of unique names corresponding to the taxonomic IDs.
    """
    taxid_name_map = {}
    for taxid in taxids:
        # Pytaxonkit list
        pytaxonkit_output = pytaxonkit.list([taxid], raw=True, data_dir="/taxonkit")

        def extract_ids_and_names(nested_dict):
            for key, value in nested_dict.items():
                # Extract taxonomic ID and name
                split_key = key.split(' ', 1)  # Splitting at the first space to separate taxid from the rest
                taxid = str(split_key[0])
                name = split_key[1] if len(split_key) > 1 else "No name"
                taxid_name_map[taxid] = name  # This ensures uniqueness of taxid
                # Recurse if there are more levels
                if isinstance(value, dict) and value:
                    extract_ids_and_names(value)

        # Starting point for extraction
        extract_ids_and_names(pytaxonkit_output)

    # Convert the dictionary back into two lists
    unique_taxids = list(taxid_name_map.keys())
    unique_names = list(taxid_name_map.values())

    return unique_taxids, unique_names

def extract_read_ids(raw_report, tax_ids, output_dir):
    matched_reads = []  # List to store tuples of (read_id, score)
    
    # Process file to collect all matched reads
    with open(raw_report, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if parts[2] in tax_ids:  # Assuming taxID is the third column
                read_id = parts[0]  # Assuming readID is the first column
                score = int(parts[3])  # Assuming score is an integer in the fourth column
                matched_reads.append((read_id, score))
    
    # Check if no matched reads were found
    if not matched_reads:
        sys.stderr.write(f"{Fore.RED}Error: No read IDs found matching the given Tax IDs\n")
        sys.exit(1)  # Exit the script with an error status code, e.g., 1
    
    # Sort matched reads by score in descending order and select the top N
    top_reads = sorted(matched_reads, key=lambda x: x[1], reverse=True)[:50]
    
    # Write the top N read_ids to the output file and count them
    read_count = 0
    read_ids = set()
    with open(os.path.join(output_dir, 'matched_read_ids.txt'), 'w') as id_file:
        for read_id, _ in top_reads:
            id_file.write(read_id + '\n')
            read_ids.add(read_id)
            read_count += 1
    
    # Output the number of reads matched with TaxIDs
    sys.stderr.write(f"{Fore.GREEN}Number of reads matched with TaxIDs (top 10 by score): {read_count}\n")
    
    return read_ids, matched_reads

def extract_read_ids_epi2me(epi2me_report, organism, output_dir):
    read_ids = set()
    read_count = 0

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(epi2me_report, 'r') as file, open(os.path.join(output_dir, 'matched_read_ids.txt'), 'w') as output_file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            if row['name'].lower() == organism.lower():
                read_id = row['readid']
                if read_id not in read_ids:
                    read_ids.add(read_id)
                    output_file.write(read_id + '\n')
                    read_count += 1

    sys.stderr.write(f"{Fore.GREEN}Number of reads matched with '{organism}': {read_count}\n")
    return read_ids

def extract_reads(fastq_dir, read_ids, output_dir):
    total_target_reads = 0
    total_files = 0
    total_reads = 0
    
    found_barcodes = set()

    with open(os.path.join(output_dir, 'concatenated_subset_reads.fasta'), 'w') as outfile:
        for fastq_file in os.listdir(fastq_dir):
            if fastq_file.endswith('.fastq.gz') or fastq_file.endswith('.fastq'):
                total_files += 1
                # Search for 'barcode' in the filename
                for part in fastq_file.split('_'):
                    if part.startswith('barcode'):
                        found_barcodes.add(part)

                if fastq_file.endswith('.fastq.gz'):
                    open_func = gzip.open
                    mode = 'rt'
                else:
                    open_func = open
                    mode = 'r'
                
                with open_func(os.path.join(fastq_dir, fastq_file), mode) as infile:
                    for line in infile:
                        total_reads += 1
                        if line.startswith('@'):
                            read_id = line.split()[0][1:]
                            if read_id in read_ids:
                                outfile.write(f">{read_id}\n")  # Write the header line in FASTA format
                                seq_line = next(infile)  # Read the next line which is the sequence
                                outfile.write(seq_line)  # Write the sequence line
                                total_target_reads += 1
                                # Skip the next two lines (quality header and quality score) in FASTQ
                                next(infile)
                                next(infile)

    # After processing, check if no target reads were found
    if total_target_reads == 0:
        sys.stderr.write(f"{Fore.RED}Error: No target reads found in the provided FASTQ files.\n")
        sys.exit(1)  # Exit the script with an error status code, e.g., 1

    # Print summary of found barcodes and total reads processed
    if found_barcodes:
        sys.stderr.write(f"Summary of barcodes found: {', '.join(sorted(found_barcodes))}\n")
    sys.stderr.write(f"Total FASTQ files processed: {total_files}\n")
    sys.stderr.write(f"{Fore.GREEN}Total reads extracted and converted to FASTA: {total_target_reads}\n")
    return total_reads

def run_blast(input_file, output_dir, blastdb):
    sys.stderr.write("BLAST analysis started...\n")
    cmd = [
        'blastn', '-evalue', '1e-5', '-max_target_seqs', '100', '-max_hsps', '1', '-outfmt', '11',
        '-query', os.path.join(output_dir, 'concatenated_subset_reads.fasta'),
        '-db', blastdb, '-out', os.path.join(output_dir, 'blast_results_11.tmp'), '-num_threads', '20'
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        # If BLAST encounters an error, print the error message and exit
        sys.stderr.write(f"{Fore.RED}Error running BLAST: {result.stderr}{Style.RESET_ALL}\n")
        sys.exit(1)
    
    sys.stderr.write(f"{Fore.GREEN}BLAST analysis completed.{Style.RESET_ALL}\n")

def run_parser(input_file, output_dir, blastdb):
    sys.stderr.write("Parsing BLAST report...\n")
    cmd = [
        'blast_formatter', '-archive', os.path.join(output_dir, 'blast_results_11.tmp'), '-num_alignments', '5', '-num_descriptions', '5', '-html', '-out', os.path.join(output_dir, 'blast_results.html')
    ]
    subprocess.run(cmd)
    
    cmd2 = [
        'blast_formatter', '-archive', os.path.join(output_dir, 'blast_results_11.tmp'), '-outfmt', '6 qseqid sseqid pident staxids sscinames pident  length mismatch gapopen qstart qend sstart send evalue bitscore qseq qlen', '-out', os.path.join(output_dir, 'blast_results_6.tmp')
    ]
    subprocess.run(cmd)
    subprocess.run(cmd2)
    sys.stderr.write("BLAST report parsing completed.\n")

def report_build(output_dir, organism, read_ids, blastdb, total_reads, fastq_dir, input_format, time_interval, matched_reads):
    file_path = os.path.join(output_dir, 'blast_results.html')
    start_string = '<b>Query=</b>'
    end_string = 'Effective search space used:'

    def extract_blast_reads(file_path, start_string, end_string):
            reads = {}
            capture = False
            read_id = 1
            read = ""

            with open(file_path, 'r') as file:
                for line in file:
                    if start_string in line:
                        capture = True
                        read = line
                        continue
                    elif end_string in line and capture:
                        reads[f"read_{read_id}"] = read
                        read_id += 1
                        capture = False
                    elif capture:
                        read += line

            return reads

    def parse_blast_output_with_stats(output_dir):
        file_path = os.path.join(output_dir, 'blast_results_6.tmp')
        
        # Check if the file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError("No BLAST output file found at the specified path.")

        try:
            # Read the TSV file
            df = pd.read_csv(file_path, sep='\t', header=None, names=[
                'qseqid', 'sseqid', 'pident', 'staxids', 'sscinames', 
                'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 
                'send', 'evalue', 'bitscore', 'qseq', 'qlen'
            ])
        except Exception as e:
            raise Exception(f"Error reading the BLAST output file: {e}")

        if df.empty:
            raise ValueError("The BLAST output file is empty or the format is incorrect.")
        
        # Group by qseqid and select the row with the lowest evalue for each group
        best_hits = df.loc[df.groupby('qseqid')['evalue'].idxmin()]

        # Count the frequency of each 'sscinames'
        sscinames_counts = Counter(best_hits['sscinames'])
        two_most_common = sscinames_counts.most_common(2)

        def get_stats(sscinames):
            hits = best_hits[best_hits['sscinames'] == sscinames]
            avg_qlen = round(hits['qlen'].mean())
            avg_pident = round(hits['pident'].mean())
            lowest_evalue = hits['evalue'].min()
            query_count = hits['qseqid'].nunique()  # Count of unique queries
            return sscinames, avg_qlen, avg_pident, lowest_evalue, query_count

        # Get stats for the most common 'sscinames'
        most_common_stats = get_stats(two_most_common[0][0]) if two_most_common else ('none', 'none', 'none', 'none', 'none')

        # Check if there is more than one organism
        if len(two_most_common) > 1:
            second_most_common_stats = get_stats(two_most_common[1][0])
        else:
            second_most_common_stats = ('none', 'none', 'none', 'none', 'none')

        return most_common_stats, second_most_common_stats

    most_common, second_most_common = parse_blast_output_with_stats(output_dir)

    reads = extract_blast_reads(file_path, start_string, end_string)
    # HTML for the accordion read view
    html_start1 = """
    <li class="page-item"><a class="page-link" href="#!" onclick="showSection({0})">{0}</a></li>
    """ 

    html_start2="""  
    <div id="section{0}" class="content-section">
    """

    html_end = """
        </div>
    """

    # Make a copy of the original dictionary
    original_reads = reads.copy()

    # First loop: Modify the original reads dictionary
    for count, i in enumerate(reads.keys()):
        reads[i] = html_start1.format(count)
    reads_list1 = '\n'.join(reads.values())

    # Second loop: Use the original_reads copy for reference
    for count, i in enumerate(original_reads.keys()):
        reads[i] = html_start2.format(count) + '<pre>' + original_reads[i] + '</pre>' + html_end
    reads_list2 = '\n'.join(reads.values())

    # Quick stats
    barcode = fastq_dir.split('/')[-1]
    sample_id = output_dir.split('/')[3]    

    report_dict = {"time": 'TEST' + " hrs",
                "title": "Clinical metagenomics report",
                "path": os.getcwd(),
                "SampleID": sample_id,
                "time": time_interval,
                "input_format": input_format,
                "Barcode": barcode,              
                "date": datetime.now(),
                "BLAST1": reads_list1,
                "BLAST2": reads_list2,
                "organism_read_count": len(matched_reads),
                "organism": organism,
                "blast_db": blastdb.split('/')[-1],
                "total_fastq_reads": total_reads,
                "most_common": most_common,
                "second_most_common": second_most_common,
                }

    # Prepare Jinja2 environment and template
    env = Environment(loader=FileSystemLoader('./template/'), 
                    autoescape = False)
    template = env.get_template('report_template.html')
    # Render the template with the extracted data
    rendered_html = template.render({"report": report_dict})
    # Output the rendered HTML to a file
    with open(os.path.join(output_dir, 'organism_report.html'), 'w', encoding='utf-8') as file:
        file.write(rendered_html)
    subprocess.Popen(["firefox" ,os.path.join(output_dir, 'organism_report.html')])

def dash_func(output_dir):
    app = dash.Dash(__name__)

    app.layout = html.Div([
        dcc.Graph(id='scatter-plot'),
        
        html.Label('Color Variable:'),
        dcc.Dropdown(id='color-dropdown', 
                    options=[
                        {'label': 'Scientific Name', 'value': 'Scientific Name'},
                        {'label': 'Alignment Length', 'value': 'Alignment Length'},
                        {'label': 'Query read ID', 'value': 'Query read ID'},
                    ],
                    value='Scientific Name')
    ])

    file_path = os.path.join(output_dir, 'blast_results_6.tmp')

    def process_blast_output(file_path):
        # Define column names as per BLAST output format 6
        col_names = ["Query read ID", "sseqid", "Percent Identity (%)", "staxids", "Scientific Name", "pident_2", 
                    "Alignment Length", "mismatch", "gapopen", "qstart", "qend", "sstart", 
                    "send", "evalue", "bitscore", "qseq", "qlen"]

        filtered_data = pd.read_csv(file_path, sep="\t", names=col_names)
        return filtered_data

    @app.callback(
        Output('scatter-plot', 'figure'),
        [Input('color-dropdown', 'value')]
    )
    
    def update_graph(color_variable):
        filtered_data = process_blast_output(file_path)
        fig = px.scatter(filtered_data, x='Percent Identity (%)', y='evalue', color=color_variable,
                        hover_data=['Scientific Name', 'Alignment Length'], log_y=True)
        fig.update_layout(xaxis_title='Percent Identity (%)', yaxis_title='E-value (log scale)')
        return fig
    app.run_server()

def cleanup(output_dir):
    path2 = os.path.join(output_dir, 'blast_results.html')
    path3 = os.path.join(output_dir, 'blast_results_11.tmp')
    os.remove(path2)
    os.remove(path3)
    path4 = os.path.join(output_dir, 'concatenated_subset_reads.fasta')
    subprocess.run(['gzip', path4])

def decompress_gzip(file_path):
    decompressed_file_path = file_path.rstrip('.gz')
    with gzip.open(file_path, 'rb') as f_in:
        with open(decompressed_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return decompressed_file_path

def compress_gzip(file_path):
    compressed_file_path = file_path + '.gz'
    with open(file_path, 'rb') as f_in:
        with gzip.open(compressed_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return compressed_file_path

ascii_art = r''' 
             ____ ____ _____ _____      ____ ___ ____  ____  
            / ___/ ___|_   _|_   _|    / ___|_ _|  _ \|  _ \ 
            | | _\___\  | |   | |_____| |    | || | | | |_) |
            | |_||___)| | |   | |_____| |___ | || |_| |  _ < 
            \____|____/ |_|   |_|      \____|___|____/|_| \_\
'''

# Initialise text colouring 
init()  # Initialize Colorama

def main():
    # Parsing arguments
    args = parse_args()
    # Subsetting for EPI2ME
    if args.epi2me_report:
        # Declaring input source
        input_format = "EPI2ME"
        time_interval = "EPI2ME"
        read_ids = extract_read_ids_epi2me(args.epi2me_report, args.organism, args.output_dir)
    # Subsetting reads CIDR workflow
    elif args.centrifuge_report_dir:
        # Declaring input source
        input_format = "CIDR metag"
        # Building CIDR pipeline centrifuge outputs
        raw_report_path = os.path.join(args.centrifuge_report_dir, 'centrifuge_raw.tsv')
        raw_report_gz_path = raw_report_path + '.gz'
        
        # Check if the gzipped file exists and decompress it
        if os.path.exists(raw_report_gz_path):
            raw_report = decompress_gzip(raw_report_gz_path)
            gzipped = True
        elif os.path.exists(raw_report_path):
            raw_report = raw_report_path
            gzipped = False
        else:
            print("Error: centrifuge_raw.tsv or centrifuge_raw.tsv.gz not found.", flush=True)
            sys.exit(1)

        time_interval = raw_report.split('/')[-3]
        centrifuge_report = os.path.join(args.centrifuge_report_dir, 'centrifuge_report.tsv')
        
        if args.organism.lower() == "unclassified":
            tax_ids = ['0']
            unique_taxids, unique_names = ['0'], ['Unclassified']
        else:
            tax_ids = extract_tax_ids(centrifuge_report, args.organism)
            unique_taxids, unique_names = get_unique_children_taxids_and_names(tax_ids)
        
        taxid_name_pairs = list(zip(unique_taxids, unique_names))
        print("Search widened for the following taxa:", flush=True)
        print(tabulate(taxid_name_pairs, headers=['Taxonomic ID', 'Name']), flush=True)
        read_ids, matched_reads = extract_read_ids(raw_report, unique_taxids, args.output_dir)
        
        # Recompress the file if it was decompressed
        if gzipped:
            compress_gzip(raw_report)
            os.remove(raw_report)
    else:
        print("Error: No workflow analysis outputs were given.", flush=True)
        sys.exit(1)

    total_reads = extract_reads(args.fastq_dir, read_ids, args.output_dir)    
    subset_reads = os.path.join(args.output_dir, 'concatenated_subset_reads.fastq')
    # Blast pipeline
    run_blast(subset_reads, args.output_dir, args.blastdb)
    run_parser(subset_reads, args.output_dir, args.blastdb)
    # Building report
    report_build(args.output_dir, args.organism, read_ids, args.blastdb, total_reads, args.fastq_dir, input_format, time_interval, matched_reads)
    cleanup(args.output_dir)
    dash_func(args.output_dir)

if __name__ == '__main__':
    main()