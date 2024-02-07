#!/usr/bin/env python3

#read parsing
import argparse
import os
import gzip
import subprocess
import sys
import pandas as pd
import csv
from collections import Counter
#interactivate plotting
import plotly.express as px
import dash
from dash import html, dcc, Input, Output
from jinja2 import Environment, FileSystemLoader
from datetime import datetime

#Argument parsing
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
    print(ascii_art)
    tax_ids = []
    with open(centrifuge_report, 'r') as file:
        for line in file:
            if organism_name in line:
                parts = line.strip().split('\t')
                tax_ids.append(parts[1])  # Assuming Tax_ID is the second column
    sys.stderr.write(f"Identified species: {organism_name}, Tax IDs: {', '.join(tax_ids)}\n")
    return tax_ids


def extract_read_ids(raw_report, tax_ids, output_dir):
    read_ids = set()
    read_count = 0
    with open(raw_report, 'r') as file, open(os.path.join(output_dir, 'matched_read_ids.txt'), 'w') as id_file:
        for line in file:
            parts = line.strip().split('\t')
            if parts[2] in tax_ids:  # Assuming taxID is the third column
                read_id = parts[0]  # Assuming readID is the first column
                read_ids.add(read_id)
                id_file.write(read_id + '\n')
                read_count += 1
    sys.stderr.write(f"Number of reads matched with TaxIDs: {read_count}\n")
#    print(read_ids)
    return read_ids



def extract_read_ids_epi2me(epi2me_report, organism, output_dir):
    read_ids = set()
    read_count = 0

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(epi2me_report, 'r') as file, open(os.path.join(output_dir, 'matched_read_ids.txt'), 'w') as output_file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            if row['name'] == organism:
                read_id = row['readid']
                if read_id not in read_ids:
                    read_ids.add(read_id)
                    output_file.write(read_id + '\n')
                    read_count += 1

    sys.stderr.write(f"Number of reads matched with '{organism}': {read_count}\n")
#    print(read_ids)
    return read_ids




def extract_reads(fastq_dir, read_ids, output_dir):
    total_target_reads = 0
    total_files = 0
    total_reads = 0
    
    found_barcodes = set()

    with open(os.path.join(output_dir, 'concatenated_subset_reads.fasta'), 'w') as outfile:
        for fastq_file in os.listdir(fastq_dir):
            if fastq_file.endswith('.fastq.gz'):
                total_files += 1
                # Search for 'barcode' in the filename
                for part in fastq_file.split('_'):
                    if part.startswith('barcode'):
                        found_barcodes.add(part)

                with gzip.open(os.path.join(fastq_dir, fastq_file), 'rt') as infile:
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

    # Print summary of found barcodes
    if found_barcodes:
        sys.stderr.write(f"Summary of barcodes found: {', '.join(sorted(found_barcodes))}\n")
    sys.stderr.write(f"Total FASTQ files processed: {total_files}\n")
    sys.stderr.write(f"Total reads extracted and converted to FASTA: {total_target_reads}\n")
    return total_reads

def run_blast(input_file, output_dir, blastdb):
    sys.stderr.write("BLAST analysis started...\n")
    cmd = [
        'blastn','-evalue', '1e-5','-max_target_seqs', '100', '-max_hsps', '1', '-outfmt', '11', '-query', output_dir + 'concatenated_subset_reads.fasta', '-db', blastdb, '-out', os.path.join(output_dir, 'blast_results_11.tmp'), '-num_threads', '20'
    ]
    subprocess.run(cmd)
    sys.stderr.write("BLAST analysis completed.\n")
    

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



##############################
##############################
##############################

def report_build(output_dir, organism, read_ids, blastdb, total_reads, fastq_dir, input_format, time_interval):
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
        # Read the TSV file without 'pident2'
        df = pd.read_csv(os.path.join(output_dir, 'blast_results_6.tmp'), sep='\t', header=None, names=[
            'qseqid', 'sseqid', 'pident', 'staxids', 'sscinames', 
            'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 
            'send', 'evalue', 'bitscore', 'qseq', 'qlen'
        ])

        # Group by qseqid and select the row with the lowest evalue for each group
        best_hits = df.loc[df.groupby('qseqid')['evalue'].idxmin()]

        # Count the frequency of each 'sscinames' in the best hits
        sscinames_counts = Counter(best_hits['sscinames'])

        # Find the two most common 'sscinames'
        two_most_common = sscinames_counts.most_common(2)

        # Function to get stats for a given 'sscinames'
        def get_stats(sscinames):
            hits = best_hits[best_hits['sscinames'] == sscinames]
            avg_qlen = round(hits['qlen'].mean())
            avg_pident = round(hits['pident'].mean())
            lowest_evalue = hits['evalue'].min()
            query_count = hits['qseqid'].nunique()  # Count of unique queries
            return sscinames, avg_qlen, avg_pident, lowest_evalue, query_count

        # Get stats for the most common and second most common 'sscinames'
        most_common_stats = get_stats(two_most_common[0][0])
        # Check if there is more than one organism
        if len(two_most_common) > 1:
            # Get stats for the second most common 'sscinames'
            second_most_common_stats = get_stats(two_most_common[1][0])
        else:
            # Return 'none' if there is only one organism
            second_most_common_stats = ('none', 'none', 'none', 'none', 'none')

        return most_common_stats, second_most_common_stats

    most_common, second_most_common = parse_blast_output_with_stats(output_dir)

    reads = extract_blast_reads(file_path, start_string, end_string)
    #HTML for the accordion read view
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

    #quick stats
    barcode = fastq_dir.split('/')[-1]
    sample_id = fastq_dir.split('/')[-4]    

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
                "organism_read_count": len(read_ids),
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
    subprocess.Popen(["firefox", os.path.join(output_dir, 'organism_report.html')])


        
    

##############################
##Dash and plotly function####
##############################

def dash_func(output_dir):
    app = dash.Dash(__name__)

    app.layout = html.Div([
        dcc.Graph(id='scatter-plot'),
        
#        html.Label('Length Threshold:'),
#        dcc.Slider(id='length-slider', min=100, max=500, step=10, value=200),
        
#        html.Label('Percent Identity Threshold:'),
#        dcc.Slider(id='pident-slider', min=50, max=100, step=5, value=80),

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

        # Read the BLAST output file
#        data = pd.read_csv(file_path, sep="\t", names=col_names)
        filtered_data = pd.read_csv(file_path, sep="\t", names=col_names)
        # Filter out rows where qlen is below the threshold and pident is below the pident threshold
#        filtered_data = data[(data["Alignment Length"] >= length_threshold) & (data["Percent Identity (%)"] >= pident_threshold)]

        return filtered_data


    @app.callback(
        Output('scatter-plot', 'figure'),
        [Input('color-dropdown', 'value')]
#        [Input('length-slider', 'value'),
#        Input('pident-slider', 'value'),
#        Input('color-dropdown', 'value')]
    )
    
    def update_graph(color_variable):
        filtered_data = process_blast_output(file_path)
        fig = px.scatter(filtered_data, x='Percent Identity (%)', y='evalue', color=color_variable,
                        hover_data=['Scientific Name', 'Alignment Length'], log_y=True)
        fig.update_layout(xaxis_title='Percent Identity (%)', yaxis_title='E-value (log scale)')
        return fig
    app.run_server()
        
        
def cleanup(output_dir):
    path1 = os.path.join(output_dir, blast_results_6.tmp)
    path2 = os.path.join(output_dir, blast_results.html)
    path2 = os.path.join(output_dir, blast_results_11.tmp)
    os.remove(path1)
    os.remove(path2)
    os.remove(path3)
    path4 = os.path.join(output_dir, concatenated_subset_reads.fasta)
    subprocess.run(['gzip', path4])
     
    
    
ascii_art = r''' 
             ____ ____ _____ _____      ____ ___ ____  ____  
            / ___/ ___|_   _|_   _|    / ___|_ _|  _ \|  _ \ 
            | | _\___\  | |   | |_____| |    | || | | | |_) |
            | |_||___)| | |   | |_____| |___ | || |_| |  _ < 
            \____|____/ |_|   |_|      \____|___|____/|_| \_\
'''


def main():
    # Parsing arguments
    args = parse_args()
    #Subsetting for EPI2ME
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
        raw_report = os.path.join(args.centrifuge_report_dir, 'centrifuge_raw.tsv')
        time_interval = raw_report.split('/')[-3]
        centrifuge_report = os.path.join(args.centrifuge_report_dir, 'centrifuge_report.tsv')
        tax_ids = extract_tax_ids(centrifuge_report, args.organism)
        read_ids = extract_read_ids(raw_report, tax_ids, args.output_dir)
    else:
        print("Error: No workflow analysis outputs were given.")
        sys.exit(1)
    total_reads = extract_reads(args.fastq_dir, read_ids, args.output_dir)    
    subset_reads = os.path.join(args.output_dir, 'concatenated_subset_reads.fastq')
    #Blast pipeline
    run_blast(subset_reads, args.output_dir, args.blastdb )
    run_parser(subset_reads, args.output_dir, args.blastdb )
    #Building report
    report_build(args.output_dir, args.organism, read_ids, args.blastdb, total_reads, args.fastq_dir, input_format, time_interval)
    dash_func(args.output_dir)
    cleanup(args.output_dir)
  
if __name__ == '__main__':
    main()

    
