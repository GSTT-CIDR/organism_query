#!/usr/bin/env python3

import argparse
import os
import sys
import gzip
import subprocess
import pandas as pd
import csv
import random
from collections import Counter
from colorama import init, Fore, Style
import plotly.express as px
from jinja2 import Environment, FileSystemLoader
from datetime import datetime
import pytaxonkit
from tabulate import tabulate
import shutil

# Could do with refactoring into classes and modules. Specifically to unify Auto Query and this script.

# Set environment variables for BLAST
os.environ['NCBI_CONFIG_OVERRIDES'] = "TRUE"
os.environ['BLASTDB'] = "/mnt/db/blastdb"

def parse_args():
    parser = argparse.ArgumentParser(description=ascii_art + '''\n    This script processes Centrifuge reports and FASTQ files for BLAST analysis and interpretation.\n    For information contact Daniel Ward - GSTT\n    ''', formatter_class=argparse.RawTextHelpFormatter)
    parser = argparse.ArgumentParser(description="For EPI2ME analysis use --epi2me_report. For CIDR metagenomics workflow use --centrifuge_report_dir.", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-e', '--epi2me_report', required=False,
                        help='Path to WIMP CSV file downloaded from EPI2ME')
    parser.add_argument('-c', '--centrifuge_report_dir', required=False,
                        help='Path to the CIDR metagenomics Centrifuge report directory.')
    parser.add_argument('-o', '--organism', required=True,
                        help='Organism name to search for in the Centrifuge report. The script extracts corresponding taxonomic IDs.')
    parser.add_argument('-f', '--fastq_dir', required=True,
                        help='Directory containing .fastq.gz files. The script processes these files to extract relevant reads.')
    parser.add_argument('-d', '--output_dir', required=True,
                        help='Output directory for storing results including extracted reads and BLAST output.')
    parser.add_argument('-b', '--blastdb', required=True,
                        help='Prefix for the blastDB. Needs to be in the $BLASTDB path')
    parser.add_argument('-n', '--top_n', type=int, default=50,
                        help='Number of reads to subset in extract_read_ids.')
    parser.add_argument('-r', '--random_sort', action='store_true',
                        help='If set, select random reads instead of top scoring reads.')
    parser.add_argument('-t','--threads', type=int, default=15,
                        help='Number of threads to use for BLAST.')
    parser.add_argument('-q','--quiet', action='store_true',
                        help='Suppress all stdout and stderr output.')

    return parser.parse_args()

def ensure_trailing_slash(path):
    return os.path.join(path, '')

def extract_tax_ids(centrifuge_report, organism_name):
    """ Extract taxonomic IDs for a given organism from a Centrifuge report."""
    print(ascii_art, flush=True)
    print("CIDR Organism Query")
    print("NHS Clinical Metagenomics Platform")
    print("Written by Daniel Ward")
    tax_ids = []
    with open(centrifuge_report, 'r') as file:
        for line in file:
            # Convert both to lowercase for a case-insensitive check
            if organism_name.lower() in line.lower():
                parts = line.strip().split('\t')
                tax_ids.append(parts[1])
    if not tax_ids:
        sys.stderr.write(f"{Fore.RED}Error: No tax IDs found for {organism_name}\n")
        sys.exit(1)
    sys.stderr.write(f"{Fore.GREEN}Identified species: {organism_name}, Tax IDs: {', '.join(tax_ids)}\n")
    return tax_ids


def get_unique_children_taxids_and_names(taxids):
    """Resolves taxonomy and finds all unique child taxids and names for given taxids."""
    taxid_name_map = {}
    for taxid in taxids:
        pytaxonkit_output = pytaxonkit.list([taxid], raw=True, data_dir="/mnt/db/ref/refseq/taxonomy")

        def extract_ids_and_names(nested_dict):
            for key, value in nested_dict.items():
                split_key = key.split(' ', 1)
                tid = str(split_key[0])
                name = split_key[1] if len(split_key) > 1 else "No name"
                taxid_name_map[tid] = name
                if isinstance(value, dict) and value:
                    extract_ids_and_names(value)

        extract_ids_and_names(pytaxonkit_output)

    unique_taxids = list(taxid_name_map.keys())
    unique_names = list(taxid_name_map.values())
    return unique_taxids, unique_names

def extract_read_ids(raw_report, tax_ids, output_dir, top_n=50, random_sort=False):
    # Dictionary to store best score for each read ID
    read_scores = {}
    
    with open(raw_report, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if parts[2] in tax_ids:
                read_id = parts[0]
                score = int(parts[3])
                
                # Keep only the highest score for each read ID
                if read_id not in read_scores or score > read_scores[read_id]:
                    read_scores[read_id] = score
    
    # Convert to list of tuples (read_id, score)
    unique_matched_reads = [(read_id, score) for read_id, score in read_scores.items()]
    
    if not unique_matched_reads:
        sys.stderr.write(f"{Fore.RED}Error: No read IDs found matching the given Tax IDs\n")
        sys.exit(1)
    
    # Report on total unique reads before subsetting
    total_unique_reads = len(unique_matched_reads)
    
    if random_sort:
        # Random selection of unique reads
        subset_size = min(top_n, total_unique_reads)
        selected_reads = random.sample(unique_matched_reads, subset_size)
    else:
        # Top N selection of unique reads by score
        selected_reads = sorted(unique_matched_reads, key=lambda x: x[1], reverse=True)[:top_n]
    
    read_count = 0
    read_ids = set()
    with open(os.path.join(output_dir, 'matched_read_ids.txt'), 'w') as id_file:
        for read_id, score in selected_reads:
            id_file.write(read_id + '\n')
            read_ids.add(read_id)
            read_count += 1
    
    sys.stderr.write(f"{Fore.GREEN}Found {total_unique_reads} unique reads matching TaxIDs\n")
    sys.stderr.write(f"{Fore.GREEN}Number of reads selected for analysis (subset size={read_count}): {read_count}\n")
    return read_ids, unique_matched_reads

def extract_read_ids_epi2me(epi2me_report, organism, output_dir):
    read_ids = set()
    read_count = 0

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(epi2me_report, 'r') as file, open(os.path.join(output_dir, 'matched_read_ids.txt'), 'w') as output_file:
        csv_reader = csv.DictReader(file)
        for row in csv_reader:
            if row['name'].lower() == organism.lower():
                rid = row['readid']
                if rid not in read_ids:
                    read_ids.add(rid)
                    output_file.write(rid + '\n')
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
                                outfile.write(f">{read_id}\n")
                                seq_line = next(infile)
                                outfile.write(seq_line)
                                total_target_reads += 1
                                next(infile)
                                next(infile)

    if total_target_reads == 0:
        sys.stderr.write(f"{Fore.RED}Error: No target reads found in the provided FASTQ files.\n")
        sys.exit(1)

    if found_barcodes:
        sys.stderr.write(f"Summary of barcodes found: {', '.join(sorted(found_barcodes))}\n")
    sys.stderr.write(f"Total FASTQ files processed: {total_files}\n")
    sys.stderr.write(f"{Fore.GREEN}Total reads extracted and converted to FASTA: {total_target_reads}\n")
    return total_reads

def run_blast(input_file, output_dir, blastdb, threads=20):
    sys.stderr.write("BLAST analysis started...\n")
    cmd = [
        'blastn', '-evalue', '1e-5', '-max_target_seqs', '100', '-max_hsps', '1', '-outfmt', '11',
        '-query', os.path.join(output_dir, 'concatenated_subset_reads.fasta'),
        '-db', blastdb, '-out', os.path.join(output_dir, 'blast_results_11.tmp'),
        '-num_threads', str(threads)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.stderr.write(f"{Fore.RED}Error running BLAST: {result.stderr}{Style.RESET_ALL}\n")
        sys.exit(1)
    sys.stderr.write(f"{Fore.GREEN}BLAST analysis completed.{Style.RESET_ALL}\n")

def run_parser(input_file, output_dir, blastdb):
    sys.stderr.write("Parsing BLAST report...\n")
    cmd = [
        'blast_formatter', '-archive', os.path.join(output_dir, 'blast_results_11.tmp'),
        '-num_alignments', '5', '-num_descriptions', '5', '-html',
        '-out', os.path.join(output_dir, 'blast_results.html')
    ]
    subprocess.run(cmd)

    cmd2 = [
        'blast_formatter', '-archive', os.path.join(output_dir, 'blast_results_11.tmp'),
        '-outfmt', '6 qseqid sseqid pident staxids sscinames length mismatch gapopen qstart qend sstart send evalue bitscore qseq qlen',
        '-out', os.path.join(output_dir, 'blast_results_6.tmp')
    ]
    subprocess.run(cmd)
    subprocess.run(cmd2)
    sys.stderr.write("BLAST report parsing completed.\n")

def report_build(output_dir, organism, read_ids, blastdb, total_reads, fastq_dir,
                 input_format, time_interval, matched_reads):
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
        file_path_6 = os.path.join(output_dir, 'blast_results_6.tmp')
        if not os.path.exists(file_path_6):
            raise FileNotFoundError("No BLAST output file found at the specified path.")
        try:
            df = pd.read_csv(file_path_6, sep='\t', header=None, names=[
                'qseqid', 'sseqid', 'pident', 'staxids', 'sscinames',
                'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
                'send', 'evalue', 'bitscore', 'qseq', 'qlen'
            ])
        except Exception as e:
            raise Exception(f"Error reading the BLAST output file: {e}")
        if df.empty:
            raise ValueError("The BLAST output file is empty or the format is incorrect.")

        # Get best hits
        # There was a bug here where best hits were erroneously sorted on evalue - usually evalue = 0.0. BLAST outfmt 6 ordered by bitscore so we got away with it. 
        best_hits = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]

        # Replace NaN sscinames with a placeholder
        best_hits['sscinames'] = best_hits['sscinames'].fillna('unknown')

        sscinames_counts = Counter(best_hits['sscinames'])
        two_most_common = sscinames_counts.most_common(2)

        def get_stats(species_name):
            # If the placeholder is 'unknown', return a safe tuple
            if species_name == 'unknown':
                return (species_name, 'N/A', 'N/A', 'N/A', 0)

            # Filter rows for this species name
            hits = best_hits[best_hits['sscinames'] == species_name]
            if hits.empty:
                # If none match, safely return placeholders
                return (species_name, 'N/A', 'N/A', 'N/A', 0)

            # Otherwise, compute your stats
            avg_qlen = round(hits['qlen'].mean())
            avg_pident = round(hits['pident'].mean())
            lowest_evalue = hits['evalue'].min()
            query_count = hits['qseqid'].nunique()
            return (species_name, avg_qlen, avg_pident, lowest_evalue, query_count)

        # Safely get the two most common species
        if len(two_most_common) > 0:
            most_common_stats = get_stats(two_most_common[0][0])
        else:
            most_common_stats = ('none', 'none', 'none', 'none', 'none')

        if len(two_most_common) > 1:
            second_most_common_stats = get_stats(two_most_common[1][0])
        else:
            second_most_common_stats = ('none', 'none', 'none', 'none', 'none')

        return most_common_stats, second_most_common_stats

    most_common, second_most_common = parse_blast_output_with_stats(output_dir)
    reads = extract_blast_reads(file_path, start_string, end_string)

    plot_html = build_embedded_plot(output_dir)

    with open(os.path.join(output_dir, 'plot_only.html'), 'w', encoding='utf-8') as fplot:
        fplot.write(plot_html)

    html_start1 = """
    <li class="page-item"><a class="page-link" href="#!" onclick="showSection({0})">{0}</a></li>
    """
    html_start2 = """  \n    <div id="section{0}" class="content-section">\n    """
    html_end = """
        </div>\n    """

    # Fixed zero indexing on HTML report

    original_reads = reads.copy()
    for count, i in enumerate(reads.keys()):
        reads[i] = html_start1.format(count + 1)
    reads_list1 = '\n'.join(reads.values())

    for count, i in enumerate(original_reads.keys()):
        reads[i] = html_start2.format(count + 1) + '<pre>' + original_reads[i] + '</pre>' + html_end
    reads_list2 = '\n'.join(reads.values())

    barcode = fastq_dir.split('/')[-1]
    sample_id = output_dir.split('/')[3]

    report_dict = {
        "time": 'TEST' + " hrs",
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
        "plot_html": plot_html,
    }

    env = Environment(loader=FileSystemLoader('/mnt/db/ref/Template/organism_query'), autoescape=False)
    template = env.get_template('report_template.html')
    rendered_html = template.render({"report": report_dict})

    with open(os.path.join(output_dir, 'organism_report.html'), 'w', encoding='utf-8') as fh:
        fh.write(rendered_html)

    subprocess.Popen([
        "chromium-browser",
        "--noerrors",
        "--disable-session-crashed-bubble",
        "--disable-infobars",
        "--no-default-browser-check",
        "--no-first-run",
        "--log-level=3",
        "--disable-logging",
        "--ozone-platform=x11",
        os.path.join(output_dir, 'organism_report.html')
    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def build_embedded_plot(output_dir):
    file_path = os.path.join(output_dir, 'blast_results_6.tmp')

    col_names = [
        "Query read ID", "sseqid", "Percent Identity (%)", "staxids", "Scientific Name",
        "Alignment Length", "mismatch", "gapopen", "qstart", "qend", "sstart",
        "send", "evalue", "bitscore", "qseq", "qlen"
    ]

    try:
        df = pd.read_csv(file_path, sep="\t", names=col_names)
    except FileNotFoundError:
        return f"<p style='color:red;'>Could not find BLAST results file at: {file_path}</p>"
    except Exception as e:
        return f"<p style='color:red;'>Error reading BLAST output file: {e}</p>"

    if df.empty:
        return "<p style='color:red;'>BLAST output is empty; cannot create offline plot.</p>"

    fig = px.scatter(
        df,
        x='Percent Identity (%)',
        y='bitscore',
        color='Scientific Name',
        hover_data=['Scientific Name', 'Alignment Length'],
        log_y=True,
        title='Offline Interactive BLAST Results'
    )
    fig.update_layout(
        xaxis_title='Percent Identity (%)',
        yaxis_title='Bit score (log scale)'
    )
    return fig.to_html(full_html=False, include_plotlyjs=True)

def cleanup(output_dir):
    path_html = os.path.join(output_dir, 'blast_results.html')
    path_archive = os.path.join(output_dir, 'blast_results_11.tmp')
    if os.path.exists(path_html):
        os.remove(path_html)
    if os.path.exists(path_archive):
        os.remove(path_archive)

    path_fasta = os.path.join(output_dir, 'concatenated_subset_reads.fasta')
    if os.path.exists(path_fasta):
        subprocess.run(['gzip', path_fasta])

def decompress_gzip(file_path):
    decompressed_file_path = file_path.rstrip('.gz')
    if os.path.exists(decompressed_file_path):
        return decompressed_file_path
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


init()

def main():
    args = parse_args()

    # If --quiet, redirect stdout/stderr to /dev/null
    if args.quiet:
        devnull = open(os.devnull, 'w')
        sys.stdout = devnull
        sys.stderr = devnull

    if args.epi2me_report:
        input_format = "EPI2ME"
        time_interval = "EPI2ME"
        read_ids = extract_read_ids_epi2me(args.epi2me_report, args.organism, args.output_dir)
        matched_reads = []
    elif args.centrifuge_report_dir:
        input_format = "CIDR metag"
        raw_report_path = os.path.join(args.centrifuge_report_dir, 'centrifuge_raw_filtered.tsv')
        raw_report_gz_path = raw_report_path + '.gz'

        if os.path.exists(raw_report_path):
            raw_report = raw_report_path
        elif os.path.exists(raw_report_gz_path):
            raw_report = decompress_gzip(raw_report_gz_path)
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

        read_ids, matched_reads = extract_read_ids(
            raw_report,
            unique_taxids,
            args.output_dir,
            top_n=args.top_n,
            random_sort=args.random_sort
        )
    else:
        print("Error: No workflow analysis outputs were given.", flush=True)
        sys.exit(1)

    total_reads = extract_reads(args.fastq_dir, read_ids, args.output_dir)
    subset_reads = os.path.join(args.output_dir, 'concatenated_subset_reads.fastq')

    # Now pass in args.threads to run_blast
    run_blast(subset_reads, args.output_dir, args.blastdb, threads=args.threads)
    run_parser(subset_reads, args.output_dir, args.blastdb)

    report_build(
        args.output_dir,
        args.organism,
        read_ids,
        args.blastdb,
        total_reads,
        args.fastq_dir,
        input_format,
        time_interval,
        matched_reads
    )

    cleanup(args.output_dir)

if __name__ == '__main__':
    main()
