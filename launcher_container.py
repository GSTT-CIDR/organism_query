import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox, scrolledtext
from tkinter import ttk  # For combobox & Notebook
from tkinter import Tk, Label
from PIL import Image, ImageTk
import subprocess
import datetime
import os
from subprocess import Popen, PIPE, STDOUT
from threading import Thread
import csv  # For reading the replacement CSV

# Globals
current_datetime = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
processes = []  # Keep track of all Popen processes

# Hard-coded CSV file path for the lookup
CSV_LOOKUP_PATH = "/mnt/db/ref/reporting_name_replacement_list.csv"

def clean_organism_name(name: str) -> str:
    """
    Lowercase, replace spaces with underscores, remove single quotes, etc.,
    to create a directory-friendly string.
    """
    return name.lower().replace(' ', '_').replace("'", "")

def add_organism_field():
    """
    Dynamically add a new row with an Entry to capture an additional organism keyword.
    """
    frame = tk.Frame(frame_organisms)
    frame.pack(fill=tk.X, pady=2)

    entry = tk.Entry(frame, width=40)
    entry.pack(side=tk.LEFT, padx=5)
    organism_entries.append(entry)

def compile_and_launch():
    """
    Launch the script in parallel for each organism keyword. Each run of the script
    gets its own tab in the Notebook, with a dedicated ScrolledText for logging.
    """
    global processes
    processes = []

    # Clear any old tabs in the notebook
    for tab_id in notebook.tabs():
        notebook.forget(tab_id)

    # Basic checks
    sample_id = entry1.get().strip()
    if not sample_id:
        messagebox.showinfo("Error", "Please fill in the 'CIDR workflow Lab/sample ID' field.")
        return

    hours_val = combobox6.get().strip()
    if not hours_val:
        messagebox.showinfo("Error", "Please select a metagenomics workflow report hour/interval.")
        return

    # Attempt to load the CSV lookup from the hard-coded path
    lookup_dict = {}
    if os.path.isfile(CSV_LOOKUP_PATH):
        try:
            with open(CSV_LOOKUP_PATH, newline='') as f:
                reader = csv.DictReader(f)
                # Ensure columns exist
                if not {'Original','Replacement'}.issubset(reader.fieldnames):
                    messagebox.showerror("Error", "CSV must contain 'original' and 'replacement' columns.")
                    return
                for row in reader:
                    # For each row, map the entire replacement -> original
                    lookup_dict[row['Replacement'].strip()] = row['Original'].strip()
        except Exception as e:
            messagebox.showerror("Error", f"Could not read CSV file '{CSV_LOOKUP_PATH}': {e}")
            return
    else:
        # If the file doesn't exist, you could show an error or just proceed.
        # Let's show an error to let the user know.
        messagebox.showerror("Error", f"Could not find CSV file: {CSV_LOOKUP_PATH}")
        return

    # For each organism entry
    launched_any = False
    for entry in organism_entries:
        org_raw = entry.get().strip()
        if not org_raw:
            continue  # Skip blank entries

        # If the user typed a known 'replacement', revert it to 'original'
        if org_raw in lookup_dict:
            org_raw = lookup_dict[org_raw]

        launched_any = True
        # Build directory path
        org_clean = clean_organism_name(org_raw)
        directory_path = f"/mnt/reports/{sample_id}/organism_query_{org_clean}_{current_datetime}/"
        os.makedirs(directory_path, exist_ok=True)
        
        # Construct command args
        org_quoted = f"'{org_raw}'"
        args = [
            f"-c /mnt/results/{sample_id}/{hours_val}_hours/centrifuge/",
            f"-f /mnt/results/{sample_id}/{hours_val}_hours/microbial/",
            f"-o {org_quoted}",
            f"-d {directory_path}",
            f"-b /mnt/db/blastdb/{dropdown_var.get()}"
        ]
        
        command = f"python -u organism_report.py {' '.join(args)}"
        print(command)
        # Create a new tab in the Notebook for this organism
        tab_frame = ttk.Frame(notebook)
        notebook.add(tab_frame, text=org_raw)  # Tab label is the raw search term
        notebook.pack(fill=tk.BOTH, expand=True)

        # Create a ScrolledText in the tab
        st = scrolledtext.ScrolledText(tab_frame, height=20)
        st.pack(fill=tk.BOTH, expand=True)

        # Set environment variables for BLAST
        os.environ['NCBI_CONFIG_OVERRIDES'] = "TRUE"
        os.environ['BLASTDB'] = "/mnt/db/blastdb"

        # Launch the process
        proc = Popen(command, stdout=PIPE, stderr=STDOUT, shell=True, text=True, bufsize=1, universal_newlines=True)
        processes.append(proc)

        # Start a thread to read this process's output
        Thread(target=read_output, args=(proc, st), daemon=True).start()

    if not launched_any:
        messagebox.showinfo("Info", "No organism keywords entered.")


def read_output(process, output_widget):
    """
    Read and display stdout/stderr from the process in the given ScrolledText widget.
    """
    for line in iter(process.stdout.readline, ''):
        output_widget.insert(tk.END, line)
        output_widget.see(tk.END)
    process.stdout.close()


def on_closing():
    """
    Terminate all subprocesses if still running, then close the GUI.
    """
    global processes
    for proc in processes:
        if proc.poll() is None:
            proc.terminate()
    root.destroy()

# Create main window
root = tk.Tk()
root.title("CIDR BLAST organism query launcher")

# Set application icon if available
icon_path = "/mnt/lib/install/CIDR_logo_square_rmg.png"
if os.path.isfile(icon_path):
    icon_img = ImageTk.PhotoImage(file=icon_path)
    root.iconphoto(True, icon_img)

# Add description and image
try:
    image = Image.open("template/CIDR_logo.png")
    resized_image = image.resize((432, 136))
    tk_image = ImageTk.PhotoImage(resized_image)
    label = Label(root, image=tk_image)
    label.pack()
except Exception:
    # If the image isn't found, skip gracefully
    label = Label(root, text="CIDR Logo Here")
    label.pack()

description1 = tk.Label(root, text="")
description1.pack()
description2 = tk.Label(
    root,
    text=(
        "Fill out the fields below with information from a CIDR Metagenomics PDF report "
        "to BLAST classified reads. Each added organism keyword will run in parallel "
        "in its own tab below."
    ),
    wraplength=430
)
description2.pack()
description3 = tk.Label(root, text="")
description3.pack()

# Frame for sample ID
frame1 = tk.Frame(root)
frame1.pack(pady=5)
label1 = tk.Label(frame1, text="CIDR workflow Lab/sample ID")
label1.pack(side=tk.TOP)
entry1 = tk.Entry(frame1, width=40)
entry1.pack(side=tk.LEFT, padx=5)

# Frame for hours/interval
frame6 = tk.Frame(root)
frame6.pack(pady=5)
label6 = tk.Label(frame6, text="Metagenomics workflow report hour/interval")
label6.pack(side=tk.TOP, padx=5)
options = ['0.5', '1', '2', '16', '24']
combobox6 = ttk.Combobox(frame6, values=options, width=37, state="readonly")
combobox6.pack(side=tk.RIGHT, padx=5)

# Frame for organism keywords + plus button
frame_organisms = tk.Frame(root)
frame_organisms.pack(padx=10, pady=5)

# Row containing label + first organism entry + plus button
first_row = tk.Frame(frame_organisms)
first_row.pack(fill=tk.X, pady=2)
second_row = tk.Frame(frame_organisms)
second_row.pack(fill=tk.X, pady=2)

label4 = tk.Label(first_row, text="Organism keyword(s):")
label4.pack(side=tk.LEFT, padx=5)

# We'll keep a list of all Entry widgets for the organism keywords
organism_entries = []

# Create the first default entry
entry_default = tk.Entry(second_row, width=40)
entry_default.pack(side=tk.LEFT, padx=5)
organism_entries.append(entry_default)

# Plus button
plus_button = tk.Button(second_row, text="+", command=add_organism_field)
plus_button.pack(side=tk.LEFT, padx=5)

# Frame for choosing BLAST database
frame_dropdown = tk.Frame(root)
frame_dropdown.pack(pady=5)
label_dropdown = tk.Label(frame_dropdown, text="Choose a BLAST database:")
label_dropdown.pack(side=tk.LEFT, padx=5)

dropdown_var = tk.StringVar()
dropdown = ttk.Combobox(frame_dropdown, textvariable=dropdown_var, state="readonly")
dropdown['values'] = ('core_nt', 'nt', 'refseq', 'ref_prok_rep_genomes', 'ref_viruses_rep_genomes', 'ref_euk_rep_genomes', 'blastdb/v14_autoquery_release')
dropdown.current(0)
dropdown.pack(side=tk.LEFT)

# 'Launch Script' button
launch_button = tk.Button(root, text="Launch Script", command=compile_and_launch)
launch_button.pack(pady=10)

# Notebook for tabbed output
notebook = ttk.Notebook(root)
notebook.pack(fill=tk.BOTH, expand=True)

version = tk.Label(root, text="Organism Query v1.7.1 (RMg network release)", font=("Helvetica", 8))
version.pack(side="bottom", anchor="w", padx=10, pady=10)

# Ensure the subprocesses are terminated when the window is closed
root.protocol("WM_DELETE_WINDOW", on_closing)

# Run the main application loop
root.mainloop()
