import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox, scrolledtext
from tkinter import ttk  # For combobox
from tkinter import Tk, Label
from tkinter import filedialog
from PIL import Image, ImageTk
import subprocess
import datetime
import os
from subprocess import Popen, PIPE, STDOUT
from threading import Thread

# Launcher for CIDR Organism Query 

current_datetime = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
# For initialising the subprocess state
process = None

# Function to choose a directory with a default path
def choose_directory_fastq(entry, default_path='/data/'):
    # Use the 'initialdir' option to open the dialog in the default path
    directory = filedialog.askdirectory(initialdir=default_path)
    entry.delete(0, tk.END)
    entry.insert(0, directory)
    
# Function to choose a file
def choose_file(entry):
    filepath = filedialog.askopenfilename()
    if filepath:
        entry.delete(0, tk.END)
        entry.insert(0, filepath)

def get_entry4_with_quotes():
    """
    Retrieves the text from entry4 and encases it in single quotes.
    """
    return "'" + entry4.get() + "'"

# The BLASTDB is determined by an environment variable in the .ncbirc
blast_db_dir = ""

# Builds arguments to pass to the organism query script.
def compile_and_launch():
    global process  # Make the process global to access it later for termination
    # Determine whether EPI2ME or CIDR input is used
    entry1_text = entry1.get()
    
    if entry1_text:
        # Clean organism name for output directory
        original_text = get_entry4_with_quotes()  # Get the text from entry4
        entry_4_clean = original_text.lower().replace(' ', '_').replace("'", "")
        # Build directory path and create directory
        directory_path = f"/mnt/reports/{entry1.get()}/organism_query_{entry_4_clean}_{current_datetime}/"
        os.makedirs(directory_path, exist_ok=True)
        # Arguments for organism query script
        args = [
            f"-c /mnt/results/{entry1.get()}/{combobox6.get()}_hours/centrifuge/",
            f"-f /mnt/results/{entry1.get()}/{combobox6.get()}_hours/microbial/",
            f"-o {get_entry4_with_quotes()}",
            f"-d /mnt/reports/{entry1.get()}/organism_query_{entry_4_clean}_{current_datetime}/",
            f"-b /mnt/db/blastdb/{dropdown_var.get()}"  # Add the selected option from dropdown
        ]
    else:
        messagebox.showinfo("Result", "Query fields incomplete.")
    # Form the final command to be executed
    # Have added here '-u' to allow for terminal output in tkinter
    command = f"python -u organism_report.py {' '.join(args)}"
    # Points to the blast DB and the taxamap
    os.environ['NCBI_CONFIG_OVERRIDES'] = f"TRUE"
    os.environ['BLASTDB'] = f"/mnt/db/blastdb"
    # Execute the command here, e.g., using os.system or subprocess.run
    # Modified part to execute command and capture output
    process = Popen(command, stdout=PIPE, stderr=STDOUT, shell=True, text=True, bufsize=1, universal_newlines=True)
    Thread(target=read_output, args=(process, terminal_output), daemon=True).start()

def read_output(process, output_widget):
    """Read and display both stdout and stderr from the process."""
    for line in iter(process.stdout.readline, ''):
        output_widget.insert(tk.END, line)
        output_widget.see(tk.END)
    process.stdout.close()

def on_closing():
    global process
    if process is not None and process.poll() is None:
        process.terminate()
    root.destroy()

# Create main window
root = tk.Tk()
root.title("CIDR BLAST organism query launcher")

# Add description and image
# Open the image file for header
image = Image.open("template/CIDR_logo.png")
# Resize the image
resized_image = image.resize((432, 136))
# Convert the image for Tkinter
tk_image = ImageTk.PhotoImage(resized_image)
label = Label(root, image=tk_image)
label.pack()
# Text section
description1 = tk.Label(root, text="")
description1.pack()
description2 = tk.Label(root, text="Fill out the fields below with the information contained within a CIDR Metagenomics PDF report to BLAST classified reads. A window will be automatically opened showing the results.", wraplength=430)
description2.pack()
description1 = tk.Label(root, text="")
description1.pack()

# Field 1
frame1 = tk.Frame(root)
frame1.pack(pady=5)
label1 = tk.Label(frame1, text="CIDR worklow Lab/sample ID")
label1.pack(side=tk.TOP)
entry1 = tk.Entry(frame1, width=40)
entry1.pack(side=tk.LEFT, padx=5)

# Field 6 - dropdown
frame6 = tk.Frame(root)
frame6.pack(pady=5)
label6 = tk.Label(frame6, text="Metagenomics workflow report hour/interval")
label6.pack(side=tk.TOP, padx=5)  # Adjust the side to LEFT to align with the dropdown
options = ['0.5', '1', '2', '16', '24']
combobox6 = ttk.Combobox(frame6, values=options, width=37)
combobox6.pack(side=tk.RIGHT, padx=5)  # Adjust according to your layout needs

# Species name input
frame4 = tk.Frame(root)
frame4.pack(padx=10)
label4 = tk.Label(frame4, text="Organism keyword")
label4.pack(side=tk.TOP)
entry4 = tk.Entry(frame4, width=40)
entry4.pack(side=tk.LEFT, padx=5)

# Dropdown (Combobox) for multiple choice input
frame_dropdown = tk.Frame(root)
frame_dropdown.pack(pady=5)
label_dropdown = tk.Label(frame_dropdown, text="Choose a BLAST database:")
label_dropdown.pack(side=tk.LEFT, padx=5)

dropdown_var = tk.StringVar()
dropdown = ttk.Combobox(frame_dropdown, textvariable=dropdown_var, state="readonly")
dropdown['values'] = ('nt', 'refseq', 'ref_prok_rep_genomes', 'ref_viruses_rep_genomes', 'ref_euk_rep_genomes')  # Hardcoded options
dropdown.current(0)  # Set the default value
dropdown.pack(side=tk.LEFT)

# 'Launch Script' button
launch_button = tk.Button(root, text="Launch Script", command=compile_and_launch)
launch_button.pack(pady=10)

# Adding a terminal-like output area at the bottom
terminal_output = scrolledtext.ScrolledText(root, height=20)
terminal_output.pack(fill=tk.BOTH, expand=True)

# Ensure the subprocess is terminated when the window is closed
root.protocol("WM_DELETE_WINDOW", on_closing)

# Run the main application loop
root.mainloop()
