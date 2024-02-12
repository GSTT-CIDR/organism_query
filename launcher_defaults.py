import tkinter as tk
from tkinter import filedialog
from tkinter import ttk  # For combobox
from tkinter import Tk, Label
from tkinter import filedialog
from PIL import Image, ImageTk
import subprocess
import datetime
import os

current_datetime = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")

# Function to choose a directory
def choose_directory(entry):
    directory = filedialog.askdirectory()
    entry.delete(0, tk.END)
    entry.insert(0, directory)

# Function to choose a file
def choose_file(entry):
    filepath = filedialog.askopenfilename()
    if filepath:
        entry.delete(0, tk.END)
        entry.insert(0, filepath)

blast_db_dir = ""

def compile_and_launch():
    
    # Determine whether EPI2ME or CIDR input is used
    entry1_text = entry1.get()
    entry2_text = entry2.get()
    
    if entry1_text and entry2_text:
        messagebox.showerror("Error", "Both report entries are filled. Please fill only one.")
    elif entry1_text:
        # Clean organism name for output directory
        original_text = entry4.get()  # Get the text from entry4
        entry_4_clean = original_text.lower().replace(' ', '_').replace("'", "")
        # Build directory path and create directory
        directory_path = f"{entry5.get()}/organism_query_{entry_4_clean}_{current_datetime}/"
        os.makedirs(directory_path, exist_ok=True)
        # Arguments for organism query script
        args = [
            f"-c {entry1.get()}/",
            f"-f {entry3.get()}",
            f"-o {entry4.get()}",
            f"-d {directory_path}",
            f"-b {entry7.get()}/{dropdown_var.get()}"  # Add the selected
        ]
    elif entry2_text:
        # Clean organism name for output directory
        original_text = entry4.get()  # Get the text from entry4
        entry_4_clean = original_text.lower().replace(' ', '_').replace("'", "")
        # Output direcotry next to where the epi2me file is
        folder_path = os.path.dirname(entry2.get)
        directory_path = f"{entry5.get()}/organism_query_{entry_4_clean}_{current_datetime}/"
        os.makedirs(directory_path, exist_ok=True)
        args = [
            f"-e {entry2.get()}",
            f"-f {entry3.get()}",
            f"-o {entry4.get()}",
            f"-d {directory_path}",
            f"-b {entry7.get()}/{dropdown_var.get()}"  # Add the selected
        ]
    else:
        messagebox.showinfo("Result", "Query fields incomplete.")
    # Form the final command to be executed
    command = f"python organism_report.py {' '.join(args)}"
    # Execute the command here, e.g., using os.system or subprocess.run
    os.environ['NCBI_CONFIG_OVERRIDES'] = f"TRUE"
    os.environ['BLASTDB'] = f"{entry7.get()}"
    subprocess.run(command, shell=True)


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
description2 = tk.Label(root, text="Fill out the fields below with the files required to launch the query.", wraplength=430)
description2.pack()
description1 = tk.Label(root, text="")
description1.pack()
description3 = tk.Label(root, text="You should only fill out one of the report input fields. If you are querying an EPI2ME report, navigate to the CSV downloaded from the EPI2ME website.", wraplength=430)
description3.pack()


# Field 1
frame1 = tk.Frame(root)
frame1.pack(pady=5)
label1 = tk.Label(frame1, text="CIDR 'workflow/results/sampleID/XXhour/centrifuge/' folder")
label1.pack(side=tk.TOP)
entry1 = tk.Entry(frame1, width=40)
entry1.pack(side=tk.LEFT, padx=5)
buttton1 = tk.Button(frame1, text="Choose directory", command=lambda: choose_directory(entry1))
buttton1.pack(side=tk.LEFT)

# Field 1 - dropdown
#frame6 = tk.Frame(root)
#frame6.pack(pady=5)
#label6 = tk.Label(frame6, text="CIDR workflow hour/interval")
#label6.pack(side=tk.TOP, padx=5)  # Adjust the side to LEFT to align #with the dropdown
#options = ['0.5', '1', '2', '16', '24']
#combobox6 = ttk.Combobox(frame6, values=options, width=37)
#combobox6.pack(side=tk.RIGHT, padx=5)  # Adjust according to your layout needs

# Field 2
frame2 = tk.Frame(root)
frame2.pack(padx=10)
label2 = tk.Label(frame2, text="EPI2ME WIMP CSV")
label2.pack(side=tk.TOP)
entry2 = tk.Entry(frame2, width=40)
entry2.pack(side=tk.LEFT, padx=5)
button2 = tk.Button(frame2, text="Choose File", command=lambda: choose_file(entry2))
button2.pack(side=tk.LEFT)

# Field 3
frame3 = tk.Frame(root)
frame3.pack(padx=10)
label3 = tk.Label(frame3, text="MinKNOW barcode FASTQ directory")
label3.pack(side=tk.TOP)
entry3 = tk.Entry(frame3, width=40)
entry3.pack(side=tk.LEFT, padx=5)
buttton3 = tk.Button(frame3, text="Choose directory", command=lambda: choose_directory(entry3))
buttton3.pack(side=tk.LEFT)

#Field 4
frame4 = tk.Frame(root)
frame4.pack(padx=10)
label4 = tk.Label(frame4, text="Species name")
label4.pack(side=tk.TOP)
entry4 = tk.Entry(frame4, width=40)
entry4.pack(side=tk.LEFT, padx=5)

# Field 7
frame7 = tk.Frame(root)
frame7.pack(padx=10)
label7 = tk.Label(frame7, text="BLASTDB directory")
label7.pack(side=tk.TOP)
entry7 = tk.Entry(frame7, width=40)
entry7.pack(side=tk.LEFT, padx=5)
buttton7 = tk.Button(frame7, text="Choose directory", command=lambda: choose_directory(entry7))
buttton7.pack(side=tk.LEFT)

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

# Field 5
frame5 = tk.Frame(root)
frame5.pack(padx=10)
label5 = tk.Label(frame5, text="Output directory")
label5.pack(side=tk.TOP)
entry5 = tk.Entry(frame5, width=40)
entry5.pack(side=tk.LEFT, padx=5)
buttton5 = tk.Button(frame5, text="Choose output directory", command=lambda: choose_directory(entry5))
buttton5.pack(side=tk.LEFT)


# 'Launch Script' button
launch_button = tk.Button(root, text="Launch Script", command=compile_and_launch)
launch_button.pack(pady=10)

# Run the main application loop
root.mainloop()
