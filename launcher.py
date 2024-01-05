import tkinter as tk
from tkinter import filedialog
from tkinter import ttk  # For combobox

def choose_directory(entry):
    directory = filedialog.askdirectory()
    entry.delete(0, tk.END)
    entry.insert(0, directory)

def compile_and_launch():
    args = [
        f"--flag1 {entry1.get()}",
        f"--flag2 {entry2.get()}",
        f"--flag3 {entry3.get()}",
        f"--flag4 {entry4.get()}",
        f"--flag5 {entry5.get()}",
        f"--option {dropdown_var.get()}"  # Add the selected option from dropdown
    ]
    command = f"python path_to_other_script.py {' '.join(args)}"
    # Execute the command here, e.g., using os.system or subprocess.run

# Create main window
root = tk.Tk()
root.title("Script Launcher")

# Add description and image
description = tk.Label(root, text="This is a launcher for your script.")
description.pack()

image = tk.PhotoImage(file="./organism_query/template/CIDR_logo.png")  # Update with the path to your image
image_label = tk.Label(root, image=image)
image_label.pack()

# Field 1
frame1 = tk.Frame(root)
frame1.pack(pady=5)
label1 = tk.Label(frame1, text="Description for Field 1")
label1.pack(side=tk.TOP)
entry1 = tk.Entry(frame1, width=40)
entry1.pack(side=tk.LEFT, padx=5)
button1 = tk.Button(frame1, text="Choose Directory", command=lambda: choose_directory(entry1))
button1.pack(side=tk.LEFT)

# Field 2
frame2 = tk.Frame(root)
frame2.pack(padx=10)
label2 = tk.Label(frame2, text="Description for Field 2")
label2.pack(side=tk.TOP)
entry2 = tk.Entry(frame2, width=40)
entry2.pack(side=tk.LEFT, padx=5)
button2 = tk.Button(frame2, text="Choose Directory", command=lambda: choose_directory(entry2))
button2.pack(side=tk.LEFT)

# [Repeat for Field 3, Field 4, and Field 5]

# Dropdown (Combobox) for multiple choice input
frame_dropdown = tk.Frame(root)
frame_dropdown.pack(pady=5)

label_dropdown = tk.Label(frame_dropdown, text="Choose an option:")
label_dropdown.pack(side=tk.LEFT, padx=5)

dropdown_var = tk.StringVar()
dropdown = ttk.Combobox(frame_dropdown, textvariable=dropdown_var, state="readonly")
dropdown['values'] = ('Option1', 'Option2', 'Option3')  # Hardcoded options
dropdown.current(0)  # Set the default value
dropdown.pack(side=tk.LEFT)

# 'Launch Script' button
launch_button = tk.Button(root, text="Launch Script", command=compile_and_launch)
launch_button.pack(pady=10)

# Run the main application loop
root.mainloop()
