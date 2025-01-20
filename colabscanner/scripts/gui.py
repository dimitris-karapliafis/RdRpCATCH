
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
import os
import subprocess
from tkinter import messagebox


class colabscanner_gui:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("ColabScanner GUI")
        self.root.geometry("1000x800")


        style = ttk.Style(self.root)
        # set theme style
        style.theme_use("clam")

        # Input file selection
        self.input_label = ttk.Label(self.root, text="Select input .fasta file:")
        self.input_label.grid(row=0, column=0, padx=5, pady=5)
        self.input_entry = ttk.Entry(self.root, width=50)
        self.input_entry.grid(row=0, column=1, padx=5, pady=5)
        self.input_browse_button = ttk.Button(self.root, text="Browse", command=self.browse_input_file)
        self.input_browse_button.grid(row=0, column=2, padx=5, pady=5)

        # Parent directory selection
        self.parent_dir_label = ttk.Label(self.root, text="Specify the directory to store output files:")
        self.parent_dir_label.grid(row=1, column=0, padx=5, pady=5)
        self.parent_dir_entry = ttk.Entry(self.root, width=50)
        self.parent_dir_entry.grid(row=1, column=1, padx=5, pady=5)
        self.parent_dir_browse_button = ttk.Button(self.root, text="Browse", command=self.browse_parent_directory)
        self.parent_dir_browse_button.grid(row=1, column=2, padx=5, pady=5)

        # Output directory name
        self.output_name_label = ttk.Label(self.root, text="Specify output Directory Name:")
        self.output_name_label.grid(row=2, column=0, padx=5, pady=5)
        self.output_name_entry = ttk.Entry(self.root, width=30)
        self.output_name_entry.grid(row=2, column=1, padx=5, pady=5)

        # Download options label
        self.download_label = ttk.Label(self.root, text="If databases are not downloaded, specify destination directory \n"
                                                       "and click 'Download Databases' (This might take a while): ")
        self.download_label.grid(row=3, column=0, padx=5, pady=5)

        # Download destination directory
        self.destination_entry = ttk.Entry(self.root, width=50)
        self.destination_entry.grid(row=3, column=1, padx=5, pady=5)
        self.destination_browse_button = ttk.Button(self.root, text="Browse", command=self.browse_destination_directory)
        self.destination_browse_button.grid(row=3, column=2, padx=5, pady=5)

        # Download button
        self.download_button = ttk.Button(self.root, text="Download Databases", command=self.download_databases)
        self.download_button.grid(row=4, column=1, padx=5, pady=5)


        # hmm_db directory selection
        self.hmm_dir_label = ttk.Label(self.root, text="Specify the directory where the HMM databases are stored:")
        self.hmm_dir_label.grid(row=5, column=0, padx=5, pady=5)
        self.hmm_dir_entry = ttk.Entry(self.root, width=50)
        self.hmm_dir_entry.grid(row=5, column=1, padx=5, pady=5)
        self.hmm_dir_browse_button = ttk.Button(self.root, text="Browse", command=self.browse_hmm_directory)
        self.hmm_dir_browse_button.grid(row=5, column=2, padx=5, pady=5)


        # Database selection frame (left side)
        self.db_frame = ttk.Frame(self.root)
        self.db_frame.grid(row=6, column=0, padx=5, pady=10)

        self.db_label = ttk.Label(self.db_frame, text="Select RdRp HMM databases to search against:")
        self.db_label.pack(anchor=tk.W)

        self.databases = {
            'RVMT': tk.BooleanVar(),
            'NeoRdRp': tk.BooleanVar(),
            'NeoRdRp.2.1': tk.BooleanVar(),
            'TSA_Olendraite': tk.BooleanVar(),
            'RDRP-scan': tk.BooleanVar(),
            'Lucaprot': tk.BooleanVar()
        }

        for db_name in self.databases:
            checkbox = ttk.Checkbutton(self.db_frame, text=db_name, variable=self.databases[db_name])
            checkbox.pack(anchor=tk.W)

        # HMMsearch parameters frame (right side)
        self.hmmsearch_frame = ttk.Frame(self.root)
        self.hmmsearch_frame.grid(row=6, column=1, padx=5, pady=10)

        self.hmmsearch_label = ttk.Label(self.hmmsearch_frame, text="HMMsearch Parameters:")
        self.hmmsearch_label.pack(anchor=tk.W)

        self.evalue_var = tk.StringVar(value='1e-05')
        self.inc_evalue_var = tk.StringVar(value='1e-05')
        self.dom_evalue_var = tk.StringVar(value='1e-05')
        self.incdom_evalue_var = tk.StringVar(value='1e-05')
        self.z_value_var = tk.IntVar(value=1000000)
        self.cpus_var = tk.IntVar(value=1)

        hmmsearch_params = [
            ('E-value threshold', self.evalue_var),
            ('Inclusion E-value threshold', self.inc_evalue_var),
            ('Domain E-value threshold', self.dom_evalue_var),
            ('Inclusion domain E-value threshold', self.incdom_evalue_var),
            ('Number of sequences to search against', self.z_value_var),
            ('Number of CPUs to use', self.cpus_var)
        ]

        for param_name, param_var in hmmsearch_params:
            frame = ttk.Frame(self.hmmsearch_frame)
            frame.pack(fill=tk.X, padx=5, pady=2)
            ttk.Label(frame, text=param_name).pack(side=tk.LEFT)
            ttk.Entry(frame, textvariable=param_var, width=15).pack(side=tk.RIGHT)

        # seqkit translate parameters
        self.seqkit_frame = ttk.Frame(self.root)
        self.seqkit_frame.grid(row=6, column=2, padx=5, pady=10)
        self.gen_code_var = tk.IntVar(value=1)
        self.frame_var = tk.IntVar(value=6)
        self.seqkit_translate_label = ttk.Label(self.seqkit_frame, text="Seqkit Translate Parameters:")
        self.seqkit_translate_label.pack(anchor=tk.W)
        seqkit_transl_params = [ ('Genetic code', self.gen_code_var), ('Frame (6: All frames)', self.frame_var)]
        for param_name, param_var in seqkit_transl_params:
            frame = ttk.Frame(self.seqkit_frame)
            frame.pack(fill=tk.X, padx=5, pady=2)
            ttk.Label(frame, text=param_name).pack(side=tk.LEFT)
            ttk.Entry(frame, textvariable=param_var, width=15).pack(side=tk.RIGHT)

        # Run button
        self.run_button = ttk.Button(self.root, text="Run Script", command=self.run_script)
        self.run_button.grid(row=7, column=1, padx=5, pady=5)

        # Status label
        self.status_label = tk.Label(self.root, text="", wraplength=400)
        self.status_label.grid(row=9, column=0, columnspan=3, padx=5, pady=5)

        # Progress text box (wider version)
        self.progress_frame = ttk.Frame(self.root)
        self.progress_frame.grid(row=8, column=0, columnspan=3, padx=5, pady=10)

        self.progress_label = ttk.Label(self.progress_frame, text="Progress:")
        self.progress_label.pack(anchor=tk.W)

        # Make the text area wider
        self.progress_text = tk.Text(self.progress_frame, height=15, width=120)
        self.progress_text.pack(fill=tk.BOTH, expand=True)
        self.progress_text.config(state='disabled')

        # Add horizontal scrollbar
        self.scrollbar_x = ttk.Scrollbar(self.progress_frame, orient=tk.HORIZONTAL)
        self.scrollbar_x.pack(side=tk.BOTTOM, fill=tk.X)
        self.progress_text.config(xscrollcommand=self.scrollbar_x.set)
        self.scrollbar_x.config(command=self.progress_text.xview)

        # Vertical scrollbar (already added earlier)
        self.scrollbar = ttk.Scrollbar(self.progress_frame)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.progress_text.config(yscrollcommand=self.scrollbar.set)
        self.scrollbar.config(command=self.progress_text.yview)

    def browse_input_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("All Files", "*.*"), ("Text Files", "*.txt")])
        self.input_entry.delete(0, tk.END)
        self.input_entry.insert(tk.END, file_path)

    def browse_parent_directory(self):
        directory = filedialog.askdirectory()
        self.parent_dir_entry.delete(0, tk.END)
        self.parent_dir_entry.insert(tk.END, directory)

    def browse_hmm_directory(self):
        directory = filedialog.askdirectory()
        self.hmm_dir_entry.delete(0, tk.END)
        self.hmm_dir_entry.insert(tk.END, directory)

    def browse_destination_directory(self):
        directory = filedialog.askdirectory()
        self.destination_entry.delete(0, tk.END)
        self.destination_entry.insert(tk.END, directory)

    def download_databases(self):
        destination = self.destination_entry.get()
        if not destination:
            self.status_label.config(text="Error: Destination directory must be specified.", fg="red")
            return

        try:
            result = subprocess.run(["python3", "colabscanner.py", "-download", "-destination", destination],
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)

            self.status_label.config(text=f"Database download completed successfully.\nOutput:\n{result.stdout}",
                                     fg="green")
        except subprocess.CalledProcessError as e:
            with open("error_log.txt", "w") as log_file:
                log_file.write(e.output)
            self.status_label.config(
                text="Database download failed. See 'error_log.txt' for details.", fg="red"
            )
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred during database download:\n{str(e)}")
            print(f"Unexpected error: {e}")

    def run_script(self):
        input_file = self.input_entry.get()
        parent_dir = self.parent_dir_entry.get()
        output_name = self.output_name_entry.get()
        hmm_dir = self.hmm_dir_entry.get()

        selected_dbs = [db for db, var in self.databases.items() if var.get()]

        if not os.path.isfile(input_file):
            self.status_label.config(text="Error: Input file does not exist.", fg="red")
            return

        if not os.path.isdir(parent_dir):
            self.status_label.config(text="Error: Parent directory does not exist.", fg="red")
            return

        output_dir = os.path.join(parent_dir, output_name)

        try:
            os.makedirs(output_dir, exist_ok=True)
        except Exception as e:
            self.status_label.config(text=f"Error creating output directory: {str(e)}", fg="red")
            return

        db_options = ",".join(selected_dbs) if selected_dbs else "all"

        hmmsearch_args = [
            f'-e', self.evalue_var.get(),
            f'-incE', self.inc_evalue_var.get(),
            f'-domE', self.dom_evalue_var.get(),
            f'-incdomE', self.incdom_evalue_var.get(),
            f'-z', str(self.z_value_var.get()),
            f'-cpus', str(self.cpus_var.get())
        ]

        try:
            # Clear previous progress
            self.clear_progress()

            # Run the script
            result = subprocess.Popen(
                ["python3", "colabscanner.py", "-i", input_file, "-o", output_dir, "-hmm_dir", hmm_dir, "-dbs",
                 db_options] + hmmsearch_args,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)

            # Read and display output line by line
            for line in iter(result.stdout.readline, ""):
                self.update_progress(line.strip())

            result.wait()
            if result.returncode == 0:
                self.status_label.config(text="Script executed successfully.", fg="green")
            else:
                self.status_label.config(text="Script execution failed.", fg="red")

        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred:\n{str(e)}")
            print(f"Unexpected error: {e}")

    def clear_progress(self):
        self.progress_text.config(state='normal')
        self.progress_text.delete('1.0', tk.END)
        self.progress_text.config(state='disabled')

    def update_progress(self, message):
        self.progress_text.config(state='normal')
        self.progress_text.insert(tk.END, message + "\n")
        self.progress_text.yview(tk.END)
        self.progress_text.config(state='disabled')
        self.root.update_idletasks()
    def run(self):
        self.root.mainloop()


if __name__ == "__main__":
    gui = colabscanner_gui()
    gui.run()


