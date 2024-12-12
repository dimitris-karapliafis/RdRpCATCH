import logging
import time
from rich.console import Console
from Bio import SeqIO
import os

class Logger:
   def __init__(self, log_file):
       self.console = Console()
       self.log_file = log_file
       self.logger = logging.getLogger('Logger')
       self.logger.setLevel(logging.INFO)
       handler = logging.FileHandler(self.log_file)
       handler.setLevel(logging.INFO)
       formatter = logging.Formatter('%(asctime)s - %(message)s')
       handler.setFormatter(formatter)
       self.logger.addHandler(handler)

   def log(self, message):
       self.console.log(message)
       self.logger.info(message)

   def silent_log(self, message):
       self.logger.info(message)

   def start_timer(self):
       self.start_time = time.time()

   def stop_timer(self):
       end_time = time.time()
       execution_time = end_time - self.start_time
       self.log(f"Execution time: {execution_time} seconds")



class fasta_checker:

    def __init__(self, fasta_file):
        self.fasta_file = fasta_file


    def check_fasta_validity(self):

        with open(self.fasta_file, 'r') as f:
            first_line = f.readline()
            if not first_line.startswith('>'):
                raise Exception(f"Invalid fasta file: {self.fasta_file}, first line is {first_line}")
            else:
                return True

    def read_fasta(self):
        with open(self.fasta_file, 'r') as f:
            fasta_dict = {}
            for line in f:
                if line.startswith('>'):
                    header = line.strip()
                    fasta_dict[header] = ''
                else:
                    fasta_dict[header] += line.strip()
            return fasta_dict


    def check_seq_type(self):

        fasta_dict = self.read_fasta()

        seq_type = ''
        dna_set = {'A', 'T', 'G', 'C'}
        dna_set_ambiguous = {'A', 'T', 'G', 'C', 'N'}
        protein_set =  {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'}
        for header, seq in fasta_dict.items():
            if set(seq).issubset(dna_set):
                seq_type = 'DNA'
            elif set(seq).issubset(set(dna_set_ambiguous)):
                seq_type = 'DNA'
            elif set(seq).issubset(protein_set):
                seq_type = 'PROTEIN'
            else:
                raise Exception(f"Invalid sequence type in fasta file: {self.fasta_file} for sequence: {header} with sequence: {set(seq)}")

        return seq_type



class fasta:

    def __init__(self, fasta_file):
        self.fasta_file = fasta_file


    def extract_contigs(self, contig_list):
        record_list = []
        with open(self.fasta_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                if record.id in contig_list:
                    record_list.append(record)

        return record_list

    def write_fasta(self, record_list, output_file):
        with open(output_file, 'w') as f:
            SeqIO.write(record_list, f, "fasta")
        return output_file



class file_handler:

    def __init__(self, file):
        self.file = file

    def check_file_exists(self):
        if not os.path.exists(self.file):
            raise Exception(f"File does not exist: {self.file}")
        return True

    def delete_file(self):
        os.remove(self.file)
        return True

    def check_file_size(self):
        return os.path.getsize(self.file)

    def check_file_extension(self):
        return os.path.splitext(self.file)[1]

    def get_file_name(self):
        return os.path.basename(self.file)

    def get_file_dir(self):
        return os.path.dirname(self.file)










