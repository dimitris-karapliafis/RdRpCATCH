import logging
import time
from rich.console import Console
from Bio import SeqIO
import os
import pandas as pd

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

   def loud_log(self, message):
       self.console.log(message)
       self.logger.info(message)

   def silent_log(self, message):
       self.logger.info(message)

   def start_timer(self):
       self.start_time = time.time()

   def stop_timer(self, verbose=None):
        end_time = time.time()
        raw_execution_time = end_time - self.start_time
        if raw_execution_time < 60:
            execution_time = f"{raw_execution_time :.2f} seconds"
        elif raw_execution_time < 3600:
            execution_time = f"{raw_execution_time/60 :.2f} minutes"
        else:
            execution_time = f"{raw_execution_time/3600 :.2f} hours"

        return execution_time


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
                seq_type = 'nuc'
            elif set(seq).issubset(set(dna_set_ambiguous)):
                seq_type = 'nuc'
            elif set(seq).issubset(protein_set):
                seq_type = 'prot'
            else:
                raise Exception(f"Invalid sequence type in fasta file: {self.fasta_file} for sequence: {header} with sequence: {set(seq)}")

        return seq_type

    def check_seq_length(self, max_len):
        """
        Check the length of sequences in a FASTA file using Biopython.

        Args:
            fasta_file_path (str): Path to the FASTA file.
            max_allowed_length (int): Maximum allowed length for sequences.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If any sequence exceeds the maximum allowed length.
        """
        ## Check if the file exists
        if not os.path.isfile(self.fasta_file):
            raise FileNotFoundError(f"The file '{self.fasta_file}' does not exist.")

        ## Parse the FASTA file
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            if len(record.seq) > max_len:
                raise ValueError(f"Sequence ID: {record.id}, Description: {record.description}. "
                                 f"Length: {len(record.seq)}, Exceeds maximum allowed length: {max_len}. Please check the input file, "
                                 f"as this will cause issues with the pyHMMER search.")

        return True



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

    def write_fasta_coords(self, coords_list, output_file, seq_type):


        if seq_type == 'nuc':
            contigs = [coords[1] for coords in coords_list]
        else:
            contigs = [coords[0] for coords in coords_list]

        record_list = self.extract_contigs(contigs)

        with open(output_file, 'w') as f:
            for record in record_list:
                for coord in coords_list:
                    if seq_type == 'nuc':
                        if record.id == coord[1]:
                            record.id = record.id + f"_RdRp_start:{coord[2]}_end:{coord[3]}"
                            SeqIO.write(record, f, "fasta")
                    elif seq_type == 'prot':
                        if record.id == coord[0]:
                            record.id = record.id + f"_RdRp_start:{coord[2]}_end:{coord[3]}"
                            SeqIO.write(record, f, "fasta")
                    else:
                        raise Exception(f"Invalid sequence type: {seq_type}")
        return output_file


class mmseqs_parser:

    def __init__(self, mmseqs_tax_out_file, mmseqs_s_out_file):

        self.mmseqs_tax_out_file = mmseqs_tax_out_file
        self.mmseqs_s_out_file = mmseqs_s_out_file


    def parse_mmseqs_tax_lca(self):
        with open(self.mmseqs_tax_out_file, 'r') as f:
            lca_dict = {}
            for line in f:
                line = line.strip().split('\t')

                contig = line[0]
                if len(line) < 5:
                    lca_lineage = line[3]
                else:
                    lca_lineage = line[4]
                lca_dict[contig] = lca_lineage
        return lca_dict

    def parse_mmseqs_e_search_tophit(self):

        with open(self.mmseqs_s_out_file, 'r') as f:
            tophit_dict = {}
            for line in f:
                line = line.strip().split('\t')
                contig = line[0]

                if contig not in tophit_dict:
                    target = line[1]
                    fident = line[2]
                    alnlen = line[3]
                    eval = line[10]
                    bits = line[11]
                    qcov = line[12]
                    lineage = line[14]
                    tophit_dict[contig] = [target,fident, alnlen, eval, bits, qcov, lineage]
                else:
                    continue

        return tophit_dict

    def tax_to_rdrpcatch(self, rdrpcatch_out, extended_rdrpcatch_out, seq_type):

        lca_dict = self.parse_mmseqs_tax_lca()
        tophit_dict = self.parse_mmseqs_e_search_tophit()

        with open(rdrpcatch_out, "r") as f_handle, open(extended_rdrpcatch_out, 'w') as out_handle:

            for line in f_handle:
                line = line.strip().split("\t")
                if line[0].startswith("#Contig_name"):
                    title_line = line + ["MMseqs_Taxonomy_2bLCA", "MMseqs_TopHit_accession", "MMseqs_TopHit_fident",
                                       "MMseqs_TopHit_alnlen", "MMseqs_TopHit_eval", "MMseqs_TopHit_bitscore",
                                       "MMseqs_TopHit_qcov", "MMseqs_TopHit_lineage"]
                    out_handle.write("\t".join(title_line) + "\n")
                else:

                    if seq_type == 'nuc':
                        c_name  = line[1]
                    else:
                        c_name = line[0]

                    full_c_name = f"{c_name}_RdRp_start:{line[9]}_end:{line[10]}"

                    if full_c_name in lca_dict:
                        lca = lca_dict[full_c_name]
                    else:
                        lca = 'NA'

                    if full_c_name in tophit_dict:
                        tophit = tophit_dict[full_c_name]
                        tophit_accession = tophit[0]
                        tophit_fident = tophit[1]
                        tophit_alnlen = tophit[2]
                        tophit_eval = tophit[3]
                        tophit_bitscore = tophit[4]
                        tophit_qcov = tophit[5]
                        tophit_lineage = tophit[6]
                    else:
                        tophit_accession = 'NA'
                        tophit_fident = 'NA'
                        tophit_alnlen = 'NA'
                        tophit_eval = 'NA'
                        tophit_bitscore = 'NA'
                        tophit_qcov = 'NA'
                        tophit_lineage = 'NA'

                    line.extend([lca, tophit_accession, tophit_fident, tophit_alnlen, tophit_eval, tophit_bitscore,
                                    tophit_qcov, tophit_lineage])

                    out_handle.write("\t".join(line) + "\n")

        df = pd.read_csv(extended_rdrpcatch_out, sep="\t")
        df = df.drop(["Best_hit_norm_bitscore_profile", "Best_hit_norm_bitscore_contig","Best_hit_ID_score", "Best_hit_aln_length"], axis=1)
        column_order = ["#Contig_name","Translated_contig_name (frame)",
                        "Sequence_length(AA)","Total_databases_that_the_contig_was_detected(No_of_Profiles)",
                        "Best_hit_Database","Best_hit_profile_name", "Best_hit_profile_length", "Best_hit_e-value",
                        "Best_hit_bitscore","RdRp_from(AA)", "RdRp_to(AA)","Best_hit_profile_coverage",
                        "Best_hit_contig_coverage", "MMseqs_Taxonomy_2bLCA", "MMseqs_TopHit_accession",
                        "MMseqs_TopHit_fident", "MMseqs_TopHit_alnlen", "MMseqs_TopHit_eval",
                        "MMseqs_TopHit_bitscore", "MMseqs_TopHit_qcov", "MMseqs_TopHit_lineage"]
        df = df[column_order]
        df.to_csv(extended_rdrpcatch_out, sep="\t", index=False)





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










