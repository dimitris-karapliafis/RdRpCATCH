import logging
import time
from rich.console import Console
import os
import polars as pl
import needletail


def write_combined_results_to_gff(output_file, combined_data,seq_type):
    with open(output_file, 'w') as f:
        f.write("##gff-version 3\n")
        for row in combined_data.iter_rows(named=True):
            record = convert_record_to_gff3_record(row, seq_type)
            f.write(f"{record}\n")

def convert_record_to_gff3_record(row,seq_type): # for dict objects expected to be coherced into a gff3
    # taken from rolypoly https://code.jgi.doe.gov/UNeri/rolypoly/-/blob/main/src/rolypoly/commands/annotation/annotate_RNA.py
    
    # try to identify a sequence_id columns (query, qseqid, contig_id, contig, id, name)
    if seq_type == 'nuc':
        sequence_id_col = "Translated_contig_name (frame)"
    else:
        sequence_id_columns = ["sequence_id",'query', 'qseqid', 'contig_id', 'contig', 'id', 'name','Contig_name']
        sequence_id_col = next((col for col in sequence_id_columns if col in row.keys()), None)
        if sequence_id_col is None:
            raise ValueError(f"No sequence ID column found in row. Available columns: {list(row.keys())}")
    
    # try to identify a score column (score, Score, bitscore, qscore, bit)
    score_columns = ["score", "Score", "bitscore", "qscore", "bit","bits"]
    score_col = next((col for col in score_columns if col in row.keys()), "score")
    
    # try to identify a source column (source, Source, db, DB)
    source_columns = ["source", "Source", "db", "DB"]
    source_col = next((col for col in source_columns if col in row.keys()), "source")
    
    # try to identify a type column (type, Type, feature, Feature)
    type_columns = ["type", "Type", "feature", "Feature"]
    type_col = next((col for col in type_columns if col in row.keys()), "type")
    
    # try to identify a strand column (strand, Strand, sense, Sense)
    strand_columns = ["strand", "Strand", "sense", "Sense"]
    strand_col = next((col for col in strand_columns if col in row.keys()), "strand")
    
    # try to identify a phase column (phase, Phase)
    phase_columns = ["phase", "Phase"]
    phase_col = next((col for col in phase_columns if col in row.keys()), "phase")
    
    # Build GFF3 attributes string
    attrs = []
    for key, value in row.items():
        if key not in [sequence_id_col, source_col, score_col, type_col, strand_col, phase_col]:
            attrs.append(f"{key}={value}")
    
    # Get values, using defaults for missing columns
    sequence_id = row[sequence_id_col]
    source = row.get(source_col, "rdrpcatch")
    score = row.get(score_col, "0")
    feature_type = row.get(type_col, "feature")
    strand = row.get(strand_col, "+")
    phase = row.get(phase_col, ".")
    
    # Format GFF3 record
    gff3_fields = [
        sequence_id,
        source,
        feature_type,
        str(row.get("RdRp_from(AA)", "1")),
        str(row.get("RdRp_to(AA)", "1")),
        str(score),
        strand,
        phase,
        ";".join(attrs) if attrs else "."
    ]
    
    return "\t".join(gff3_fields)





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

       return self.start_time

   def stop_timer(self, start_time, verbose=None):
        end_time = time.time()

        raw_execution_time = end_time - start_time

        # Calculate hours, minutes, and seconds
        hours = int(raw_execution_time // 3600)
        minutes = int((raw_execution_time % 3600) // 60)
        seconds = int(raw_execution_time % 60)
        milliseconds = int((raw_execution_time % 1) * 1000)

        # Format the output
        execution_time = f"{hours} Hours {minutes} Minutes {seconds} Seconds {milliseconds} ms"

        return execution_time


class fasta_checker:

    def __init__(self, fasta_file, logger=None):
        self.fasta_file = fasta_file
        self.logger = logger

    def check_fasta_validity(self):
        reader = needletail.parse_fastx_file(self.fasta_file)
        try:
            first_record = next(reader)
            if self.logger:
                self.logger.silent_log(f"Successfully validated fasta file: {self.fasta_file}")
            return True
        except StopIteration:
            error_msg = f"Invalid or empty fasta file: {self.fasta_file}"
            if self.logger:
                self.logger.silent_log(error_msg)
            raise Exception(error_msg)
        except Exception as e:
            error_msg = f"Invalid fasta file: {self.fasta_file}, error: {str(e)}"
            if self.logger:
                self.logger.silent_log(error_msg)
            raise Exception(error_msg)

    def read_fasta(self):
        fasta_dict = {}
        reader = needletail.parse_fastx_file(self.fasta_file)
        for record in reader:
            header = f">{record.id}"
            fasta_dict[header] = record.seq
        if self.logger:
            self.logger.silent_log(f"Read {len(fasta_dict)} sequences from {self.fasta_file}")
        return fasta_dict

    def check_seq_type(self):
        reader = needletail.parse_fastx_file(self.fasta_file)
        dna_set = {'A', 'T', 'G', 'C'}
        dna_set_ambiguous = {'A', 'T', 'G', 'C', 'N'}
        protein_set = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'}
        
        for record in reader:
            seq = record.seq.upper()
            if set(seq).issubset(dna_set):
                if self.logger:
                    self.logger.silent_log(f"Detected nucleotide sequence (strict DNA alphabet)")
                return 'nuc'
            elif set(seq).issubset(dna_set_ambiguous):
                if self.logger:
                    self.logger.silent_log(f"Detected nucleotide sequence (ambiguous DNA alphabet)")
                return 'nuc'
            elif set(seq).issubset(protein_set):
                if self.logger:
                    self.logger.silent_log(f"Detected protein sequence")
                return 'prot'
            else:
                error_msg = f"Invalid sequence type in fasta file: {self.fasta_file} for sequence: {record.id.encode()} with sequence: {set(seq)}"
                if self.logger:
                    self.logger.silent_log(error_msg)
                raise Exception(error_msg)

    def check_seq_length(self, max_len):
        if not os.path.isfile(self.fasta_file):
            error_msg = f"The file '{self.fasta_file}' does not exist."
            if self.logger:
                self.logger.silent_log(error_msg)
            raise FileNotFoundError(error_msg)

        reader = needletail.parse_fastx_file(self.fasta_file)
        for record in reader:
            if len(record.seq) > max_len:
                error_msg = f"Sequence ID: {record.id}, Length: {len(record.seq)}, " \
                           f"Exceeds maximum allowed length: {max_len}. Please check the input file, " \
                           f"as this will cause issues with the pyHMMER search."
                if self.logger:
                    self.logger.silent_log(error_msg)
                raise ValueError(error_msg)
        if self.logger:
            self.logger.silent_log(f"All sequences are within length limit of {max_len}")
        return True



class fasta:

    def __init__(self, fasta_file, logger=None):
        self.fasta_file = fasta_file
        self.logger = logger


    def extract_contigs(self, contig_list):
        """
        Extract contigs from a fasta file based on a list of contig names.

        :param contig_list: List of contig names to extract.
        :type contig_list: list
        :return: Dictionary with contig names as keys and sequences as values.
        :rtype: dict
        """
        contig_dict = {}
        reader = needletail.parse_fastx_file(self.fasta_file)
        for record in reader:
            # pyhmmer uses the first word of the header as the ID, so split on whitespace
            if record.id.strip().split(" ")[0] in contig_list:
                contig_dict[record.id] = record.seq
        return contig_dict

    def write_fasta(self, contig_dict, outfile):
        """
        Write a dictionary of contigs to a fasta file.

        :param contig_dict: Dictionary with contig names as keys and sequences as values.
        :type contig_dict: dict
        :param outfile: Path to the output file.
        :type outfile: str
        :return: None
        """
        with open(outfile, 'w') as out_handle:
            for contig_name, seq in contig_dict.items():
                out_handle.write(f">{contig_name}\n{seq}\n")

    def write_fasta_coords(self, rdrp_coords_list, outfile, seq_type):
        """
        Write a list of RdRp coordinates to a fasta file.

        :param rdrp_coords_list: List of tuples containing contig name and RdRp coordinates.
        :type rdrp_coords_list: list
        :param outfile: Path to the output file.
        :type outfile: str
        :param seq_type: Type of sequence (prot or nuc).
        :type seq_type: str
        :return: None
        """
        if self.logger:
            self.logger.silent_log(f"Processing {len(rdrp_coords_list)} coordinates")
            self.logger.silent_log(f"First few coordinates: {rdrp_coords_list[:3]}")

        reader = needletail.parse_fastx_file(self.fasta_file)
        matches_found = 0
        with open(outfile, 'w') as out_handle:
            for record in reader:
                # pyhmmer uses the first word of the header as the ID, so split on whitespace
                record_id = record.id.strip().split(" ")[0]
                if self.logger:
                    self.logger.silent_log(f"Processing record with ID: '{record_id}'")
                for contig_name, rdrp_from, rdrp_to in rdrp_coords_list:
                    contig_name = str(contig_name).strip()
                    if self.logger:
                        self.logger.silent_log(f"Comparing record '{record_id}' with contig '{contig_name}'")
                    if record_id == contig_name:
                        matches_found += 1
                        seq = record.seq[rdrp_from-1:rdrp_to]
                        fasta_header = f"{record_id}_RdRp_{rdrp_from}-{rdrp_to}"
                        out_handle.write(f">{fasta_header}\n{seq}\n")
                        if self.logger:
                            self.logger.silent_log(f"Match found! Writing sequence of length {len(seq)}")
                    else:
                        if self.logger:
                            self.logger.silent_log(f"No match - lengths: {len(record_id)}|{len(contig_name)}, "
                                                 f"record_id bytes: {record_id.encode()}, contig bytes: {contig_name.encode()}")
        
        if self.logger:
            self.logger.silent_log(f"Total matches found: {matches_found}")


class mmseqs_parser:

    def __init__(self, mmseqs_tax_out_file, mmseqs_s_out_file):
        self.mmseqs_tax_out_file = mmseqs_tax_out_file
        self.mmseqs_s_out_file = mmseqs_s_out_file


    def parse_mmseqs_tax_lca(self):
        """
        Parse the MMseqs2 taxonomy output file.

        :return: Dictionary with contig names as keys and taxonomy lineages as values.
        :rtype: dict
        """
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
        """
        Parse the MMseqs2 easy-search output file.

        :return: Dictionary with contig names as keys and lists of hit information as values.
        :rtype: dict
        """
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
                    tophit_dict[contig] = [target, fident, alnlen, eval, bits, qcov, lineage]
                else:
                    continue

        return tophit_dict

    def tax_to_rdrpcatch(self, rdrpcatch_out, extended_rdrpcatch_out, seq_type):
        """
        Add taxonomy information to the RdRpCATCH output file.

        :param rdrpcatch_out: Path to the RdRpCATCH output file.
        :type rdrpcatch_out: str
        :param extended_rdrpcatch_out: Path to the extended RdRpCATCH output file.
        :type extended_rdrpcatch_out: str
        :param seq_type: Type of sequence (prot or nuc).
        :type seq_type: str
        :return: None
        """
        lca_dict = self.parse_mmseqs_tax_lca()
        tophit_dict = self.parse_mmseqs_e_search_tophit()

        df = pl.read_csv(rdrpcatch_out, separator='\t')

        # drop columns that are not needed
        df = df.drop(["Best_hit_norm_bitscore_profile", "Best_hit_norm_bitscore_contig",
                        "Best_hit_ID_score"])
        
        # Create new columns for taxonomy information
        # For translated sequences, use the frame-specific name
        lookup_col = 'Translated_contig_name (frame)' if seq_type == 'nuc' else 'Contig_name'
        
        df = df.with_columns([
            pl.Series(name='MMseqs_Taxonomy_2bLCA', values=[lca_dict.get(row[lookup_col], '') for row in df.iter_rows(named=True)]),
            pl.Series(name='MMseqs_TopHit_accession', values=[tophit_dict.get(row[lookup_col], ['', '', '', '', '', '', ''])[0] for row in df.iter_rows(named=True)]),
            pl.Series(name='MMseqs_TopHit_fident', values=[tophit_dict.get(row[lookup_col], ['', '', '', '', '', '', ''])[1] for row in df.iter_rows(named=True)]),
            pl.Series(name='MMseqs_TopHit_alnlen', values=[tophit_dict.get(row[lookup_col], ['', '', '', '', '', '', ''])[2] for row in df.iter_rows(named=True)]),
            pl.Series(name='MMseqs_TopHit_eval', values=[tophit_dict.get(row[lookup_col], ['', '', '', '', '', '', ''])[3] for row in df.iter_rows(named=True)]),
            pl.Series(name='MMseqs_TopHit_bitscore', values=[tophit_dict.get(row[lookup_col], ['', '', '', '', '', '', ''])[4] for row in df.iter_rows(named=True)]),
            pl.Series(name='MMseqs_TopHit_qcov', values=[tophit_dict.get(row[lookup_col], ['', '', '', '', '', '', ''])[5] for row in df.iter_rows(named=True)]),
            pl.Series(name='MMseqs_TopHit_lineage', values=[tophit_dict.get(row[lookup_col], ['', '', '', '', '', '', ''])[6] for row in df.iter_rows(named=True)])
        ])

        # Sort by Best_hit_bitscore
        sorted_df = df.sort("Best_hit_bitscore", descending=True)

        sorted_df.write_csv(extended_rdrpcatch_out, separator='\t')


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









