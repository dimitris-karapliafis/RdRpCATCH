from dataclasses import dataclass
from pathlib import Path


class classproperty(property):
    def __get__(self, cls, owner):
        return classmethod(self.fget).__get__(None, owner)()

@dataclass
class rdrpcatch_input:

    #TODO: Change this line for final version

    source_dir : Path = Path(__file__).parents[0].parents[0].parents[0]


    @classproperty
    def db_dir(cls):
        return cls.source_dir / 'DBs'

    @classproperty
    def hmm_dbs_dir(cls):
        return cls.source_dir / 'DBs'/ 'hmm_dbs'


    @classproperty
    def test_dir(cls):
        return cls.source_dir / 'test'


    @classproperty
    def input_fasta(cls):
        return cls.source_dir / "input.fasta"


@dataclass
class rdrpcatch_output:

    prefix: str
    output_dir: Path

    @property
    def tmp_dir(self):
        return self.output_dir / "tmp"

    @property
    def hmm_output_dir (self):
        return self.tmp_dir /"hmm_output"

    def hmm_output_path(self, db_name):
        return self.hmm_output_dir / f"{self.prefix}_{db_name}_hmmsearch_output.txt"

    @property
    def formatted_hmm_output_dir(self):
        return self.tmp_dir / "formatted_hmm_output"

    def formatted_hmm_output_path(self, db_name):
        return self.formatted_hmm_output_dir / f"{self.prefix}_{db_name}_hmm_output_formatted.txt"

    @property
    def best_hit_dir(self):
        return self.tmp_dir / "best_hit_hmm_output"

    def best_hit_path(self, db_name):
        return self.best_hit_dir / f"{self.prefix}_{db_name}_hmm_output_best_hit.txt"

    @property
    def seqkit_seq_output_dir(self):
        return self.tmp_dir/ "seqkit_seq_output"

    @property
    def seqkit_seq_output_path(self):
        return self.seqkit_seq_output_dir / f"{self.prefix}_seqkit_seq_output.fasta"

    @property
    def seqkit_translate_output_dir(self):
        return self.tmp_dir/ "seqkit_translate_output"

    @property
    def seqkit_translate_output_path(self):
        return self.seqkit_translate_output_dir / f"{self.prefix}_seqkit_translate_output.fasta"


    @property
    def tsv_outdir(self):
        return self.tmp_dir/ "tsv_files"
    @property
    def combined_tsv_path(self):
        return self.tsv_outdir / f"{self.prefix}_combined.tsv"

    @property
    def mmseqs_tax_output_dir(self ):
        return self.tmp_dir/ "mmseqs_tax_output"

    @property
    def mmseqs_tax_output_prefix(self):
        return self.mmseqs_tax_output_dir / f"{self.prefix}_mmseqs_tax"

    @property
    def mmseqs_tax_log_path(self):
        return self.log_dir / f"{self.prefix}_mmseqs_tax.log"

    @property
    def mmseqs_tax_output_lca_path(self):
        return self.mmseqs_tax_output_dir / f"{self.prefix}_mmseqs_tax_lca.tsv"
    @property
    def mmseqs_e_search_output_dir(self):
        return self.tmp_dir/ "mmseqs_e_search_output"
    @property
    def mmseqs_e_search_log_path(self):
        return self.log_dir/ f"{self.prefix}_mmseqs_e_search.log"
    @property
    def mmseqs_e_search_output_prefix(self):
        return self.mmseqs_e_search_output_dir / f"{self.prefix}_mmseqs_e_search"

    @property
    def mmseqs_e_search_output_path(self):
        return self.mmseqs_e_search_output_dir / f"{self.prefix}_mmseqs_e_search.tsv"


    @property
    def plot_outdir(self):
        return self.output_dir / f"{self.prefix}_rdrpcatch_plots"

    @property
    def fasta_output_dir(self):
        return self.output_dir / f"{self.prefix}_rdrpcatch_fasta"

    @property
    def fasta_nuc_out_path(self):
        return self.fasta_output_dir / f"{self.prefix}_contigs.fasta"
    @property
    def fasta_prot_out_path(self):
        return self.fasta_output_dir / f"{self.prefix}_transl_contigs.fasta"

    @property
    def rdrpcatch_output(self):
        return self.tsv_outdir / f"{self.prefix}_rdrpcatch_output.tsv"

    @property
    def extended_rdrpcatch_output(self):
        return self.output_dir / f"{self.prefix}_rdrpcatch_output_annotated.tsv"
    @property
    def log_dir(self):
        return self.tmp_dir / f"{self.prefix}_logs"
    @property
    def log_file(self):
        return self.log_dir / f"{self.prefix}_rdrpcatch.log"
    @property
    def gff_output_dir(self):
        return self.output_dir / f"{self.prefix}_gff_files"
    @property
    def gff_output_path(self):
        return self.gff_output_dir / f"{self.prefix}_rdrpcatch.gff"



