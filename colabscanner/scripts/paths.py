from dataclasses import dataclass
from pathlib import Path


class classproperty(property):
    def __get__(self, cls, owner):
        return classmethod(self.fget).__get__(None, owner)()

@dataclass
class colabscan_input:

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
class colabscan_output:

    prefix: str
    output_dir: Path
    @property
    def hmm_output_dir (self):
        return self.output_dir / "hmm_output"

    # @property
    # def hmm_output_path(self,db_name, fasta_name):
    #     return self.hmm_output_dir / f"{self.prefix}_{db_name}_{fasta_name}_hmmsearch.txt"
    @property
    def formatted_hmm_output_dir(self):
        return self.output_dir / "formatted_hmm_output"

    @property
    def lowest_evalue_dir(self):
        return self.output_dir / "lowest_evalue_hmm_output"

    @property
    def seqkit_output_dir(self):
        return self.output_dir / "seqkit_output"
    @property
    def seqkit_translate_output_dir(self):
        return self.output_dir / "seqkit_translate_output"
    @property
    def transeq_output_dir (self):
        return self.output_dir / "transeq_output"
    @property
    def plot_outdir(self):
        return self.output_dir / "plots"
    @property
    def tsv_outdir(self):
        return self.output_dir / "tsv_files"

    @property
    def fasta_output_dir(self):
        return self.output_dir / "fasta_output"

    @property
    def log_file(self):
        return self.output_dir / f"{self.prefix}_colabscan.log"
    @property
    def output_file(self):
        return self.output_dir / f"{self.prefix}.out"

