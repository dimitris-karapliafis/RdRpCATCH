import subprocess
import os



class mmseqs:

    def __init__(self,fasta_fn, mmseqs_db, mmseqs_out_prefix,outdir_path, sens, cpus, log_file):
        self.fasta_fn = fasta_fn
        self.mmseqs_db = mmseqs_db
        self.out_prefix = mmseqs_out_prefix
        self.outdir_path = outdir_path
        self.sens = sens
        self.cpus = cpus
        self.log_file = log_file

    def run_mmseqs_easy_tax_lca(self):
        """Run mmseqs easy-tax command."""

        mmseqs_easy_tax_cmd = ["mmseqs",
                               "easy-taxonomy",
                               str(self.fasta_fn),
                               str(self.mmseqs_db),
                               str(self.out_prefix),
                               f"{str(self.outdir_path)}/tmp",
                               "--tax-lineage",
                               "1",
                               "--alignment-mode",
                               "3",
                               "--start-sens",
                               str(self.sens),
                               "--threads",
                               str(self.cpus)
                               ]

        try:
            with open(self.log_file, 'w') as fout:
                subprocess.run(mmseqs_easy_tax_cmd, stdout=fout, stderr=fout, shell=False, check=True)

        except subprocess.CalledProcessError as e:
            cmd_str = ' '.join(mmseqs_easy_tax_cmd)
            raise Exception(f"Error running mmseqs easy-tax command: {cmd_str}")


    def run_mmseqs_easy_tax_tophit(self):


        mmseqs_easy_tax_cmd = ["mmseqs",
                               "easy-taxonomy",
                                 self.fasta_fn,
                                    self.mmseqs_db,
                                    self.out_prefix,
                                    "tmp",
                                    "--tax-lineage",
                                    "1",
                                    self.sens,
                                    "--threads",
                                    str(self.cpus),
                                    "--lca-mode",
                                    4]
        try:
            with open(self.log_file, 'w') as fout:
                subprocess.run(mmseqs_easy_tax_cmd, stdout=fout, stderr=fout, shell=False, check=True)

        except subprocess.CalledProcessError as e:
            cmd_str = ' '.join(mmseqs_easy_tax_cmd)
            raise Exception(f"Error running mmseqs easy-tax command: {cmd_str}")


    def run_mmseqs_e_search(self):
        mmseqs_e_search_cmd = ["mmseqs",
                               "easy-search",
                               str(self.fasta_fn),
                               str(self.mmseqs_db),
                               str(self.outdir_path),
                               f"{str(self.out_prefix)}/tmp",
                               "--start-sens",
                               str(self.sens),
                               "--threads",
                               str(self.cpus),
                               "--format-output",
                               "query,target,fident,alnlen,mismatch,gapopen,"
                               "qstart,qend,tstart,tend,evalue,bits,qcov,tcov,taxlineage",
                               "--sort-results",
                               "1"]
        try:
            with open(self.log_file, 'w') as fout:
                subprocess.run(mmseqs_e_search_cmd, stdout=fout, stderr=fout, shell=False, check=True)

        except subprocess.CalledProcessError as e:
            cmd_str = ' '.join(mmseqs_e_search_cmd)
            raise Exception(f"Error running mmseqs easy-search command: {cmd_str}")








