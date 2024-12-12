import os
import subprocess
import sys
from pathlib import Path


class seqkit:

    def __init__(self, input_file: Path ,
                 output_file: Path ,
                 log_file: Path ,
                 threads: int=4,
                 length_thr: int=400):

            self.input_file = input_file
            self.output_file = output_file
            self.log_file = log_file
            self.threads = threads
            self.length_thr = length_thr



    def run_seqkit(self):

        seqkit_cmd = ["seqkit",
                      "seq",
                      "--threads",
                      str(self.threads),
                      "-m",
                      str(self.length_thr),
                      self.input_file,
                      "-o",
                      self.output_file]

        with open(self.log_file, 'w') as fout:

            try:
                subprocess.run(seqkit_cmd, stdout=fout, stderr=fout, shell=True, check=True)

            except subprocess.CalledProcessError as e:

                cmd_str = ' '.join(seqkit_cmd)
                raise Exception(f"Error running seqkit command: {cmd_str}")




