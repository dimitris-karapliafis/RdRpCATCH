import os
import subprocess
import sys
from pathlib import Path


class transeq:

    def __init__(self, input_file: Path ,
                 output_file: Path ,
                 log_file: Path ,
                 gen_code: int,
                 frame: int):

            self.input_file = input_file
            self.output_file = output_file
            self.log_file = log_file
            self.gen_code = gen_code
            self.frame = frame

    def run_transeq(self):

        transeq_cmd = ["transeq",
                      str(self.input_file),
                      str(self.output_file),
                      f"-table {self.gen_code}",
                      f"-frame {self.frame}",
                      "-clean"]

        with open(self.log_file, 'w') as fout:
            try:
                subprocess.run(transeq_cmd, stdout=fout, stderr=fout, shell=True, check=True)

            except subprocess.CalledProcessError as e:
                cmd_str = ' '.join(transeq_cmd)
                raise Exception(f"Error running transeq command: {cmd_str}")

        return self.output_file


