import os
import subprocess
import sys
from pathlib import Path


class seqkit:

    def __init__(self, input_file: Path ,
                 output_file: Path ,
                 log_file: Path ,
                 threads: int=4):

            self.input_file = input_file
            self.output_file = output_file
            self.log_file = log_file
            self.threads = threads

    def run_seqkit_seq(self, length_thr: int=400):

        seqkit_cmd = ["seqkit",
                      "seq",
                      "--threads",
                      str(self.threads),
                      "-m",
                      str(length_thr),
                      str(self.input_file),
                      "-o",
                      str(self.output_file)]

        with open(self.log_file, 'w') as fout:

            try:
                subprocess.run(seqkit_cmd, stdout=fout, stderr=fout, shell=False, check=True)

            except subprocess.CalledProcessError as e:
                cmd_str = ' '.join(seqkit_cmd)
                raise Exception(f"Error running seqkit command: {cmd_str}")

        return str(self.output_file)


    def run_seqkit_translate(self, gen_code: int=1, frame: int=6):

            seqkit_cmd = ["seqkit",
                        "translate",
                        "--threads",
                        str(self.threads),
                        "--clean",
                        "--append-frame",
                        "-f",
                        f"{frame}",
                        "-T",
                        f"{gen_code}",
                        str(self.input_file),
                        "-o",
                        str(self.output_file)]

            with open(self.log_file, 'w') as fout:

                try:
                    subprocess.run(seqkit_cmd, stdout=fout, stderr=fout, shell=False, check=True)

                except subprocess.CalledProcessError as e:

                    cmd_str = ' '.join(seqkit_cmd)
                    raise Exception(f"Error running seqkit command: {cmd_str}")

            return str(self.output_file)




