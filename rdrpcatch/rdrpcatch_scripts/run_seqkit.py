class seqkit:

    def __init__(self, input_file,
                 output_file,
                 log_file,
                 threads=4,
                 logger=None):

            self.input_file = input_file
            self.output_file = output_file
            self.log_file = log_file
            self.threads = threads
            self.logger = logger

    def run_seqkit_seq(self, length_thr=400):
        import os
        import subprocess
        import sys
        from pathlib import Path

        if self.logger:
            self.logger.silent_log(f"Running seqkit seq on {self.input_file}")
            self.logger.silent_log(f"Length threshold: {length_thr}")

        seqkit_cmd = ["seqkit",
                      "seq",
                      "--threads",
                      str(self.threads),
                      "-m",
                      str(length_thr),
                      str(self.input_file),
                      "-o",
                      str(self.output_file)]

        if self.logger:
            self.logger.silent_log(f"Running command: {' '.join(seqkit_cmd)}")

        with open(self.log_file, 'w') as fout:
            try:
                subprocess.run(seqkit_cmd, stdout=fout, stderr=fout, shell=False, check=True)
                if self.logger:
                    self.logger.silent_log(f"Successfully filtered sequences to {self.output_file}")

            except subprocess.CalledProcessError as e:
                cmd_str = ' '.join(seqkit_cmd)
                error_msg = f"Error running seqkit command: {cmd_str}"
                if self.logger:
                    self.logger.silent_log(error_msg)
                raise Exception(error_msg)

        return str(self.output_file)


    def run_seqkit_translate(self, gen_code=1, frame=6):
            import os
            import subprocess
            import sys
            from pathlib import Path

            if self.logger:
                self.logger.silent_log(f"Running seqkit translate on {self.input_file}")
                self.logger.silent_log(f"Output will be written to {self.output_file}")
                self.logger.silent_log(f"Using genetic code {gen_code} and frame {frame}")

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

            if self.logger:
                self.logger.silent_log(f"Running command: {' '.join(seqkit_cmd)}")

            with open(self.log_file, 'w') as fout:
                try:
                    subprocess.run(seqkit_cmd, stdout=fout, stderr=fout, shell=False, check=True)
                    # Check the output file exists and has content
                    if os.path.exists(self.output_file):
                        with open(self.output_file, 'r') as f:
                            first_few_lines = [next(f) for _ in range(6)]
                        if self.logger:
                            self.logger.silent_log("First few lines of output:")
                            for line in first_few_lines:
                                self.logger.silent_log(f"{line.strip()}")
                    else:
                        error_msg = f"Output file {self.output_file} was not created!"
                        if self.logger:
                            self.logger.silent_log(error_msg)
                        raise Exception(error_msg)

                except subprocess.CalledProcessError as e:
                    cmd_str = ' '.join(seqkit_cmd)
                    error_msg = f"Error running seqkit command: {cmd_str}"
                    error_details = f"Error details: {str(e)}"
                    if self.logger:
                        self.logger.silent_log(error_msg)
                        self.logger.silent_log(error_details)
                    raise Exception(error_msg)

            return str(self.output_file)




