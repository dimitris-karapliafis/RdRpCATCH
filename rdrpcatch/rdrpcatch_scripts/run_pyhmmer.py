import pyhmmer
import os

class pyhmmsearch:

    def __init__(self, hmmsearch_out_path, seq_file, hmm_file, cpus, e, incdomE, domE, incE, z):

        self.hmmsearch_out_path = hmmsearch_out_path
        self.seq_file = seq_file
        self.hmm_file = hmm_file
        self.cpus = cpus
        self.e = e
        self.incdomE = incdomE
        self.domE = domE
        self.incE = incE
        self.z = z



    def run_pyhmmsearch(self):
        """
            TODO: 1. Add option to run hmmsearch on long sequences (longer than 100kb) as pyhmmer.Pipeline is not able to handle
            TODO: long sequences.  See: https://pyhmmer.readthedocs.io/en/latest/api/plan7.html#pyhmmer.plan7.LongTargetsPipeline
            TODO: 2. Parameters are now hardcoded, add option to change them
        """

        if not os.path.exists(self.hmmsearch_out_path):

            with pyhmmer.plan7.HMMPressedFile(self.hmm_file) as handle:
                hmms = list(handle)

            with pyhmmer.easel.SequenceFile(self.seq_file, digital=True) as handle:
                db = list(handle)

            with open(self.hmmsearch_out_path, 'wb') as handle:
                title_line = ["#t_name", "t_acc", "tlen", "q_name", "q_acc", "qlen", "E-value",
                              "score", "bias", "dom_num", "dom_total", "dom_c_value", "dom_i_value", "dom_score",
                              "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
                              "description of target"]
                handle.write("\t".join(title_line).encode("utf-8") + b"\n")

                for result in pyhmmer.hmmer.hmmsearch(hmms,
                                                      db,
                                                      cpus=self.cpus,
                                                      E=self.e,
                                                      incdomE=self.incdomE,
                                                      domE=self.domE,
                                                      incE=self.incE,
                                                      Z=self.z):

                    result.write(handle, format="domains", header=False)

        return self.hmmsearch_out_path

