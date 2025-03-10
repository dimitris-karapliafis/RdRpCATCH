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

    def run_pyhmmsearch_long_sequences(self):
        """
        Run hmmsearch for sequences longer than 100,000 residues.
        """
        if not os.path.exists(self.hmmsearch_out_path):
            with pyhmmer.plan7.HMMPressedFile(self.hmm_file) as handle:
                hmms = list(handle)

            with pyhmmer.easel.SequenceFile(self.seq_file, digital=True) as handle:
                db = list(handle)

            # Create a LongTargetsPipeline instance
            alphabet = pyhmmer.easel.Alphabet.amino()
            pipeline = pyhmmer.plan7.LongTargetsPipeline(alphabet,
                                                         block_length=262144,  # Default block length
                                                         F1=0.02, F2=0.003, F3=3e-05)

            with open(self.hmmsearch_out_path, 'wb') as handle:
                title_line = ["#t_name", "t_acc", "tlen", "q_name", "q_acc", "qlen", "E-value",
                              "score", "bias", "dom_num", "dom_total", "dom_c_value", "dom_i_value", "dom_score",
                              "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
                              "description of target"]
                handle.write("\t".join(title_line).encode("utf-8") + b"\n")

                for hmm in hmms:
                    iterator = pipeline.iterate_seq(hmm, db)
                    max_iterations = 10  # Prevent infinite loop
                    for n in range(max_iterations):
                        _, hits, _, converged, _ = next(iterator)
                        if converged:
                            break

                    # Process hits and write to file
                    for hit in hits:
                        # Assuming hit is a plan7.Hit object
                        # Extract relevant information and write to file
                        # Note: This part might need adjustment based on actual hit structure
                        handle.write(f"{hit.target_name}\t{hit.target_accession}\t{hit.target_length}\t"
                                     f"{hit.query_name}\t{hit.query_accession}\t{hit.query_length}\t"
                                     f"{hit.evalue}\t{hit.score}\t{hit.bias}\t"
                                     f"{hit.domain_num}\t{hit.domain_total}\t{hit.domain_cvalue}\t"
                                     f"{hit.domain_ivalue}\t{hit.domain_score}\t{hit.domain_bias}\t"
                                     f"{hit.hmm_from}\t{hit.hmm_to}\t{hit.ali_from}\t{hit.ali_to}\t"
                                     f"{hit.env_from}\t{hit.env_to}\t{hit.acc}\t"
                                     f"{hit.description}\n".encode("utf-8"))

        return self.hmmsearch_out_path

