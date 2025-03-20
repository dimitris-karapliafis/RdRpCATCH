import os
import pyhmmer 
class pyhmmsearch:

    def __init__(self, hmmsearch_out_path, seq_file, hmm_file, cpus, e, incdomE, domE, incE, z):
        self.hmmsearch_out_path = hmmsearch_out_path
        self.hmmsearch_out_path_custom = str(self.hmmsearch_out_path.with_suffix('.custom.tsv'))
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
        # import pyhmmer

        if not os.path.exists(self.hmmsearch_out_path):

            with pyhmmer.plan7.HMMPressedFile(self.hmm_file) as handle:
                hmms = list(handle)

            with pyhmmer.easel.SequenceFile(self.seq_file, digital=True) as handle:
                db = list(handle)

            with open(self.hmmsearch_out_path, 'wb') as raw_out:
                title_line = ["t_name", "t_acc", "tlen", "q_name", "q_acc", "qlen", "E-value",
                              "score", "bias", "dom_num", "dom_total", "dom_c_value", "dom_i_value", "dom_score",
                              "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
                              "description of target"]
                raw_out.write("\t".join(title_line).encode("utf-8") + b"\n")
                
                with open(self.hmmsearch_out_path_custom, 'wb') as custom_out:
                    custom_out.write("\t".join(title_line).encode("utf-8") + b"\n")
                    for result in pyhmmer.hmmer.hmmsearch(hmms,
                                                      db,
                                                      cpus=self.cpus,
                                                      E=self.e,
                                                      incdomE=self.incdomE,
                                                      domE=self.domE,
                                                      incE=self.incE,
                                                      Z=self.z):

                        result.write(raw_out, format="domains", header=False)
                        if len(result) >= 1:
                            # result.reported.
                        # print(hits.query_name.decode())
                            for hit in result:
                                hit_desc = hit.accession or bytes("", "utf-8")
                                t_desc = hit.description or bytes("-", "utf-8")

                                # print(dir(hit.domains.ex))
                                # hit_name  = hit.name.decode()
                                # join the prot name and acc into a single string because God knows why there are spaces in fasta headers
                                # full_prot_name = f"{hit_name} {hit_desc.decode()}"
                                
                                total_domains = len(hit.domains.included)
                                dom_desc = result.query.description or bytes("", "utf-8")
                                
                                for i, domain in enumerate(hit.domains.included):
                                    domain_num = i + 1
                                    
                                    # print(dir(domain.alignment))
                                    # remove the non-numeric characters from the posterior_probabilities string, then convert to int
                                    import re
                                    # print(domain.alignment.posterior_probabilities)
                                    aligned_probs = (re.sub(r'[^0-9]', '', domain.alignment.posterior_probabilities))
                                    mean_aligned_prob = sum(int(digit) for digit in aligned_probs) / len(domain.alignment.posterior_probabilities)
                                    MEA = mean_aligned_prob 
                                    # print(MEA)
                                    outputline = [
                                        f"{hit.name.decode()}",  # t_name (protein)
                                        f"{hit_desc.decode()}",  # t_acc (empty if none)
                                        f"{hit.length}",  # tlen (protein length)
                                        f"{result.query.name.decode()}", # q_name (HMM name)
                                        f"{dom_desc.decode()}",  # q_acc (empty if none)
                                        f"{domain.alignment.hmm_length}", # qlen (HMM length)
                                        f"{hit.evalue}", # E-value
                                        f"{hit.score}", # score
                                        f"{hit.bias}", # bias
                                        f"{domain_num}", # dom_num (number of this domain)
                                        f"{total_domains}", # dom_total (total number of domains)
                                        f"{domain.c_evalue}", # dom_c_value
                                        f"{domain.i_evalue}", # dom_i_value
                                        f"{domain.score}", # dom_score
                                        f"{domain.bias}", # dom_bias
                                        f"{domain.alignment.hmm_from}", # hmm_from (query from)
                                        f"{domain.alignment.hmm_to}", # hmm_to (query to)
                                        f"{domain.alignment.target_from}", # ali_from (target from)
                                        f"{domain.alignment.target_to}", # ali_to (target to)
                                        f"{domain.env_from}", # env_from
                                        f"{domain.env_to}", # env_to
                                        f"{MEA}", # acc
                                        f"{t_desc.decode()}" # description of target
                                    ]
                                    custom_out.write(("\t".join(outputline) + "\n").encode())
        return self.hmmsearch_out_path

    def run_pyhmmsearch_long_sequences(self):
        """
        Run hmmsearch for sequences longer than 100,000 residues.
        """
        import pyhmmer
        
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

