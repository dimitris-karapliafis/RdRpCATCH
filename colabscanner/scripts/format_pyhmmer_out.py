import pandas as pd

class hmmsearch_formatter:
    """
    Class for parsing hmmsearch output files.

    Attributes:
        data (dict): A dictionary containing the parsed data from the hmmscan output file.
        hmm_output_file (str): Path to the hmmscan output file.

    Methods:
        parse_output(hmm_output_file): Parses the hmmsearch output file and returns a dictionary.
        calculate_coverage(data): Calculates the coverage of all domains in a profile.
        get_contig(contig_name): Returns all profiles and domains for a given contig.
        export_processed_file(data, outfile, p_cov_threshold=0): Exports the processed hmmscan output file.
    """

    def __init__(self, hmm_raw, hmm_processed, seq_type):
        """
        Constructor for the hmmsearch_parser class.

        :param hmm_raw: Path to the raw hmmsearch output file.
        :type hmm_raw: str
        :param hmm_processed: Path to the processed output file.
        :type hmm_processed: str

        If PROTEIN: contig name is the first column
        If DNA: contig name is the last column, first column is the translated sequence name (e.g. contig_name_frame)

        """
        self.data = {}
        self.hmm_output_file = hmm_raw
        parsed_data = self.parse_output(self.hmm_output_file)
        parsed_data = self.calculate_norm_bitscore_profile(parsed_data)
        parsed_data = self.calculate_norm_bitscore_contig(parsed_data)
        parsed_data = self.calculate_norm_bitscore_custom(parsed_data)
        self.data = self.calculate_coverage(parsed_data)

        if seq_type == 'PROTEIN':
            self.export_processed_file_aa(self.data, hmm_processed)

        elif seq_type == 'DNA':
            self.export_processed_file_dna(self.data, hmm_processed)

    def parse_output(self, hmm_raw_out):
        """
        Parse hmmsearch output file and return a dictionary with the following structure:
        {contig_name: {profile_name: [[list of domain data]]}}

       :param hmm_output_file: Path to the hmmscan output file.
       :type hmm_output_file: str
       :return: Dictionary with parsed data.
       :rtype: dict
       """

        with open(hmm_raw_out, 'r') as hmm_file:

            for line in hmm_file:
                if line.startswith('#') or line.startswith('t_name'):
                    continue

                tmp_line = line.strip().split()
                for element in tmp_line:
                    element.replace(' ', '')

                desc = tmp_line[22:]
                tmp_line = tmp_line[:22] + [' '.join(desc)]

                if tmp_line[0] not in self.data:
                    self.data[tmp_line[0]] = {tmp_line[3]: [[col for col in tmp_line[1:]]]}
                else:
                    if tmp_line[3] not in self.data[tmp_line[0]]:
                        self.data[tmp_line[0]][tmp_line[3]] = [[col for col in tmp_line[1:]]]
                    else:
                        self.data[tmp_line[0]][tmp_line[3]].append([col for col in tmp_line[1:]])

        return self.data

    def calculate_norm_bitscore_profile(self, data):
        """
        Calculates the BitScore/Length for all domains in a profile.
        :param data:
        :return:
        """

        for contig, profiles in data.items():
                for profile, domains in profiles.items():
                    for domain in domains:
                        model_length = float(domain[4])
                        domain_bitscore = float(domain[6])
                        norm_bitscore = domain_bitscore/model_length
                        domain.append(norm_bitscore)

        return data

    def calculate_norm_bitscore_contig(self, data):
        """
        Calculates the BitScore/Length for the contig size.
        :param data:
        :return:
        """

        for contig, profiles in data.items():
                for profile, domains in profiles.items():
                    for domain in domains:
                        contig_length = float(domain[1])
                        domain_bitscore = float(domain[6])
                        norm_bitscore = domain_bitscore/contig_length
                        domain.append(norm_bitscore)

        return data

    def calculate_norm_bitscore_custom(self, data):
        """
        Calculates the 2*BitScore/Length of contig + Length of profile for all domains in a profile.
        :param data:
        :return:
        """

        for contig, profiles in data.items():
                for profile, domains in profiles.items():
                    for domain in domains:
                        model_length = float(domain[4])
                        contig_length = float(domain[1])
                        domain_bitscore = float(domain[6])
                        norm_bitscore = (domain_bitscore/model_length) + (domain_bitscore/contig_length)
                        domain.append(norm_bitscore)

        return data

    def calculate_coverage(self, data):
        """
        Calculates the % coverage of all domains in a profile.

        :param data: Dictionary with parsed data.
        :type data: dict
        :return: Dictionary with parsed data and domain coverage.
        :rtype: dict
        """
        # TODO: This needs to be thoroughly tested, the %perc coverage is known to be difficult to calculate due to
        # TODO: many diffent scenarios of overlap between domains. See this issue for more info:
        # TODO: COntig coverage takes into account the hmm from-to positions, not the ali_from-to positions, so it can be
        # TODO: greater than 1. See this issue for more info:
        # https://github.com/althonos/pyhmmer/issues/27
        overlap_set = set()
        cov_set = set()
        for contig, profiles in data.items():

            for profile, domains in profiles.items():


                query_length = int(domains[0][4])
                contig_length = int(domains[0][1])
                bitscore = float(domains[0][6])
                for i in range(len(domains)):  # iterate over domains
                    prev_dom_start = int(domains[i][14])
                    prev_dom_end = int(domains[i][15])
                    cov_range = list(range(prev_dom_start, prev_dom_end + 1))
                    for pos in cov_range:
                        cov_set.add(pos)

                    if i == len(domains)-1:
                        domain_coverage = len(cov_set) / query_length
                        contig_coverage = len(cov_set) / contig_length
                        domains[i].append(bitscore/len(cov_set))
                        domains[i].append(len(cov_set))
                        domains[i].append(domain_coverage)
                        domains[i].append(contig_coverage)

                    cov_set = set()

        return data  # return data with domain coverage added

    def get_contig(self, contig_name):
        """
        Returns all profiles and domains for a given contig.

        :param contig_name: Name of the contig.
        :type contig_name: str
        :return: Dictionary with profiles and domains.
        :rtype: dict
        """
        return self.data.get(contig_name, {})



    def export_processed_file_aa(self, data, outfile, p_cov_threshold=0):
        """
        Exports the processed hmmsearch output file.

        :param data: Dictionary with parsed data.
        :type data: dict
        :param outfile: Path to the output file.
        :type outfile: str
        :param p_cov_threshold: Minimum profile coverage threshold, defaults to 0.
        :type p_cov_threshold: int, optional
        :return: Path to the output file.
        :rtype: str
        """
        title_line = ["#t_name", "t_acc", "tlen", "q_name", "q_acc", "qlen", "E-value",
                      "score", "bias", "dom_num", "dom_total", "dom_c_value", "dom_i_value", "dom_score",
                      "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
                      "description of target", "norm_bitscore_profile", "norm_bitscore_contig", 'norm_bitscore_custom',
                      "ID_score",'aln_length','profile_coverage', "contig_coverage"]

        line_list = []
        with open(outfile, 'w') as out:
            for contig, profiles in data.items():
                for profile, domains in profiles.items():
                    if len(profiles[profile]) > 1:
                        line_list.append([contig] + domains[-1])
                    else:
                        for domain in domains:
                            line_list.append([contig] + domain)
            out.write("\t".join(title_line) + '\n')

            for line in line_list:
                if line[-1] > p_cov_threshold:
                    j_line = "\t".join(str(l) for l in line)
                    out.write(j_line + '\n')
        return outfile


    def export_processed_file_dna(self, data, outfile, p_cov_threshold=0):
        """
        Exports the processed hmmsearch output file.

        :param data: Dictionary with parsed data.
        :type data: dict
        :param outfile: Path to the output file.
        :type outfile: str
        :param p_cov_threshold: Minimum profile coverage threshold, defaults to 0.
        :type p_cov_threshold: int, optional
        :return: Path to the output file.
        :rtype: str
        """
        title_line = ["#t_name", "t_acc", "tlen", "q_name", "q_acc", "qlen", "E-value",
                      "score", "bias", "dom_num", "dom_total", "dom_c_value", "dom_i_value", "dom_score",
                      "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
                      "description of target", "norm_bitscore_profile", "norm_bitscore_contig", 'norm_bitscore_custom',
                      "ID_score",'aln_length','profile_coverage', "contig_coverage", 'contig_name']

        line_list = []
        with open(outfile, 'w') as out:
            for contig, profiles in data.items():
                for profile, domains in profiles.items():
                    if len(profiles[profile]) > 1:
                        line_list.append([contig] + domains[-1] + [contig.rsplit('_', 1)[0]])

                    else:
                        for domain in domains:
                            line_list.append([contig] + domain + [contig.rsplit('_', 1)[0]])
            out.write("\t".join(title_line) + '\n')

            for line in line_list:
                if line[-1] > p_cov_threshold:
                    j_line = "\t".join(str(l) for l in line)
                    out.write(j_line + '\n')
        return outfile


class hmmsearch_format_helpers:

    def __init__(self, hmm_outfn, seq_type):

        self.hmm_outfn = hmm_outfn
        self.seq_type = seq_type

    def hmm_to_contig_set(self):
        """
        Returns a set of all contig names in the data.

        :return: Set of contig names.
        :rtype: set
        """

        contig_set = set()
        with open(self.hmm_outfn) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    if self.seq_type == 'PROTEIN':
                        contig_set.add(line.split()[0])
                    elif self.seq_type == 'DNA':
                        contig_set.add(line.split()[-1])

        return contig_set


    def lowest_evalue_hits(self,filtered_file):
        """
        Filters the hmmsearch output file based on the lowest E-value for each contig.

        :param evalue_threshold: E-value threshold for filtering.
        :type evalue_threshold: float
        :return: Path to the filtered output file.
        :rtype: str

        If PROTEIN: contig name is the first column
        If DNA: contig name is the last column, first column is the translated sequence name (e.g. contig_name_frame)
        """
        hmm_dict = {}

        with open(self.hmm_outfn) as in_handle:
            for line in in_handle:
                if line.startswith('#'):
                    title_line = line
                    continue
                line = line.strip().split('\t')
                if self.seq_type == 'PROTEIN':
                    contig_name = line[0]
                elif self.seq_type == 'DNA':
                    contig_name = line[-1]
                if contig_name not in hmm_dict:
                    hmm_dict[contig_name] = [line]
                else:
                    hmm_dict[contig_name].append(line)

        with open(filtered_file, 'w') as out_handle:
            out_handle.write(title_line)
            for contig, hits in hmm_dict.items():
                best_hit = min(hits, key=lambda x: float(x[6]))
                out_handle.write('\t'.join(best_hit) + '\n')


    def extract_col(self, index):
        """Tranforms the hmmsearch output file to a pandas dataframe. Then extracts the column of interest based on index. Outputs a list of the column values.

        :return:
        """
        df = pd.read_csv(self.hmm_outfn, sep='\t')

        return df.iloc[:, index].tolist()



class hmmsearch_combiner:

    def __init__(self, hmmsearch_files, combined_file):
        """
        Constructor for the hmmsearch_combiner class.

        :param hmmsearch_files: List of paths to the hmmsearch output files.
        :type hmmsearch_files: list
        :param combined_file: Path to the combined output file.
        :type combined_file: str
        """
        self.hmmsearch_files = hmmsearch_files
        self.combined_file = combined_file
        self.combine_files(self.hmmsearch_files, self.combined_file)

    def combine_files(self, hmmsearch_files, combined_file):
        """
        Combines multiple hmmsearch output files into a single file.

        :param hmmsearch_files: List of paths to the hmmsearch output files.
        :type hmmsearch_files: list
        :param combined_file: Path to the combined output file.
        :type combined_file: str
        :return: Path to the combined output file.
        :rtype: str
        """
        with open(combined_file, 'w') as out:
            for file in hmmsearch_files:
                with open(file, 'r') as f:
                    for line in f:
                        out.write(line)
        return combined_file