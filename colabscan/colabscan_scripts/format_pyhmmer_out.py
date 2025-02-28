import pandas as pd
import re

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
        self.data = self.calculate_coverage_stats(parsed_data)

        if seq_type == 'prot':
            self.export_processed_file_aa(self.data, hmm_processed)

        elif seq_type == 'nuc':
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
                        domain.append(round(norm_bitscore,5))

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
                        domain.append(round(norm_bitscore,5))

        return data


    def calculate_coverage_stats(self, data):
        """
        Calculates the % coverage of all domains in a profile.

        :param data: Dictionary with parsed data.
        :type data: dict
        :return: Dictionary with parsed data and domain coverage.
        :rtype: dict
        """
        # TODO: This needs to be thoroughly tested, the %perc coverage is known to be difficult to calculate due to
        # TODO: many diffent scenarios of overlap between domains. See this issue for more info:
        # TODO: Contig coverage takes into account the hmm from-to positions, not the ali_from-to positions, so it can be
        # TODO: greater than 1. See this issue for more info:
        # https://github.com/althonos/pyhmmer/issues/27
        overlap_set = set()
        hmm_cov_set = set()
        env_seq_cov_set = set()
        ali_seq_cov_set = set()
        for contig, profiles in data.items():

            for profile, domains in profiles.items():


                query_length = int(domains[0][4])
                contig_length = int(domains[0][1])
                bitscore = float(domains[0][6])
                for i in range(len(domains)):  # iterate over domains
                    hmm_prev_dom_start = int(domains[i][14])
                    hmm_prev_dom_end = int(domains[i][15])
                    hmm_cov_range = list(range(hmm_prev_dom_start, hmm_prev_dom_end + 1))

                    env_seq_prev_dom_start = int(domains[i][18])
                    env_seq_prev_dom_end = int(domains[i][19])
                    seq_cov_range = list(range(env_seq_prev_dom_start, env_seq_prev_dom_end + 1))

                    ali_seq_prev_dom_start = int(domains[i][16])
                    ali_seq_prev_dom_end = int(domains[i][17])
                    ali_seq_cov_range = list(range(ali_seq_prev_dom_start, ali_seq_prev_dom_end + 1))



                    for pos in hmm_cov_range:
                        hmm_cov_set.add(pos)

                    for pos in seq_cov_range:
                        env_seq_cov_set.add(pos)

                    for pos in ali_seq_cov_range:
                        ali_seq_cov_set.add(pos)

                    if i == len(domains)-1:
                        domain_coverage = len(hmm_cov_set) / query_length
                        contig_coverage = len(env_seq_cov_set) / contig_length
                        id_score = bitscore/len(ali_seq_cov_set)
                        domains[i].append(round(id_score,5))
                        domains[i].append(len(ali_seq_cov_set))
                        domains[i].append(round(domain_coverage,5))
                        domains[i].append(round(contig_coverage,5))
                        domains[i].append(min(env_seq_cov_set))
                        domains[i].append(max(env_seq_cov_set))

                hmm_cov_set = set()
                env_seq_cov_set = set()
                ali_seq_cov_set = set()

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
                      "description of target", "norm_bitscore_profile", "norm_bitscore_contig",
                      "ID_score",'aln_length','profile_coverage', "contig_coverage", "RdRp_start", "RdRp_end"]

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
                      "description of target", "norm_bitscore_profile", "norm_bitscore_contig",
                      "ID_score",'aln_length','profile_coverage', "contig_coverage", "RdRp_start", "RdRp_end", 'contig_name']

        line_list = []
        with open(outfile, 'w') as out:
            for contig, profiles in data.items():
                for profile, domains in profiles.items():
                    if len(profiles[profile]) > 1:
                        line_list.append([contig] + domains[-1] + [re.sub(r'_frame=[+-]?\d+', '', contig)])

                    else:
                        for domain in domains:
                            line_list.append([contig] + domain + [re.sub(r'_frame=[+-]?\d+', '', contig)])
            out.write("\t".join(title_line) + '\n')

            for line in line_list:
                if line[-2] > p_cov_threshold:
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
                    if self.seq_type == 'prot':
                        contig_set.add(line.split()[0])
                    elif self.seq_type == 'nuc':
                        contig_set.add(line.split()[-1])

        return contig_set


    def highest_bitscore_hits(self, filtered_file):
        """
        Filters the hmmsearch output file based on the highest bitscore for each contig.

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
                    title_line = line.strip().split('\t')
                    title_line = title_line[:27] + ['Total_positive_profiles'] + title_line[27:]
                    title_line = '\t'.join(title_line) + '\n'
                    continue
                line = line.strip().split('\t')
                if self.seq_type == 'prot':
                    contig_name = line[0]
                elif self.seq_type == 'nuc':
                    contig_name = line[-1]
                if contig_name not in hmm_dict:
                    hmm_dict[contig_name] = [line]
                else:
                    hmm_dict[contig_name].append(line)

        with open(filtered_file, 'w') as out_handle:
            out_handle.write(title_line)
            for contig, hits in hmm_dict.items():
                total_hits = len(hits)
                best_hit = max(hits, key=lambda x: float(x[7]))
                best_hit = best_hit[:27] + [str(total_hits)] + best_hit[27:]
                out_handle.write('\t'.join(best_hit) + '\n')

    def highest_norm_bit_prof_hits(self, filtered_file):
        """
        Filters the hmmsearch output file based on the highest normalized bitscore for each contig.

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
                if self.seq_type == 'prot':
                    contig_name = line[0]
                elif self.seq_type == 'nuc':
                    contig_name = line[-1]
                if contig_name not in hmm_dict:
                    hmm_dict[contig_name] = [line]
                else:
                    hmm_dict[contig_name].append(line)

        with open(filtered_file, 'w') as out_handle:
            out_handle.write(title_line)
            for contig, hits in hmm_dict.items():
                best_hit = max(hits, key=lambda x: float(x[23]))
                out_handle.write('\t'.join(best_hit) + '\n')


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
                if self.seq_type == 'prot':
                    contig_name = line[0]
                elif self.seq_type == 'nuc':
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
        """Tranforms the hmmsearch output file to a pandas dataframe.
         Then extracts the column of interest based on index. Outputs a list of the column values.

        :return:
        """
        df = pd.read_csv(self.hmm_outfn, sep='\t')

        return df.iloc[:, index].tolist()


class hmmsearch_output_writter:

    def __init__(self):
        """
        Constructor for the hmmsearch_output_writter class.

        :
        """
    def get_hmmsearch_hits(self, hmmsearch_combined_fn, seq_type):

        hmm_dict = {}

        with open(hmmsearch_combined_fn) as in_handle:
            for line in in_handle:
                if line.startswith('#'):
                    title_line = line
                    continue
                line = line.strip().split('\t')
                if seq_type == 'prot':
                    contig_name = line[0]
                elif seq_type == 'nuc':
                    contig_name = line[-2]
                if contig_name not in hmm_dict:
                    hmm_dict[contig_name] = [line]
                else:
                    hmm_dict[contig_name].append(line)

        return hmm_dict

    def write_hmmsearch_hits(self, hmmsearch_combined_fn, seq_type, out_fn,gff_out_fn):

        hmm_dict = self.get_hmmsearch_hits(hmmsearch_combined_fn, seq_type)
        db_list = []
        output_list = []

        for contig, hits in hmm_dict.items():
            best_hit = max(hits, key=lambda x: float(x[7]))
            total_hits = ""


            for hit in hits:
                db_list.append(hit[-1])
                total_hits += f"{hit[-1]}={hit[27]};"

                if seq_type == 'nuc':
                    translated_seq_name = hit[0]
                else:
                    translated_seq_name = "-"

            if seq_type == 'nuc':
                hit_line = [contig, translated_seq_name, best_hit[-1], str(total_hits), best_hit[2], best_hit[3],
                            best_hit[5], best_hit[6], best_hit[7], best_hit[-4], best_hit[-3], best_hit[22],
                            best_hit[23], best_hit[24], best_hit[25], best_hit[26], best_hit[28], best_hit[29]]
            else:
                hit_line = [contig, translated_seq_name, best_hit[-1], str(total_hits), best_hit[2], best_hit[3],
                            best_hit[5], best_hit[6], best_hit[7], best_hit[-3], best_hit[-2], best_hit[22],
                            best_hit[23], best_hit[24], best_hit[25], best_hit[26], best_hit[28], best_hit[29]]

            output_list.append(hit_line)
            db_list = []


        title_line = ["#Contig_name","Translated_contig_name (frame)", "Best_hit_Database", "Total_databases_that_the_contig_was_detected(No_of_Profiles)", "Sequence_length(AA)", "Best_hit_profile_name", "Best_hit_profile_length",
                      "Best_hit_e-value", "Best_hit_bitscore", "RdRp_from(AA)", "RdRp_to(AA)", "Best_hit_description", "Best_hit_norm_bitscore_profile",
                      "Best_hit_norm_bitscore_contig", "Best_hit_ID_score", "Best_hit_aln_length","Best_hit_profile_coverage", "Best_hit_contig_coverage"]

        with open(out_fn, 'w') as out_handle:
            out_handle.write("\t".join(title_line) + '\n')
            for line in output_list:
                out_handle.write('\t'.join(line) + '\n')

        self.write_gff_file(output_list, seq_type, gff_out_fn)

        return out_fn

    def write_gff_file(self, output_list, seq_type, gff_fn):

        with open(gff_fn, 'w') as out_handle:
            out_handle.write("##gff-version 3\n")

            for line in output_list:
                if seq_type == 'nuc':
                    contig_name = line[0]
                    translated_seq_name = line[1]
                    db_name = line[2]
                    db_list = line[3]
                    seq_length = line[4]
                    profile_name = line[5]
                    profile_length = line[6]
                    e_value = line[7]
                    bitscore = line[8]
                    rdrp_start = line[9]
                    rdrp_end = line[10]
                    description = line[11]
                    norm_bitscore_profile = line[12]
                    norm_bitscore_contig = line[13]
                    ID_score = line[14]
                    aln_length = line[15]
                    profile_coverage = line[16]
                    contig_coverage = line[17]

                    out_handle.write(f"{translated_seq_name}\tColabScan\tRdRp_domain\t{rdrp_start}\t{rdrp_end}"
                                     f"\t{bitscore}\t+\t.\tID={contig_name};Profile_Name={profile_name};"
                                     f" Profile_Db={db_name};E-value={e_value};Bitscore={bitscore};"
                                     f"Norm_bitscore_profile={norm_bitscore_profile};"
                                     f"Norm_bitscore_contig={norm_bitscore_contig};"
                                     f"ID_score={ID_score};Aln_length={aln_length};Profile_coverage={profile_coverage};"
                                     f"Contig_coverage={contig_coverage};Description={description}\n")
                else:
                    contig_name = line[0]
                    translated_seq_name = line[1]
                    db_name = line[2]
                    db_list = line[3]
                    seq_length = line[4]
                    profile_name = line[5]
                    profile_length = line[6]
                    e_value = line[7]
                    bitscore = line[8]
                    rdrp_start = line[9]
                    rdrp_end = line[10]
                    description = line[11]
                    norm_bitscore_profile = line[12]
                    norm_bitscore_contig = line[13]
                    ID_score = line[14]
                    aln_length = line[15]
                    profile_coverage = line[16]
                    contig_coverage = line[17]
                    out_handle.write(f"{contig_name}\tColabScan\tRdRp_domain\t{rdrp_start}\t{rdrp_end}\t{bitscore}\t."
                                     f"\t.\tID={contig_name};Profile_Name={profile_name}; Profile_Db={db_name};"
                                     f"E-value={e_value};Bitscore={bitscore};"
                                     f"Norm_bitscore_profile={norm_bitscore_profile};"
                                     f"Norm_bitscore_contig={norm_bitscore_contig};ID_score={ID_score};"
                                     f"Aln_length={aln_length};Profile_coverage={profile_coverage};"
                                     f"Contig_coverage={contig_coverage};Description={description}\n")


        return gff_fn

    def get_rdrp_coords(self, out_fn):

        rdrp_coords = []

        with open(out_fn) as in_handle:
            for line in in_handle:
                if line.startswith('Contig_name'):
                    line = line.strip().split('\t')

                    continue
                line = line.strip().split('\t')

                contig_name = line[0]
                translated_seq_name = line[1]
                rdrp_coords.append([contig_name, translated_seq_name,  line[9], line[10]])
        return rdrp_coords

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

