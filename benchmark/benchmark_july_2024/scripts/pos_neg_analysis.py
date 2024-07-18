"""
Script for running benchmark analysis for RdRp pHMMs on positive and negative ground truth datasets.

********* IMPORTANT *********

Before running the script, make sure to specify the paths to the databases in the script.
Under the section "DATABASES", specify the paths to the following databases:
- hmm_rdrp
- hmm_neordrp
- hmm_olendraite
- hmm_rvmt
- hmm_neordrp_2
Check lines 78-82 for the paths.
The databases are located in the "DBs" directory. The paths should be specified as absolute paths.

Usage:
    python3 pos_neg_analysis.py /path/to/positive_negative_dataset

    Inside the positive_negative_dataset directory, there should be two fasta files with protein sequences:
    - all_dataset_positive.fasta
    - all_dataset_negative.fasta

    The files were retrieved from:
    https://github.com/alibaba/LucaProt

    - all_dataset_positive.fasta:
    Viral RdRP(Positive: 5,979)
    http://47.93.21.181/LucaProt/data/rdrp/all_dataset_positive.fasta.zip

    - all_dataset_negative.fasta:
    Non-viral RdRP(Negative: 229434)
     the non-virus RdRPs contained proteins from Eukaryota DNA dependent RNA polymerase (Eu DdRP, N=1,184),
     Eukaryota RNA dependent RNA polymerase (Eu RdRP, N=2,233),
     Reverse Transcriptase (RT, N=48,490),
     proteins obtained from DNA viruses (N=1,533),
     non-RdRP proteins obtained from RNA viruses (N=1,574),
     and a wide array of cellular proteins from different functional categories (N=174,420)

    http://47.93.21.181/LucaProt/data/rdrp/all_dataset_negative.fasta.zip

    The files were unzipped and put inside the positive_negative_dataset directory.

    The script will run pHMM searches for RdRp pHMMs on the positive and negative datasets and calculate the following metrics:
    - True Positives (TP)
    - False Positives (FP)
    - True Negatives (TN)
    - False Negatives (FN)
    - Sensitivity
    - Specificity
    - Precision
    - Accuracy
    - F1 score
    - False Positive Rate (FPR)

    The script will output a summary file with the metrics for each RdRp pHMM and E-value threshold.
    The script will also output plots for Sensitivity, Specificity, Precision, Accuracy, F1 score and ROC curve.
    The script will output the results in a directory called "results" inside the input directory.
"""

##########        IMPORTS     ##########
from sys import argv
import os
from Bio import SeqIO
import pandas as pd
from matplotlib import pyplot as plt
from rich.console import Console
import pyhmmer

##########        DATABASES     ##########
########## SPECIFY THE PATH TO THE DATABASES BEFORE RUNNING ##########

hmm_rdrp = "/mnt/c/Users/karso/PycharmProjects/ColaB-Scan/DBs/RDRP-scan/RdRp_HMM_profile_CLUSTALO.db"
hmm_neordrp = "/mnt/c/Users/karso/PycharmProjects/ColaB-Scan/DBs/neo-rdrp/NeoRdRp-HMM.v1.1.hmm"
hmm_olendraite = "/mnt/c/Users/karso/PycharmProjects/ColaB-Scan/DBs/olendraite/conc_prof.hmm"
hmm_rvmt = "/mnt/c/Users/karso/PycharmProjects/ColaB-Scan/DBs/rvmt/RVMT.hmm"
hmm_neordrp_2 = "/mnt/c/Users/karso/PycharmProjects/ColaB-Scan/DBs/NeoRdRp.2.1/NeoRdRp.2.1.hmm"

##########        CLASSES     ##########
class hmmsearch_parser:
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

    def __init__(self, hmm_raw, hmm_processed):
        """
        Constructor for the hmmsearch_parser class.

        :param hmm_raw: Path to the raw hmmsearch output file.
        :type hmm_raw: str
        :param hmm_processed: Path to the processed output file.
        :type hmm_processed: str
        """
        self.data = {}
        self.hmm_output_file = hmm_raw
        parsed_data = self.parse_output(self.hmm_output_file)
        parsed_data = self.calculate_norm_bitscore_profile(parsed_data)
        parsed_data = self.calculate_norm_bitscore_contig(parsed_data)
        parsed_data = self.calculate_norm_bitscore_custom(parsed_data)
        self.data = self.calculate_coverage_stats(parsed_data)
        self.export_processed_file(self.data, hmm_processed)

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


    def export_processed_file(self, data, outfile, p_cov_threshold=0):
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

############        FUNCTIONS     ##########
console = Console()


############          MAIN        ##########


def extract_dir_names(in_path):
    """Extracts directory names from a given path

    :param in_path:
    :return: list of directory names
    """
    dir_names = []

    for content in os.listdir(in_path):
        dir_names.append(content)

    return dir_names


def create_outdir(in_path, dir_names):
    """Creates output directories

    :param in_path:
    :param dir_names:
    :return: list of output paths
    """
    parent_dir = os.path.split(in_path)[0]
    out_paths = []

    for filename in dir_names:

        out_path = os.path.join(parent_dir, f"results/{filename}_results")

        if not os.path.exists(out_path):
            os.makedirs(out_path)
        out_paths.append(out_path)

    return out_paths


def run_pyhmmer_search(hmmsearch_out, seq_file, hmm_file, eval):
    """Runs pyhmmer search

    :param hmmsearch_out: str, path to the output file
    :param seq_file: str, path to the sequence file
    :param hmm_file: str, path to the hmm file
    :param eval: float, e-value threshold
    :return: hmmsearch_out, path to the output file
    """
    if not os.path.exists(hmmsearch_out):

        with pyhmmer.plan7.HMMPressedFile(hmm_file) as handle:
            hmms = list(handle)

        with pyhmmer.easel.SequenceFile(seq_file, digital=True) as handle:
            db = list(handle)

        with open(hmmsearch_out, 'wb') as handle:
            title_line = ["#t_name", "t_acc", "tlen", "q_name", "q_acc", "qlen", "E-value",
                          "score", "bias", "dom_num", "dom_total", "dom_c_value", "dom_i_value", "dom_score",
                          "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc",
                          "description of target"]
            handle.write("\t".join(title_line).encode("utf-8") + b"\n")

            for result in pyhmmer.hmmer.hmmsearch(hmms, db, cpus=20, E=eval, incdomE=eval, domE=eval, incE=eval):
                result.write(handle, format="domains", header=False)

    return hmmsearch_out



def parse_to_set_hmm(hmm_fn):
    """Parses hmmsearch output file and returns a set of contigs with hits

    :param hmm_fn: str, path to the hmmsearch output file
    :return: set of contigs with hits
    """

    contigs_set = set()
    with open(hmm_fn, 'r') as in_handle:
        for line in in_handle:
            if line.startswith('Contig') or line.startswith('#'):
                continue
            contig = line.split('\t')[0]
            contigs_set.add(contig)

    return contigs_set


def count_seqs(fasta_file):
    """Counts the number of sequences in a fasta file

    :param fasta_file: str, path to the fasta file
    :return: int, number of sequences
    """

    seq_count = 0
    for seq in SeqIO.parse(fasta_file, 'fasta'):
        seq_count += 1

    return seq_count


def parse_fasta_heads(fasta_file):
    """Parses fasta file and returns a list of sequence headers

    :param fasta_file: str, path to the fasta file
    :return: list of sequence headers
    """

    heads = []
    for seq in SeqIO.parse(fasta_file, 'fasta'):
        heads.append(seq.id)

    return heads


def write_sum_file(overview_dict, general_outdir, out_fn='pos_neg_summary.tsv'):
    """Writes a summary file with the metrics for each RdRp pHMM and E-value threshold

    :param overview_dict: dict, dictionary with the metrics
    :param general_outdir: str, path to the output directory
    :return: str, path to the output file
    """

    out_path =os.path.join(general_outdir, out_fn)

    with open(out_path, 'w') as out_handle:

        out_handle.write(f"eval\thmm_db_name\tTP\tFP\tTN\tFN\tsensitivity\tspecificity\t"
                         f"precision\taccuracy\tf1_score\tfpr\n")

        for eval, hmm_dict in overview_dict.items():

            for hmm_name, values in hmm_dict.items():
                tp = values['tp']
                fp = values['fp']
                tn = values['tn']
                fn = values['fn']
                sensitivity = values['sensitivity']
                specificity = values['specificity']
                precision = values['precision']
                accuracy = values['accuracy']
                f1_score = values['f1_score']
                fpr = values['fpr']
                out_handle.write(f"{eval}\t{hmm_name}\t{tp}\t{fp}\t{tn}\t{fn}\t{sensitivity}\t{specificity}\t"
                                 f"{precision}\t{accuracy}\t{f1_score}\t{fpr}\n")

    return out_path


def calculate_tp_fp(analysis_dict, total_pos,total_negs):
    """Calculates the True Positives (TP) and False Positives (FP) for each RdRp pHMM and E-value threshold

    :param analysis_dict: dict, dictionary with the analysis results
    :param total_pos: list, list of positive contigs
    :param total_negs: list, list of negative contigs
    :return: dict, dictionary with the metrics
    """

    tp = 0
    fp = 0
    overview_dict = {}

    for eval, sample_dict in analysis_dict.items():
        overview_dict[eval] = {}
        for sample, hmm_dict in sample_dict.items():

            if sample == 'all_dataset_positive':
                for hmm_name, hits in hmm_dict.items():
                    if hmm_name not in overview_dict[eval]:
                        overview_dict[eval][hmm_name] = {}

                    for hit in hits:
                        if hit in total_pos:
                            tp += 1
                    overview_dict[eval][hmm_name]['tp'] = tp
                    tp = 0


            elif sample == 'all_dataset_negative':
                for hmm_name, hits in hmm_dict.items():
                    if hmm_name not in overview_dict[eval]:
                        overview_dict[eval][hmm_name] = {}
                    for hit in hits:
                        if hit in total_negs:
                            fp += 1
                    overview_dict[eval][hmm_name]['fp'] = fp
                    fp = 0


    for eval, hmm_dict in overview_dict.items():
        for hmm_name, counts in hmm_dict.items():
            tn = len(total_negs) - overview_dict[eval][hmm_name]['fp']
            fn = len(total_pos) - overview_dict[eval][hmm_name]['tp']
            overview_dict[eval][hmm_name]['tn'] = tn
            overview_dict[eval][hmm_name]['fn'] = fn

            # Calculate sensitivity (recall)
            sensitivity = overview_dict[eval][hmm_name]['tp'] / (overview_dict[eval][hmm_name]['tp']
                                                                 + overview_dict[eval][hmm_name]['fn'])

            # Calculate specificity
            specificity = overview_dict[eval][hmm_name]['tn'] / (overview_dict[eval][hmm_name]['tn']
                                                                 + overview_dict[eval][hmm_name]['fp'])

            # Calculate precision
            precision = overview_dict[eval][hmm_name]['tp'] / (overview_dict[eval][hmm_name]['tp']
                                                               + overview_dict[eval][hmm_name]['fp'])

            # Calculate accuracy
            accuracy = ((overview_dict[eval][hmm_name]['tp'] + overview_dict[eval][hmm_name]['tn'])
                        / (overview_dict[eval][hmm_name]['tp'] + overview_dict[eval][hmm_name]['tn']
                           + overview_dict[eval][hmm_name]['fp'] + overview_dict[eval][hmm_name]['fn']))

            # Calculate F1 score
            f1_score = 2 * (precision * sensitivity) / (precision + sensitivity)

            # Calculate False Positive Rate
            fpr = overview_dict[eval][hmm_name]['fp'] / (overview_dict[eval][hmm_name]['fp']
                                                         + overview_dict[eval][hmm_name]['tn'])

            overview_dict[eval][hmm_name]['sensitivity'] = sensitivity
            overview_dict[eval][hmm_name]['specificity'] = specificity
            overview_dict[eval][hmm_name]['precision'] = precision
            overview_dict[eval][hmm_name]['accuracy'] = accuracy
            overview_dict[eval][hmm_name]['f1_score'] = f1_score
            overview_dict[eval][hmm_name]['fpr'] = fpr

    return overview_dict


def calculate_and_plot_roc(data):
    """Calculates and plots the Receiver Operating Characteristic (ROC) curve

    :param data: dict, dictionary with the metrics
    :return: None
    """

    fpr = data['fp'] / (data['fp'] + data['tn'])
    tpr = data['tp'] / (data['tp'] + data['fn'])

    plt.figure()
    plt.plot(fpr, tpr, label='ROC curve')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic')
    plt.legend(loc="lower right")
    plt.savefig('roc.png')

# Usage



def plotter(summ_path, general_outdir):
    """Plots the Sensitivity, Specificity, Precision, Accuracy, F1 score and ROC curve

    :param summ_path: str, path to the summary file
    :param general_outdir: str, path to the output directory
    :return: None
    """
    
    plt.rcParams['xtick.labelsize']=14
    plt.rcParams['ytick.labelsize']=14
    # backend to Agg to avoid display
    plt.switch_backend('Agg')

    df = pd.read_csv(summ_path, sep='\t')
    pivot_df = df.pivot(index='eval', columns='hmm_db_name', values='sensitivity')
    plt.figure(figsize=(10, 6))

    # Loop through each hmm_db_name and plot its data
    for hmm_db_name, data in pivot_df.items():
        plt.semilogx(data.index, data.values, marker='o', label=hmm_db_name)  # Use semilogx for logarithmic x-axis
    plt.xlabel('E-values',fontsize=15)
    plt.ylabel('Sensitivity',fontsize=15)
    plt.title('Sensitivity',fontsize=22,weight='bold')
    plt.grid(True)
    plt.savefig(os.path.join(general_outdir, "sensitivity.png"))


    pivot_df = df.pivot(index='eval', columns='hmm_db_name', values='specificity')
    plt.figure(figsize=(10, 6))

    # Loop through each hmm_db_name and plot its data
    for hmm_db_name, data in pivot_df.items():
        plt.semilogx(data.index, data.values, marker='o', label=hmm_db_name)
    plt.xlabel('E-values',fontsize=15)
    plt.ylabel('Specificity',fontsize=15)
    plt.title('Specificity',fontsize=22,weight='bold')
    plt.grid(True)
    plt.savefig(os.path.join(general_outdir, "specificity.png"))

    pivot_df = df.pivot(index='eval', columns='hmm_db_name', values='precision')
    plt.figure(figsize=(10, 6))

    # Loop through each hmm_db_name and plot its data
    for hmm_db_name, data in pivot_df.items():
        plt.semilogx(data.index, data.values, marker='o', label=hmm_db_name)

    plt.xlabel('E-values',fontsize=15)
    plt.ylabel('Precision',fontsize=15)
    plt.title('Precision',fontsize=22,weight='bold')
    plt.legend(fontsize=18)
    plt.grid(True)
    plt.savefig(os.path.join(general_outdir, "precision.png"))

    pivot_df = df.pivot(index='eval', columns='hmm_db_name', values='accuracy')
    plt.figure(figsize=(10, 6))

    # Loop through each hmm_db_name and plot its data
    for hmm_db_name, data in pivot_df.items():
        plt.semilogx(data.index, data.values, marker='o', label=hmm_db_name)

    plt.xlabel('E-values',fontsize=15)
    plt.ylabel('Accuracy',fontsize=15)
    plt.title('Accuracy',fontsize=22,weight='bold')
    plt.legend(fontsize=18)
    plt.grid(True)
    plt.savefig(os.path.join(general_outdir, "accuracy.png"))

    pivot_df = df.pivot(index='eval', columns='hmm_db_name', values='f1_score')
    plt.figure(figsize=(10, 6))

    # Loop through each hmm_db_name and plot its data
    for hmm_db_name, data in pivot_df.items():
        plt.semilogx(data.index, data.values, marker='o', label=hmm_db_name)

    plt.xlabel('E-values',fontsize=15)
    plt.ylabel('F1 score',fontsize=15)
    plt.title('F1 score',fontsize=22,weight='bold')
    plt.legend(fontsize=18)
    plt.grid(True)
    plt.savefig(os.path.join(general_outdir, "f1_score.png"))

    pivot_df_fpr = df.pivot(index='eval', columns='hmm_db_name', values="fpr")
    pivot_df_sens = df.pivot(index='eval', columns='hmm_db_name', values="sensitivity")
    # Create ROC plot
    plt.figure(figsize=(10, 6))

    # Loop through each hmm_db_name and plot its data
    for hmm_db_name in pivot_df_fpr.keys():
        plt.plot(pivot_df_fpr[hmm_db_name], pivot_df_sens[hmm_db_name], marker='o', label=hmm_db_name, markersize=5)

    plt.xlabel('False Positive Rate',fontsize=15)
    plt.ylabel('True Positive Rate',fontsize=15)
    plt.title('ROC',fontsize=22,weight='bold')
    plt.legend(fontsize=18)
    plt.grid(True)
    # plt.plot([0, 1], [0, 1], 'k--')
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.05])
    # plt.legend(loc="lower right")
    plt.savefig(os.path.join(general_outdir, "roc.png"))
    plt.close('all')

def parse_hmmsearch(hmmsearch_out, eval, parsed_outfn):
    """Parses hmmsearch output file and writes the results to a new file

    :param hmmsearch_out: str, path to the hmmsearch output file
    :param eval: float, e-value threshold
    :param parsed_outfn: str, path to the output file
    :return: None
    """

    hmmsearch_processed_fn = os.path.split(hmmsearch_out)[1].split('.')[0] + '_processed.txt'

    hmmsearch_processed_path = os.path.join(os.path.split(hmmsearch_out)[0], hmmsearch_processed_fn)

    if not os.path.exists(hmmsearch_processed_fn):
        hmmsearch_parser(hmmsearch_out, hmmsearch_processed_path)

    with open(hmmsearch_processed_path, 'r') as in_handle, open(parsed_outfn, 'w') as out_handle:
        for line in in_handle:
            if line.startswith('#'):
                out_handle.write(line)
                continue

            line = line.strip().split('\t')
            new_line = []
            for element in line:
                if element != '' or element != ' ':
                    new_line.append(element)
            evalue = float(new_line[6])
            if evalue <= eval:
                out_handle.write("\t".join(new_line) + '\n')




def pos_neg_analysis(in_path):
    """ Wrapper function for running analysis

    :param in_path: str, path to the input directory
    :return: None
    """

    analysis_dict = {}
    evals = [1, 1e-1,1e-2,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15]

    samples = extract_dir_names(in_path)
    console.log(f"Running positive and negative control analysis for samples: {' '.join(samples)}")
    out_dirs = create_outdir(in_path, samples)
    console.log(f"Output directories created: {' '.join(out_dirs)}")
    general_outdir = os.path.split(out_dirs[0])[0]

    for sample, out_dir in zip(samples,out_dirs):
        hmm_outdir = os.path.join(out_dir,"hmm")
        if not os.path.exists(hmm_outdir):
            os.mkdir(hmm_outdir)
        # Define input directory and contig path
        in_dir = os.path.join(in_path, sample)
        input = os.listdir(in_dir)[0]
        data_path = os.path.join(in_dir, input)

        console.log(f"Running analysis for sample {sample} | rdrpscan")
        out_rdrpscan = os.path.join(hmm_outdir, f"{sample}_full_rdrpscan_hmmsearch_out.txt")
        if not os.path.exists(out_rdrpscan):
            run_pyhmmer_search(out_rdrpscan, data_path, hmm_rdrp, 1)

        console.log(f"Running analysis for sample {sample} | olendraite")
        out_olendraite = os.path.join(hmm_outdir, f"{sample}_full_olendraite_hmmsearch_out.txt")
        if not os.path.exists(out_olendraite):
            run_pyhmmer_search(out_olendraite, data_path, hmm_olendraite, 1)

        console.log(f"Running analysis for sample {sample} | neordrp")
        out_neordrp = os.path.join(hmm_outdir, f"{sample}_full_neordrp_hmmsearch_out.txt")
        if not os.path.exists(out_neordrp):
            run_pyhmmer_search(out_neordrp, data_path, hmm_neordrp, 1)

        console.log(f"Running analysis for sample {sample} | rvmt")
        out_rvmt = os.path.join(hmm_outdir, f"{sample}_full_rvmt_hmmsearch_out.txt")
        if not os.path.exists(out_rvmt):
            run_pyhmmer_search(out_rvmt, data_path, hmm_rvmt, 1)

        console.log(f"Running analysis for sample {sample} | neordrp_2")
        out_neordrp_2 = os.path.join(hmm_outdir, f"{sample}_full_neordrp_2_hmmsearch_out.txt")
        if not os.path.exists(out_neordrp_2):
            run_pyhmmer_search(out_neordrp_2, data_path, hmm_neordrp_2, 1)




    for eval in evals:

        console.log(f"Running analysis for e-value {eval}")

        analysis_dict[eval] = {}

        for sample, out_dir in zip(samples, out_dirs):
            # Create hmm output directory
            hmm_outdir = os.path.join(out_dir, 'hmm')
            if not os.path.exists(hmm_outdir):
                os.mkdir(hmm_outdir)

            analysis_dict[eval][sample] = {}

            # Define input directory and contig path
            in_dir = os.path.join(in_path, sample)
            input = os.listdir(in_dir)[0]
            data_path = os.path.join(in_dir, input)

            if input == 'all_dataset_negative.fasta':
                tn = count_seqs(data_path)
                total_negs= parse_fasta_heads(data_path)

            elif input == 'all_dataset_positive.fasta':
                tp = count_seqs(data_path)
                total_pos = parse_fasta_heads(data_path)

            else:
                raise ValueError(f"Input file {input} does not match any of the expected input files")

            console.log(f"Running analysis for sample {sample} | rdrpscan")
            out_rdrpscan = os.path.join(hmm_outdir, f"{sample}_{eval}_rdrpscan_hmmsearch_out.txt")
            in_rdrpscan = os.path.join(hmm_outdir,f"{sample}_full_rdrpscan_hmmsearch_out.txt" )
            if not os.path.exists(out_rdrpscan):
                parse_hmmsearch(in_rdrpscan,eval,out_rdrpscan)
            rdrpscan_set = parse_to_set_hmm(out_rdrpscan)
            analysis_dict[eval][sample]['rdrpscan'] = rdrpscan_set

            console.log(f"Running analysis for sample {sample} | olendraite")
            out_olendraite = os.path.join(hmm_outdir, f"{sample}_{eval}_olendraite_hmmsearch_out.txt")
            in_olendraite = os.path.join(hmm_outdir, f"{sample}_full_olendraite_hmmsearch_out.txt" )
            if not os.path.exists(out_olendraite):
                parse_hmmsearch(in_olendraite, eval, out_olendraite)
            olendraite_set = parse_to_set_hmm(out_olendraite)
            analysis_dict[eval][sample]['tsa_olendraite'] = olendraite_set

            console.log(f"Running analysis for sample {sample} | neordrp")
            out_neordrp= os.path.join(hmm_outdir, f"{sample}_{eval}_neordrp_hmmsearch_out.txt")
            in_neordrp = os.path.join(hmm_outdir, f"{sample}_full_neordrp_hmmsearch_out.txt")
            if not os.path.exists(out_neordrp):
                parse_hmmsearch(in_neordrp,eval,out_neordrp)
            neordrp_set = parse_to_set_hmm(out_neordrp)
            analysis_dict[eval][sample]['neordrp'] = neordrp_set

            console.log(f"Running analysis for sample {sample} | rvmt")
            out_rvmt = os.path.join(hmm_outdir, f"{sample}_{eval}_rvmt_hmmsearch_out.txt")
            in_rvmt = os.path.join(hmm_outdir,f"{sample}_full_rvmt_hmmsearch_out.txt")

            if not os.path.exists(out_rvmt):
                parse_hmmsearch(in_rvmt,eval,out_rvmt)
            rvmt_set = parse_to_set_hmm(out_rvmt)
            analysis_dict[eval][sample]['rvmt'] = rvmt_set

            console.log(f"Running analysis for sample {sample} | neordrp_2")
            out_neordrp_2 = os.path.join(hmm_outdir, f"{sample}_{eval}_neordrp_2_hmmsearch_out.txt")
            in_neordrp_2 = os.path.join(hmm_outdir,f"{sample}_full_neordrp_2_hmmsearch_out.txt")

            if not os.path.exists(out_neordrp_2):
                parse_hmmsearch(in_neordrp_2,eval,out_neordrp_2)
            neordrp_2_set = parse_to_set_hmm(out_neordrp_2)
            analysis_dict[eval][sample]['neordrp_2'] = neordrp_2_set

    print("total_pos", len(total_pos))
    print("total_negs", len(total_negs))
    overview_dict = calculate_tp_fp(analysis_dict, total_pos, total_negs)
    sum_path = write_sum_file(overview_dict, general_outdir)
    plotter(sum_path, general_outdir)


if __name__ == "__main__":
    in_path = argv[1]
    pos_neg_analysis(in_path)


