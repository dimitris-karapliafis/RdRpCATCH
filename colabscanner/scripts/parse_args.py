import argparse


def parse_comma_separated_options(options_str):

    allowed_choices = ['RVMT', 'NeoRdRp', 'NeoRdRp.2.1', 'TSA_Olendraite', 'RDRP-scan', 'Lucaprot', 'all']
    lower_choices = [choice.lower() for choice in allowed_choices]

    options = options_str.split(',')
    lower_options = [option.lower() for option in options]
    for option in options:
        if option.lower() not in lower_choices:
            raise argparse.ArgumentTypeError(f"Invalid choice: '{option}' (choose from {', '.join(allowed_choices)})")
    return lower_options



def parse_args():
    parser = argparse.ArgumentParser(description="ColabScan: A package for scanning metagenomic sequences for RdRps.")
    scanner_group = parser.add_argument_group('Scanner options')

    scanner_group.add_argument("-i", "--input", help="Path to the input FASTA file.")
    scanner_group.add_argument("-o", "--output", help="Path to the output directory.")
    scanner_group.add_argument("-dbs", "--db_options", type=parse_comma_separated_options,
                               default= 'all',
                               help="Comma-separated list of databases to search against. Choose from: RVMT,"
                                    " NeoRdRp,"
                                    " NeoRdRp.2.1,"
                                    " TSA_Olendraite,"
                                    " RDRP-scan,"
                                    " Lucaprot,"
                                    " all")

    download_group = parser.add_argument_group('Database download options')
    download_group.add_argument('-download', '--download_flag', action='store_true', help="Download HMM databases.")

    hmmsearch_group = parser.add_argument_group('HMMsearch options')
    hmmsearch_group.add_argument('-e', '--evalue', type=float, default=1e-5, help="E-value threshold for HMMsearch.")
    hmmsearch_group.add_argument('-incE', '--incEvalue', type=float, default=1e-5, help="Inclusion E-value threshold for HMMsearch.")
    hmmsearch_group.add_argument('-domE', '--domEvalue', type=float, default=1e-5, help="Domain E-value threshold for HMMsearch.")
    hmmsearch_group.add_argument('-incdomE', '--incdomEvalue', type=float, default=1e-5, help="Inclusion domain E-value threshold for HMMsearch.")
    hmmsearch_group.add_argument('-z', '--zvalue', type=int, default=1, help="Number of sequences to search against.")
    hmmsearch_group.add_argument('-cpus', '--cpus', type=int, default=1, help="Number of CPUs to use for HMMsearch.")

    transeq_group = parser.add_argument_group('Transeq options')
    transeq_group.add_argument('-gen_code', '--gen_code', type=int, default=1, help="Genetic code to use for translation.")
    transeq_group.add_argument('-frame', '--frame', type=int, default=6, help="Frame to use for translation.")

    return parser.parse_args()




def main():

    args = parse_args()
    input_file = args.input
    output_dir = args.output

    print(f"Input file: {input_file}")
    print(f"Output directory: {output_dir}")

    return input_file, output_dir




