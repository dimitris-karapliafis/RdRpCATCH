import argparse




def parse_args():

    parser = argparse.ArgumentParser(description="ColabScan: A package for scanning metagenomic sequences for RdRps.")
    scanner_group= parser.add_argument_group('Scanner options')

    scanner_group.add_argument("-i", "--input",  help="Path to the input FASTA file.")
    scanner_group.add_argument("-o", "--output", help="Path to the output directory.")

    download_group = parser.add_argument_group('Database download options')
    download_group.add_argument('-download','--download_flag', action='store_true',  help="Download HMM databases.")


    return parser.parse_args()



def main():

    args = parse_args()
    input_file = args.input
    output_dir = args.output

    print(f"Input file: {input_file}")
    print(f"Output directory: {output_dir}")

    return input_file, output_dir




