"""
Wrapper for the ColabScan package.


"""
import parse_args
import utils
import paths
import run_pyhmmer
import fetch_dbs
import format_pyhmmer_out
import os
from pathlib import Path
import run_seqkit
import plot
import pandas as pd
import warnings
import gui


def main():
    args = parse_args.parse_args()

    commands= {'scan': run_scan,
        'download': run_download,
               "gui": run_gui}

    if args.command in commands:
        commands[args.command](args)
    else:
        print(f"Unknown command: {args.command}")
        exit(1)


def run_download(args):

    if not args.destination_dir:
        fetch_dbs.db_downloader(paths.colabscan_input.db_dir).download_db()
        fetch_dbs.db_downloader(paths.colabscan_input.db_dir).extract_db()
        fetch_dbs.db_downloader(paths.colabscan_input.db_dir).del_tar()
    else:
        fetch_dbs.db_downloader(Path(args.destination_dir)).download_db()
        fetch_dbs.db_downloader(Path(args.destination_dir)).extract_db()
        fetch_dbs.db_downloader(Path(args.destination_dir)).del_tar()


def run_gui(args):

    gui_runner = gui.colabscanner_gui()
    gui_runner.run()


def run_scan(args):


    ## Ignore warnings
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)

    ## Parse arguments
    input_file = args.input
    output_dir = args.output
    db_options = args.db_options
    seq_type = args.seq_type
    verbose = args.verbose

    ## HMMsearch parameters
    e = args.evalue
    cpus = args.cpus
    incdomE = args.incdomEvalue
    domE = args.domEvalue
    incE = args.incEvalue
    z = args.zvalue

    ## seqkit seq parameters
    length_thr = args.length_thr

    ## seqkit translate parameters
    gen_code = args.gen_code
    frame = args.frame


    ## Set output directories
    prefix = Path(input_file).stem
    outputs = paths.colabscan_output(prefix, Path(output_dir))

    ## Set up logger
    log_file = outputs.log_file
    if not os.path.exists(outputs.output_dir):
        os.makedirs(outputs.output_dir)
    logger = utils.Logger(log_file)

    ## Start time
    start_time = logger.start_timer()

    ## Check fasta validity
    if not utils.fasta_checker(input_file).check_fasta_validity():
        raise Exception("Invalid fasta file.")
    else:
        if verbose:
            logger.loud_log(f"Valid fasta file: {input_file}")
        else:
            logger.silent_log(f"Valid fasta file: {input_file}")

    ## Check sequence type
    if not seq_type:
        seq_type = utils.fasta_checker(input_file).check_seq_type()
    if verbose:
        logger.loud_log(f"Sequence type: {seq_type}")
    else:
        logger.silent_log(f"Sequence type: {seq_type}")

    ## Create hmm database directory object
    if args.hmm_dir:
        hmm_db_dir = Path(args.hmm_dir)
    else:
        hmm_db_dir = paths.colabscan_input.hmm_dbs_dir

    if verbose:
        logger.loud_log(f"HMM database directory: {hmm_db_dir}")
    else:
        logger.silent_log(f"HMM database directory: {hmm_db_dir}")


    ## Fetch HMM databases- RVMT, NeoRdRp, NeoRdRp.2.1, TSA_Olendraite, RDRP-scan, Lucaprot
    rvmt_hmm_db = fetch_dbs.db_fetcher(hmm_db_dir, "RVMT").fetch_db_path()
    if verbose:
        logger.loud_log(f"RVMT HMM database fetched from: {rvmt_hmm_db}")
    else:
        logger.silent_log(f"RVMT HMM database fetched from: {rvmt_hmm_db}")
    neordrp_hmm_db = fetch_dbs.db_fetcher(hmm_db_dir, "NeoRdRp").fetch_db_path()
    if verbose:
        logger.loud_log(f"NeoRdRp HMM database fetched from: {neordrp_hmm_db}")
    else:
        logger.silent_log(f"NeoRdRp HMM database fetched from: {neordrp_hmm_db}")
    neordrp_2_hmm_db = fetch_dbs.db_fetcher(hmm_db_dir, "NeoRdRp.2.1").fetch_db_path()
    if verbose:
        logger.loud_log(f"NeoRdRp.2.1 HMM database fetched from: {neordrp_2_hmm_db}")
    else:
        logger.silent_log(f"NeoRdRp.2.1 HMM database fetched from: {neordrp_2_hmm_db}")
    tsa_olen_hmm_db = fetch_dbs.db_fetcher(hmm_db_dir, "TSA_Olendraite").fetch_db_path()
    if verbose:
        logger.loud_log(f"TSA_Olendraite HMM database fetched from: {tsa_olen_hmm_db}")
    else:
        logger.silent_log(f"TSA_Olendraite HMM database fetched from: {tsa_olen_hmm_db}")
    rdrpscan_hmm_db = fetch_dbs.db_fetcher(hmm_db_dir, "RDRP-scan").fetch_db_path()
    if verbose:
        logger.loud_log(f"RDRP-scan HMM database fetched from: {rdrpscan_hmm_db}")
    else:
        logger.silent_log(f"RDRP-scan HMM database fetched from: {rdrpscan_hmm_db}")
    lucaprot_hmm_db = fetch_dbs.db_fetcher(hmm_db_dir, "Lucaprot").fetch_db_path()
    if verbose:
        logger.loud_log(f"Lucaprot HMM database fetched from: {lucaprot_hmm_db}")
    else:
        logger.silent_log(f"Lucaprot HMM database fetched from: {lucaprot_hmm_db}")

    db_name_list = []
    db_path_list = []

    ## Set up HMM databases
    if db_options == ['all']:
        db_name_list = ["RVMT", "NeoRdRp", "NeoRdRp.2.1", "TSA_Olendraite", "RDRP-scan", "Lucaprot"]
        db_path_list = [rvmt_hmm_db, neordrp_hmm_db, neordrp_2_hmm_db, tsa_olen_hmm_db, rdrpscan_hmm_db, lucaprot_hmm_db]

    else:
        for db in db_options:
            if db == "RVMT".lower():
                db_name_list.append("RVMT")
                db_path_list.append(rvmt_hmm_db)
            elif db == "NeoRdRp".lower():
                db_name_list.append("NeoRdRp")
                db_path_list.append(neordrp_hmm_db)
            elif db == "NeoRdRp.2.1":
                db_name_list.append("NeoRdRp.2.1".lower())
                db_path_list.append(neordrp_2_hmm_db)
            elif db == "TSA_Olendraite".lower():
                db_name_list.append("TSA_Olendraite")
                db_path_list.append(tsa_olen_hmm_db)
            elif db == "RDRP-scan".lower():
                db_name_list.append("RDRP-scan")
                db_path_list.append(rdrpscan_hmm_db)
            elif db == "Lucaprot".lower():
                db_name_list.append("Lucaprot")
                db_path_list.append(lucaprot_hmm_db)
            else:
                raise Exception(f"Invalid database option: {db}")

    if not os.path.exists(outputs.hmm_output_dir):
        outputs.hmm_output_dir.mkdir(parents=True)

    if seq_type == 'nuc':
        if verbose:
            logger.loud_log("Nucleotide sequence detected.")
        else:
            logger.silent_log("Nucleotide sequence detected.")

        set_dict = {}
        translated_set_dict = {}
        df_list = []


        ## Filter out sequences with length less than 400 bp with seqkit
        if verbose:
            logger.loud_log("Filtering out sequences with length less than 400 bp.")
        else:
            logger.silent_log("Filtering out sequences with length less than 400 bp.")

        if not os.path.exists(outputs.seqkit_seq_output_dir):
            outputs.seqkit_seq_output_dir.mkdir(parents=True)

        seqkit_out = outputs.seqkit_seq_output_dir / f"{prefix}_filtered.fasta"
        run_seqkit.seqkit(input_file, seqkit_out, log_file, threads=4).run_seqkit_seq(length_thr)
        if verbose:
            logger.loud_log(f"Filtered sequence written to: {seqkit_out}")
        else:
            logger.silent_log(f"Filtered sequence written to: {seqkit_out}")

        ## Translate nucleotide sequences to protein sequences with seqkit
        if verbose:
            logger.loud_log("Translating nucleotide sequences to protein sequences.")
        else:
            logger.silent_log("Translating nucleotide sequences to protein sequences.")

        if not os.path.exists(outputs.seqkit_translate_output_dir):
            outputs.seqkit_translate_output_dir.mkdir(parents=True)
        seqkit_translate_out = outputs.seqkit_translate_output_dir / f"{prefix}_seqkit_translate.fasta"

        seqkit_translate_out = run_seqkit.seqkit(seqkit_out, seqkit_translate_out, log_file, threads=4).run_seqkit_translate(gen_code, frame)

        if verbose:
            logger.loud_log(f"Translated sequence written to: {seqkit_translate_out}")
        else:
            logger.silent_log(f"Translated sequence written to: {seqkit_translate_out}")

        for db_name,db_path in zip(db_name_list, db_path_list):
            hmmsearch_out_path = outputs.hmm_output_dir / f"{prefix}_{db_name}_hmmsearch_out.txt"
            if verbose:
                logger.loud_log(f"HMM output path: {hmmsearch_out_path}")
            else:
                logger.silent_log(f"HMM output path: {hmmsearch_out_path}")

            start_hmmsearch_time = logger.start_timer()
            hmm_out = run_pyhmmer.pyhmmsearch(hmmsearch_out_path, seqkit_translate_out, db_path, cpus, e, incdomE, domE, incE,
                                              z).run_pyhmmsearch()
            end_hmmsearch_time = logger.stop_timer(verbose)
            if verbose:
                logger.loud_log(f"{db_name} HMMsearch Runtime: {end_hmmsearch_time}")
            else:
                logger.silent_log(f"{db_name} HMMsearch Runtime: {end_hmmsearch_time}")

            if verbose:
                logger.loud_log(f"Pyhmmer output written to: {hmm_out}")
            else:
                logger.silent_log(f"Pyhmmer output written to: {hmm_out}")
            if not os.path.exists(outputs.formatted_hmm_output_dir):
                outputs.formatted_hmm_output_dir.mkdir(parents=True)

            form_hmm_out = outputs.formatted_hmm_output_dir / f"{prefix}_{db_name}_hmmsearch_out_formatted.txt"
            format_pyhmmer_out.hmmsearch_formatter(hmm_out, form_hmm_out, seq_type)
            if verbose:
                logger.loud_log(f"Formatted Pyhmmer output written to: {form_hmm_out}")
            else:
                logger.silent_log(f"Formatted Pyhmmer output written to: {form_hmm_out}")
            if not os.path.exists(outputs.lowest_evalue_dir):
                outputs.lowest_evalue_dir.mkdir(parents=True)

            lowest_evalue_out = outputs.lowest_evalue_dir / f"{prefix}_{db_name}_lowest_evalue_hits.txt"
            format_pyhmmer_out.hmmsearch_format_helpers(form_hmm_out, seq_type).lowest_evalue_hits(
                lowest_evalue_out)
            if verbose:
                logger.loud_log(f"Lowest e-value hits written to: {lowest_evalue_out}")
            else:
                logger.silent_log(f"Lowest e-value hits written to: {lowest_evalue_out}")

            set_dict[db_name] = format_pyhmmer_out.hmmsearch_format_helpers(form_hmm_out,
                                                                            seq_type).hmm_to_contig_set()
            translated_set_dict[db_name] = format_pyhmmer_out.hmmsearch_format_helpers(form_hmm_out,
                                                                                       'prot').hmm_to_contig_set()

            # Convert to pandas dataframe, add db_name column and append to df_list
            df = pd.read_csv(lowest_evalue_out, sep='\t')
            df['db_name'] = db_name
            df_list.append(df)


        if not os.path.exists(outputs.plot_outdir):
            outputs.plot_outdir.mkdir(parents=True)

        if not os.path.exists(outputs.tsv_outdir):
            outputs.tsv_outdir.mkdir(parents=True)

        if len(db_name_list) > 1:
            if verbose:
                logger.loud_log("Generating upset plot.")
            else:
                logger.silent_log("Generating upset plot.")

            plot.Plotter(outputs.plot_outdir, outputs.tsv_outdir, prefix).upset_plotter(set_dict)

        # Combine all the dataframes in the list
        combined_df = pd.concat(df_list, ignore_index=True)
        # Write the combined dataframe to a tsv file
        for col in ['E-value', 'score', 'norm_bitscore_profile', 'norm_bitscore_contig',
                    'ID_score', 'profile_coverage', 'contig_coverage']:
            combined_df[col] = pd.to_numeric(combined_df[col], errors='coerce')

        combined_df_path = outputs.tsv_outdir / f"{prefix}_combined_df.tsv"
        combined_df.to_csv(outputs.tsv_outdir / combined_df_path, sep="\t", index=False)

        if verbose:
            logger.loud_log(f"Combined dataframe written to: {outputs.tsv_outdir / f'{prefix}_combined_df.tsv'}")
        else:
            logger.silent_log(f"Combined dataframe written to: {outputs.tsv_outdir / f'{prefix}_combined_df.tsv'}")
        # Generate e-value plot
        plot.Plotter(outputs.plot_outdir, outputs.tsv_outdir, prefix).plot_evalue(combined_df)
        # Generate score plot
        plot.Plotter(outputs.plot_outdir, outputs.tsv_outdir, prefix).plot_score(combined_df)
        # Generate normalized bitscore plot
        plot.Plotter(outputs.plot_outdir, outputs.tsv_outdir, prefix).plot_norm_bitscore_profile(combined_df)
        # Generate normalized bitscore contig plot
        plot.Plotter(outputs.plot_outdir, outputs.tsv_outdir, prefix).plot_norm_bitscore_contig(combined_df)
        # Generate ID score plot
        plot.Plotter(outputs.plot_outdir, outputs.tsv_outdir, prefix).plot_ID_score(combined_df)
        # Generate Profile coverage plot
        plot.Plotter(outputs.plot_outdir, outputs.tsv_outdir, prefix).plot_profile_coverage(combined_df)
        # Generate contig coverage plot
        plot.Plotter(outputs.plot_outdir, outputs.tsv_outdir, prefix).plot_contig_coverage(combined_df)

        # Extract all the contigs
        combined_set = set.union(*[value for value in set_dict.values()])
        translated_combined_set = set.union(*[value for value in translated_set_dict.values()])
        # Write a fasta file with all the contigs
        if not os.path.exists(outputs.fasta_output_dir):
            outputs.fasta_output_dir.mkdir(parents=True)

        fasta_out = outputs.fasta_output_dir / f"{prefix}_contigs.fasta"
        trans_fasta_out = outputs.fasta_output_dir / f"{prefix}_translated_contigs.fasta"
        utils.fasta(input_file).write_fasta(utils.fasta(input_file).extract_contigs(combined_set), fasta_out)
        utils.fasta(seqkit_translate_out).write_fasta(utils.fasta(seqkit_translate_out).extract_contigs(translated_combined_set),trans_fasta_out)
        colabscan_output= outputs.tsv_outdir / f"{prefix}_colabscan_output.tsv"
        format_pyhmmer_out.hmmsearch_output_writter().write_hmmsearch_hits(combined_df_path, seq_type, colabscan_output)
        if verbose:
            logger.loud_log(f"Contigs written to: {fasta_out}")
            logger.loud_log(f"Translated contigs written to: {trans_fasta_out}")
        else:
            logger.silent_log(f"Contigs written to: {fasta_out}")
            logger.silent_log(f"Translated contigs written to: {trans_fasta_out}")


    elif seq_type == 'prot':

        if verbose:
            logger.loud_log("Protein sequence detected.")
        else:
            logger.silent_log("Protein sequence detected.")

        set_dict = {}
        df_list = []

        for db_name,db_path in zip (db_name_list, db_path_list):
            hmmsearch_out_path = outputs.hmm_output_dir / f"{prefix}_{db_name}_hmmsearch_out.txt"

            if verbose:
                logger.loud_log(f"HMM output path: {hmmsearch_out_path}")
            else:
                logger.silent_log(f"HMM output path: {hmmsearch_out_path}")
            hmm_out = run_pyhmmer.pyhmmsearch(hmmsearch_out_path, input_file, db_path, cpus, e, incdomE, domE, incE, z).run_pyhmmsearch()
            if verbose:
                logger.loud_log(f"Pyhmmer output written to: {hmm_out}")
            else:
                logger.silent_log(f"Pyhmmer output written to: {hmm_out}")
            if not os.path.exists(outputs.formatted_hmm_output_dir):
                outputs.formatted_hmm_output_dir.mkdir(parents=True)
            form_hmm_out = outputs.formatted_hmm_output_dir / f"{prefix}_{db_name}_hmmsearch_out_formatted.txt"
            format_pyhmmer_out.hmmsearch_formatter(hmm_out, form_hmm_out, seq_type)
            if verbose:
                logger.loud_log(f"Formatted Pyhmmer output written to: {form_hmm_out}")
            else:
                logger.silent_log(f"Formatted Pyhmmer output written to: {form_hmm_out}")

            # Extract lowest e-value hits from the formatted hmm output

            if not os.path.exists(outputs.lowest_evalue_dir):
                outputs.lowest_evalue_dir.mkdir(parents=True)

            lowest_evalue_out = outputs.lowest_evalue_dir / f"{prefix}_{db_name}_lowest_evalue_hits.txt"
            format_pyhmmer_out.hmmsearch_format_helpers(form_hmm_out,seq_type).lowest_evalue_hits(lowest_evalue_out)
            if verbose:
                logger.loud_log(f"Lowest e-value hits written to: {lowest_evalue_out}")
            else:
                logger.silent_log(f"Lowest e-value hits written to: {lowest_evalue_out}")

            set_dict[db_name] = format_pyhmmer_out.hmmsearch_format_helpers(form_hmm_out,seq_type).hmm_to_contig_set()

            # Convert to pandas dataframe, add db_name column and append to df_list
            df = pd.read_csv(lowest_evalue_out, sep='\t')
            df['db_name'] = db_name
            df_list.append(df)

        if not os.path.exists(outputs.plot_outdir):
            outputs.plot_outdir.mkdir(parents=True)

        if not os.path.exists(outputs.tsv_outdir):
            outputs.tsv_outdir.mkdir(parents=True)

        if len(db_name_list) > 1:
            if verbose:
                logger.loud_log("Generating upset plot.")
            else:
                logger.silent_log("Generating upset plot.")

            plot.Plotter(outputs.plot_outdir,outputs.tsv_outdir, prefix).upset_plotter(set_dict)

        # Combine all the dataframes in the list
        combined_df = pd.concat(df_list, ignore_index=True)
        # Write the combined dataframe to a tsv file
        for col in ['E-value', 'score', 'norm_bitscore_profile', 'norm_bitscore_contig',
                    'ID_score', 'profile_coverage', 'contig_coverage']:
            combined_df[col] = pd.to_numeric(combined_df[col], errors='coerce')

        combined_df_path = outputs.tsv_outdir / f"{prefix}_combined_df.tsv"
        combined_df.to_csv(outputs.tsv_outdir / combined_df_path, sep="\t", index=False)

        if verbose:
            logger.loud_log(f"Combined dataframe written to: {outputs.tsv_outdir / f'{prefix}_combined_df.tsv'}")
        else:
            logger.silent_log(f"Combined dataframe written to: {outputs.tsv_outdir / f'{prefix}_combined_df.tsv'}")

        # Generate e-value plot
        plot.Plotter(outputs.plot_outdir,outputs.tsv_outdir, prefix).plot_evalue(combined_df)
        # Generate score plot
        plot.Plotter(outputs.plot_outdir,outputs.tsv_outdir, prefix).plot_score(combined_df)
        # Generate normalized bitscore plot
        plot.Plotter(outputs.plot_outdir,outputs.tsv_outdir, prefix).plot_norm_bitscore_profile(combined_df)
        # Generate normalized bitscore contig plot
        plot.Plotter(outputs.plot_outdir,outputs.tsv_outdir, prefix).plot_norm_bitscore_contig(combined_df)
        # Generate ID score plot
        plot.Plotter(outputs.plot_outdir,outputs.tsv_outdir, prefix).plot_ID_score(combined_df)
        # Generate Profile coverage plot
        plot.Plotter(outputs.plot_outdir,outputs.tsv_outdir, prefix).plot_profile_coverage(combined_df)
        # Generate contig coverage plot
        plot.Plotter(outputs.plot_outdir,outputs.tsv_outdir, prefix).plot_contig_coverage(combined_df)

        # Extract all the contigs
        combined_set = set.union(*[value for value in set_dict.values()])
        # Write a fasta file with all the contigs
        if not os.path.exists(outputs.fasta_output_dir):
            outputs.fasta_output_dir.mkdir(parents=True)

        fasta_out = outputs.fasta_output_dir / f"{prefix}_contigs.fasta"
        utils.fasta(input_file).write_fasta(utils.fasta(input_file).extract_contigs(combined_set), fasta_out)
        if verbose:
            logger.loud_log(f"Contigs written to: {fasta_out}")
        else:
            logger.silent_log(f"Contigs written to: {fasta_out}")
        colabscan_output= outputs.tsv_outdir / f"{prefix}_colabscan_output.tsv"
        format_pyhmmer_out.hmmsearch_output_writter().write_hmmsearch_hits(combined_df_path, seq_type, colabscan_output)
        if verbose:
            logger.loud_log(f"ColabScan output file written to: {fasta_out}")
        else:
            logger.silent_log(f"ColabScan output file written to: {fasta_out}")



if __name__ == "__main__":
    main()












