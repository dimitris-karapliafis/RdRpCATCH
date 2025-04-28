"""
Wrapper for the RdRpCATCH package.

"""
import os
from pathlib import Path
from rich.console import Console
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy") # see https://moyix.blogspot.com/2022/09/someones-been-messing-with-my-subnormals.html

def main():
    pass


# def run_gui():
#
#     gui_runner = gui.colabscanner_gui()
#     gui_runner.run()


def bundle_results(output_dir, prefix):
    """
    Bundle the results into a tar.gz file.

    :param output_dir: Path to the output directory.
    :type output_dir: str
    :param prefix: Prefix for the output files.
    :type prefix: str
    :return: Path to the bundled file
    :rtype: str
    """
    import tarfile
    import datetime
    
    # Create timestamp for the archive name
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_name = f"{prefix}_rdrpcatch_results_{timestamp}.tar.gz"
    archive_path = os.path.join(output_dir, archive_name)
    
    # Create tar.gz archive
    with tarfile.open(archive_path, "w:gz") as tar:
        # Add all relevant directories
        for dir_name in [f"{prefix}_rdrpcatch_fasta", f"{prefix}_rdrpcatch_plots", 
                        f"{prefix}_gff_files", "tmp"]:
            dir_path = os.path.join(output_dir, dir_name)
            if os.path.exists(dir_path):
                tar.add(dir_path, arcname=dir_name)
        
        # Add the main output file
        output_file = os.path.join(output_dir, f"{prefix}_rdrpcatch_output_annotated.tsv")
        if os.path.exists(output_file):
            tar.add(output_file, arcname=os.path.basename(output_file))
    
    return archive_path

def run_scan(input_file, output_dir, db_options, db_dir, seq_type, verbose, e,incdomE,domE,incE,z, cpus, length_thr, gen_code, bundle, keep_tmp, overwrite):
    """
    Run RdRpCATCH scan.

    :param input_file: Path to the input FASTA file.
    :type input_file: str
    :param output_dir: Path to the output directory.
    :type output_dir: str
    :param db_options: List of databases to search against.
    :type db_options: list
    :param db_dir: Path to the directory containing RdRpCATCH databases.
    :type db_dir: str
    :param seq_type: Type of sequence (prot or nuc).
    :type seq_type: str
    :param verbose: Whether to print verbose output.
    :type verbose: bool
    :param e: E-value threshold for HMMsearch.
    :type e: float
    :param incdomE: Inclusion domain E-value threshold for HMMsearch.
    :type incdomE: float
    :param domE: Domain E-value threshold for HMMsearch.
    :type domE: float
    :param incE: Inclusion E-value threshold for HMMsearch.
    :type incE: float
    :param z: Number of sequences to search against.
    :type z: int
    :param cpus: Number of CPUs to use for HMMsearch.
    :type cpus: int
    :param length_thr: Minimum length threshold for seqkit seq.
    :type length_thr: int
    :param gen_code: Genetic code to use for translation.
    :type gen_code: int
    :return: None
    """
    from .rdrpcatch_scripts import utils
    from .rdrpcatch_scripts import paths
    from .rdrpcatch_scripts import run_pyhmmer
    from .rdrpcatch_scripts import fetch_dbs
    from .rdrpcatch_scripts import format_pyhmmer_out
    from .rdrpcatch_scripts import run_seqkit
    from .rdrpcatch_scripts import plot
    import polars as pl
    from .rdrpcatch_scripts import mmseqs_tax
    import datetime
    
    ## Ignore warnings
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)

    ## Set output directories
    prefix = Path(input_file).stem
    outputs = paths.rdrpcatch_output(prefix, Path(output_dir))

    ## Set up logger
    log_file = outputs.log_file
    if not os.path.exists(outputs.output_dir):
        os.makedirs(outputs.output_dir)
    elif os.path.exists(outputs.output_dir) and overwrite:
        # If the output directory already exists and force_overwrite is True, remove the existing directory
        import shutil
        shutil.rmtree(outputs.output_dir)
        os.makedirs(outputs.output_dir)
        outputs = paths.rdrpcatch_output(prefix, Path(output_dir))
    else:
        raise FileExistsError(f"Output directory already exists: {outputs.output_dir}, Please choose a different directory"
                              f" or activate the -overwrite flag to overwrite the contents of the directory.")

    if not os.path.exists(outputs.log_dir):
        os.makedirs(outputs.log_dir)

    logger = utils.Logger(log_file)

    logger.silent_log(f"Input File: {input_file}")
    logger.silent_log(f"Output Directory: {output_dir}")
    logger.silent_log(f"Databases: {db_options}")
    logger.silent_log(f"Database Directory: {db_dir}")
    logger.silent_log(f"Sequence Type: {seq_type}")
    logger.silent_log(f"Verbose Mode: {'ON' if verbose else 'OFF'}")
    logger.silent_log(f"E-value: {e}")
    logger.silent_log(f"Inclusion E-value: {incE}")
    logger.silent_log(f"Domain E-value: {domE}")
    logger.silent_log(f"Inclusion Domain E-value: {incdomE}")
    logger.silent_log(f"Z-value: {z}")
    logger.silent_log(f"CPUs: {cpus}")
    logger.silent_log(f"Length Threshold: {length_thr}")
    logger.silent_log(f"Genetic Code: {gen_code}")
    logger.silent_log(f"Bundle Results: {'ON' if bundle else 'OFF'}")
    logger.silent_log(f"Save Temporary Files: {'ON' if keep_tmp else 'OFF'}")

    ## Start time
    start_time = logger.start_timer()

    ## Check fasta validity
    if not utils.fasta_checker(input_file, logger).check_fasta_validity():
        raise Exception("Invalid fasta file.")
    else:
        if verbose:
            logger.loud_log(f"Valid fasta file: {input_file}")
        else:
            logger.silent_log(f"Valid fasta file: {input_file}")

    ## Check sequence type
    if not seq_type:
        seq_type = utils.fasta_checker(input_file, logger).check_seq_type()
    if verbose:
        logger.loud_log(f"Sequence type: {seq_type}")
    else:
        logger.silent_log(f"Sequence type: {seq_type}")

    ## Check sequence length in .fasta files, if >100000, pyHMMER breaks
    if seq_type == 'nuc':
        utils.fasta_checker(input_file, logger).check_seq_length(300000)
    if seq_type == 'prot':
        utils.fasta_checker(input_file, logger).check_seq_length(100000)

    logger.loud_log("Fetching HMM databases...")

    ## Fetch HMM databases- RVMT, NeoRdRp, NeoRdRp.2.1, TSA_Olendraite, RDRP-scan, Lucaprot
    rvmt_hmm_db = fetch_dbs.db_fetcher(db_dir).fetch_hmm_db_path("RVMT")
    if verbose:
        logger.loud_log(f"RVMT HMM database fetched from: {rvmt_hmm_db}")
    else:
        logger.silent_log(f"RVMT HMM database fetched from: {rvmt_hmm_db}")
    neordrp_hmm_db = fetch_dbs.db_fetcher(db_dir).fetch_hmm_db_path("NeoRdRp")
    if verbose:
        logger.loud_log(f"NeoRdRp HMM database fetched from: {neordrp_hmm_db}")
    else:
        logger.silent_log(f"NeoRdRp HMM database fetched from: {neordrp_hmm_db}")
    neordrp_2_hmm_db = fetch_dbs.db_fetcher(db_dir).fetch_hmm_db_path("NeoRdRp.2.1")
    if verbose:
        logger.loud_log(f"NeoRdRp.2.1 HMM database fetched from: {neordrp_2_hmm_db}")
    else:
        logger.silent_log(f"NeoRdRp.2.1 HMM database fetched from: {neordrp_2_hmm_db}")
    tsa_olen_fam_hmm_db = fetch_dbs.db_fetcher(db_dir).fetch_hmm_db_path("TSA_Olendraite_fam")
    if verbose:
        logger.loud_log(f"TSA_Olendraite_fam HMM database fetched from: {tsa_olen_fam_hmm_db}")
    else:
        logger.silent_log(f"TSA_Olendraite_fam HMM database fetched from: {tsa_olen_fam_hmm_db}")

    tsa_olen_gen_hmm_db = fetch_dbs.db_fetcher(db_dir).fetch_hmm_db_path("TSA_Olendraite_gen")
    if verbose:
        logger.loud_log(f"TSA_Olendraite HMM database fetched from: {tsa_olen_gen_hmm_db}")
    else:
        logger.silent_log(f"TSA_Olendraite HMM database fetched from: {tsa_olen_gen_hmm_db}")
    rdrpscan_hmm_db = fetch_dbs.db_fetcher(db_dir).fetch_hmm_db_path("RDRP-scan")
    if verbose:
        logger.loud_log(f"RDRP-scan HMM database fetched from: {rdrpscan_hmm_db}")
    else:
        logger.silent_log(f"RDRP-scan HMM database fetched from: {rdrpscan_hmm_db}")
    lucaprot_hmm_db = fetch_dbs.db_fetcher(db_dir).fetch_hmm_db_path("Lucaprot")
    if verbose:
        logger.loud_log(f"Lucaprot HMM database fetched from: {lucaprot_hmm_db}")
    else:
        logger.silent_log(f"Lucaprot HMM database fetched from: {lucaprot_hmm_db}")

    db_name_list = []
    db_path_list = []

    ## Set up HMM databases
    if db_options == ['all']:
        db_name_list = ["RVMT", "NeoRdRp", "NeoRdRp.2.1", "TSA_Olendraite_fam","TSA_Olendraite_gen", "RDRP-scan", "Lucaprot"]
        db_path_list = [rvmt_hmm_db, neordrp_hmm_db, neordrp_2_hmm_db, tsa_olen_fam_hmm_db,tsa_olen_gen_hmm_db, rdrpscan_hmm_db, lucaprot_hmm_db]

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
            elif db == "TSA_Olendraite_fam".lower():
                db_name_list.append("TSA_Olendraite_fam")
                db_path_list.append(tsa_olen_fam_hmm_db)
            elif db == "TSA_Olendraite_gen".lower():
                db_name_list.append("TSA_Olendraite_gen")
                db_path_list.append(tsa_olen_gen_hmm_db)
            elif db == "RDRP-scan".lower():
                db_name_list.append("RDRP-scan")
                db_path_list.append(rdrpscan_hmm_db)
            elif db == "Lucaprot".lower():
                db_name_list.append("Lucaprot")
                db_path_list.append(lucaprot_hmm_db)
            else:
                raise Exception(f"Invalid database option: {db}")

    # Fetch mmseqs database


    logger.loud_log("Fetching Mmseqs2 databases...")

    mmseqs_db_path = fetch_dbs.db_fetcher(db_dir).fetch_mmseqs_db_path("mmseqs_refseq_riboviria_20250211")

    if verbose:
        logger.loud_log(f"mmseqs database fetched from: {mmseqs_db_path}")
    else:
        logger.silent_log(f"mmseqs database fetched from: {mmseqs_db_path}")

    if not os.path.exists(outputs.hmm_output_dir):
        outputs.hmm_output_dir.mkdir(parents=True)

    if not os.path.exists(outputs.formatted_hmm_output_dir):
        outputs.formatted_hmm_output_dir.mkdir(parents=True)

    if not os.path.exists(outputs.tsv_outdir):
        outputs.tsv_outdir.mkdir(parents=True)

    if not os.path.exists(outputs.plot_outdir):
        outputs.plot_outdir.mkdir(parents=True)

    if not os.path.exists(outputs.tmp_dir):
        outputs.tmp_dir.mkdir(parents=True)

    logger.loud_log("Databases fetched successfully.")

    if seq_type == 'nuc':
        logger.loud_log("Nucleotide sequence detected.")

        set_dict = {}
        translated_set_dict = {}
        df_list = []

        ## Filter out sequences with length less than 400 bp with seqkit
        logger.loud_log("Filtering out sequences with length less than 400 bp.")

        if not os.path.exists(outputs.seqkit_seq_output_dir):
            outputs.seqkit_seq_output_dir.mkdir(parents=True)

        run_seqkit.seqkit(input_file, outputs.seqkit_seq_output_path, log_file, threads=cpus, logger=logger).run_seqkit_seq(length_thr)
        if verbose:
            logger.loud_log(f"Filtered sequence written to: { outputs.seqkit_seq_output_path}")
        else:
            logger.silent_log(f"Filtered sequence written to: { outputs.seqkit_seq_output_path}")

        ## Translate nucleotide sequences to protein sequences with seqkit
        logger.loud_log("Translating nucleotide sequences to protein sequences.")

        if not os.path.exists(outputs.seqkit_translate_output_dir):
            outputs.seqkit_translate_output_dir.mkdir(parents=True)

        run_seqkit.seqkit(outputs.seqkit_seq_output_path, outputs.seqkit_translate_output_path, log_file, threads=cpus, logger=logger).run_seqkit_translate(gen_code, 6)

        if verbose:
            logger.loud_log(f"Translated sequence written to: {outputs.seqkit_translate_output_path}")
        else:
            logger.silent_log(f"Translated sequence written to: {outputs.seqkit_translate_output_path}")

        for db_name,db_path in zip(db_name_list, db_path_list):
            logger.loud_log(f"Running HMMsearch for {db_name} database.")

            if verbose:
                logger.loud_log(f"HMM output path: {outputs.hmm_output_path(db_name)}")
            else:
                logger.silent_log(f"HMM output path: {outputs.hmm_output_path(db_name)}")

            start_hmmsearch_time = logger.start_timer()
            run_pyhmmer.pyhmmsearch(outputs.hmm_output_path(db_name), outputs.seqkit_translate_output_path, db_path, cpus, e, incdomE, domE, incE,
                                              z).run_pyhmmsearch()
            end_hmmsearch_time = logger.stop_timer(start_hmmsearch_time, verbose)
            if verbose:
                logger.loud_log(f"{db_name} HMMsearch Runtime: {end_hmmsearch_time}")
            else:
                logger.silent_log(f"{db_name} HMMsearch Runtime: {end_hmmsearch_time}")

            if verbose:
                logger.loud_log(f"Pyhmmer output written to: {outputs.hmm_output_path(db_name)}")
            else:
                logger.silent_log(f"Pyhmmer output written to: {outputs.hmm_output_path(db_name)}")

            if not os.path.exists(outputs.formatted_hmm_output_dir):
                outputs.formatted_hmm_output_dir.mkdir(parents=True)

            format_pyhmmer_out.hmmsearch_formatter(outputs.hmm_output_path(db_name), outputs.formatted_hmm_output_path(db_name), seq_type)

            if verbose:
                logger.loud_log(f"Formatted Pyhmmer output written to: {outputs.formatted_hmm_output_path(db_name)}")
            else:
                logger.silent_log(f"Formatted Pyhmmer output written to: {outputs.formatted_hmm_output_path(db_name)}")
            if not os.path.exists(outputs.best_hit_dir):
                outputs.best_hit_dir.mkdir(parents=True)

            format_pyhmmer_out.hmmsearch_format_helpers(outputs.formatted_hmm_output_path(db_name), seq_type, logger).highest_bitscore_hits(
                outputs.best_hit_path(db_name))
            if verbose:
                logger.loud_log(f"Highest Bitscore hits written to: {outputs.best_hit_path(db_name)}")
            else:
                logger.silent_log(f"Highest Bitscore hits written to: {outputs.best_hit_path(db_name)}")

            set_dict[db_name] = format_pyhmmer_out.hmmsearch_format_helpers(outputs.formatted_hmm_output_path(db_name),
                                                                            seq_type, logger).hmm_to_contig_set()
            translated_set_dict[db_name] = format_pyhmmer_out.hmmsearch_format_helpers(outputs.formatted_hmm_output_path(db_name),
                                                                                       'prot', logger).hmm_to_contig_set()

            # Convert to dataframe, add db_name column and append to df_list
            df = pl.read_csv(outputs.best_hit_path(db_name), separator='\t')
            df = df.with_columns([
                pl.lit(db_name).alias('db_name')
            ])
            df_list.append(df)

            logger.loud_log(f"HMMsearch for {db_name} completed.")

        logger.loud_log("HMMsearch completed.")

        if not os.path.exists(outputs.plot_outdir):
            outputs.plot_outdir.mkdir(parents=True)

        if not os.path.exists(outputs.tsv_outdir):
            outputs.tsv_outdir.mkdir(parents=True)

        logger.loud_log("Consolidating results.")

        # Combine all the dataframes in the list
        combined_df = pl.concat(df_list, how='vertical_relaxed')
        # Write the combined dataframe to a tsv file
        for col in ['E-value', 'score', 'norm_bitscore_profile', 'norm_bitscore_contig',
                    'ID_score', 'profile_coverage', 'contig_coverage']:
            combined_df = combined_df.with_columns([
                pl.col(col).cast(pl.Float64)
            ])


        combined_df.write_csv(outputs.combined_tsv_path, separator="\t")

        # Check if the combined dataframe is empty
        if combined_df.is_empty():
            logger.loud_log("No hits found by RdRpCATCH. Exiting.")
            return None

        # Generate upset plot
        logger.loud_log("Generating plots.")

        if len(db_name_list) > 1:
            if verbose:
                logger.loud_log("Generating upset plot.")
            else:
                logger.silent_log("Generating upset plot.")

            plot.Plotter(outputs.plot_outdir, outputs.tsv_outdir, prefix).upset_plotter(set_dict)


        if verbose:
            logger.loud_log(f"Combined dataframe written to: {outputs.combined_tsv_path}")
        else:
            logger.silent_log(f"Combined dataframe written to: {outputs.combined_tsv_path}")
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

        logger.loud_log("Extracting RdRp contigs from the input file.")

        # Write a fasta file with all the contigs
        if not os.path.exists(outputs.fasta_output_dir):
            outputs.fasta_output_dir.mkdir(parents=True)

        utils.fasta(input_file).write_fasta(utils.fasta(input_file).extract_contigs(combined_set), outputs.fasta_nuc_out_path)

        utils.fasta(outputs.seqkit_translate_output_path).write_fasta(utils.fasta(outputs.seqkit_translate_output_path).extract_contigs(translated_combined_set),
                                            outputs.fasta_prot_out_path)

        if not os.path.exists(outputs.gff_output_dir):
            outputs.gff_output_dir.mkdir(parents=True)
        hmm_writer = format_pyhmmer_out.hmmsearch_output_writter(logger)
        hmm_writer.write_hmmsearch_hits(outputs.combined_tsv_path, seq_type, outputs.rdrpcatch_output_tsv, outputs.gff_output_path)
        rdrp_coords_list = hmm_writer.get_rdrp_coords(outputs.rdrpcatch_output_tsv,seq_type)
        utils.fasta(outputs.seqkit_translate_output_path, logger).write_fasta_coords(rdrp_coords_list,outputs.fasta_trimmed_out_path, seq_type)

        if verbose:
            logger.loud_log(f"Contigs written to: {outputs.fasta_nuc_out_path}")
            logger.loud_log(f"Translated contigs written to: {outputs.fasta_prot_out_path}")
            logger.loud_log(f"Trimmed contigs written to: {outputs.fasta_trimmed_out_path}")
        else:
            logger.silent_log(f"Contigs written to: {outputs.fasta_nuc_out_path}")
            logger.silent_log(f"Translated contigs written to: {outputs.fasta_prot_out_path}")
            logger.silent_log(f"Trimmed contigs written to: {outputs.fasta_trimmed_out_path}")

        if not os.path.exists(outputs.mmseqs_tax_output_dir):
            outputs.mmseqs_tax_output_dir.mkdir(parents=True)

        logger.loud_log("Running mmseqs easy-taxonomy for taxonomic annotation.")

        mmseqs_tax.mmseqs(outputs.fasta_prot_out_path, mmseqs_db_path, outputs.mmseqs_tax_output_prefix,
                          outputs.mmseqs_tax_output_dir, 7, cpus, outputs.mmseqs_tax_log_path).run_mmseqs_easy_tax_lca()

        logger.loud_log("Running mmseqs easy-search for taxonomic annotation.")

        if not os.path.exists(outputs.mmseqs_e_search_output_dir):
            outputs.mmseqs_e_search_output_dir.mkdir(parents=True)


        mmseqs_tax.mmseqs(outputs.fasta_prot_out_path, mmseqs_db_path, outputs.mmseqs_e_search_output_dir,
                          outputs.mmseqs_e_search_output_path, 7, cpus, outputs.mmseqs_e_search_log_path).run_mmseqs_e_search()

        utils.mmseqs_parser(outputs.mmseqs_tax_output_lca_path, outputs.mmseqs_e_search_output_path).tax_to_rdrpcatch(
            outputs.rdrpcatch_output_tsv, outputs.extended_rdrpcatch_output, seq_type)

        logger.loud_log("Taxonomic annotation completed.")

    elif seq_type == 'prot':

        logger.loud_log("Protein sequence detected.")

        set_dict = {}
        df_list = []

        for db_name,db_path in zip (db_name_list, db_path_list):
            logger.loud_log(f"Running HMMsearch for {db_name} database.")

            if verbose:
                logger.loud_log(f"HMM output path: {outputs.hmm_output_path(db_name)}")
            else:
                logger.silent_log(f"HMM output path: {outputs.hmm_output_path(db_name)}")
            start_hmmsearch_time = logger.start_timer()
            hmm_out = run_pyhmmer.pyhmmsearch(outputs.hmm_output_path(db_name), input_file, db_path, cpus, e, incdomE, domE, incE, z).run_pyhmmsearch()
            end_hmmsearch_time = logger.stop_timer(start_hmmsearch_time,verbose)
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

            format_pyhmmer_out.hmmsearch_formatter(hmm_out, outputs.formatted_hmm_output_path(db_name), seq_type)
            if verbose:
                logger.loud_log(f"Formatted Pyhmmer output written to: {outputs.formatted_hmm_output_path(db_name)}")
            else:
                logger.silent_log(f"Formatted Pyhmmer output written to: {outputs.formatted_hmm_output_path(db_name)}")

            # Extract Highest Bitscore hits from the formatted hmm output

            if not os.path.exists(outputs.best_hit_dir):
                outputs.best_hit_dir.mkdir(parents=True)

            format_pyhmmer_out.hmmsearch_format_helpers(outputs.formatted_hmm_output_path(db_name),seq_type, logger).highest_bitscore_hits(outputs.best_hit_path(db_name))

            if verbose:
                logger.loud_log(f"Highest Bitscore hits written to: {outputs.best_hit_path(db_name)}")
            else:
                logger.silent_log(f"Highest Bitscore hits written to: {outputs.best_hit_path(db_name)}")
            # Here I overwrite prot to nuc, because I need the contig name to extract the contigs
            set_dict[db_name] = format_pyhmmer_out.hmmsearch_format_helpers(outputs.formatted_hmm_output_path(db_name),"nuc", logger).hmm_to_contig_set()

            # Convert to  dataframe, add db_name column and append to df_list
            df = pl.read_csv(outputs.best_hit_path(db_name), separator='\t')
            df = df.with_columns([
                pl.lit(db_name).alias('db_name')
            ])
            df_list.append(df)

            logger.loud_log(f"HMMsearch for {db_name} completed.")

        logger.loud_log("HMMsearch completed.")

        if not os.path.exists(outputs.plot_outdir):
            outputs.plot_outdir.mkdir(parents=True)

        if not os.path.exists(outputs.tsv_outdir):
            outputs.tsv_outdir.mkdir(parents=True)

        logger.loud_log("Consolidating results.")

        # Combine all the dataframes in the list
        combined_df = pl.concat(df_list, how='vertical_relaxed')
        # Write the combined dataframe to a tsv file
        for col in ['E-value', 'score', 'norm_bitscore_profile', 'norm_bitscore_contig',
                    'ID_score', 'profile_coverage', 'contig_coverage']:
            combined_df = combined_df.with_columns([
                pl.col(col).cast(pl.Float64)
            ])

        combined_df.write_csv(outputs.combined_tsv_path, separator="\t")

        # Check if the combined dataframe is empty
        if combined_df.is_empty():
            logger.loud_log("No hits found by RdRpCATCH. Exiting.")
            return None

        # Generate upset plot
        logger.loud_log("Generating plots.")

        if len(db_name_list) > 1:
            if verbose:
                logger.loud_log("Generating upset plot.")
            else:
                logger.silent_log("Generating upset plot.")

            plot.Plotter(outputs.plot_outdir,outputs.tsv_outdir, prefix).upset_plotter(set_dict)


        if verbose:
            logger.loud_log(f"Combined dataframe written to: {outputs.combined_tsv_path}")
        else:
            logger.silent_log(f"Combined dataframe written to: {outputs.combined_tsv_path}")

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

        logger.loud_log("Extracting RdRp contigs from the input file.")

        if not os.path.exists(outputs.fasta_output_dir):
            outputs.fasta_output_dir.mkdir(parents=True)

        utils.fasta(input_file).write_fasta(utils.fasta(input_file).extract_contigs(combined_set), outputs.fasta_prot_out_path)

        if verbose:
            logger.loud_log(f"Full aminoacid contigs written to: {outputs.fasta_prot_out_path}")
        else:
            logger.silent_log(f" Full aminoacid contigs written to: {outputs.fasta_prot_out_path}")

        if not os.path.exists(outputs.gff_output_dir):
            outputs.gff_output_dir.mkdir(parents=True)

        hmm_writer = format_pyhmmer_out.hmmsearch_output_writter(logger)
        hmm_writer.write_hmmsearch_hits(outputs.combined_tsv_path, seq_type, outputs.rdrpcatch_output_tsv, outputs.gff_output_path)
        rdrp_coords_list = hmm_writer.get_rdrp_coords(outputs.rdrpcatch_output_tsv,seq_type)
        utils.fasta(input_file, logger).write_fasta_coords(rdrp_coords_list,outputs.fasta_trimmed_out_path, seq_type)

        if verbose:
            logger.loud_log(f"Trimmed contigs written to: {outputs.fasta_trimmed_out_path}")
        else:
            logger.silent_log(f"Trimmed contigs written to: {outputs.fasta_trimmed_out_path}")

        if not os.path.exists(outputs.mmseqs_tax_output_dir):
            outputs.mmseqs_tax_output_dir.mkdir(parents=True)

        logger.loud_log("Running mmseqs easy-taxonomy for taxonomic annotation.")

        mmseqs_tax.mmseqs(outputs.fasta_prot_out_path, mmseqs_db_path, outputs.mmseqs_tax_output_prefix,
                          outputs.mmseqs_tax_output_dir, 7, cpus, outputs.mmseqs_tax_log_path).run_mmseqs_easy_tax_lca()

        if not os.path.exists(outputs.mmseqs_e_search_output_dir):
            outputs.mmseqs_e_search_output_dir.mkdir(parents=True)

        logger.loud_log("Running mmseqs easy-search for taxonomic annotation.")

        mmseqs_tax.mmseqs(outputs.fasta_prot_out_path, mmseqs_db_path, outputs.mmseqs_e_search_output_dir,
                          outputs.mmseqs_e_search_output_path, 7, cpus, outputs.mmseqs_e_search_log_path).run_mmseqs_e_search()

        utils.mmseqs_parser(outputs.mmseqs_tax_output_lca_path, outputs.mmseqs_e_search_output_path).tax_to_rdrpcatch(
            outputs.rdrpcatch_output_tsv, outputs.extended_rdrpcatch_output, seq_type)






    if not keep_tmp:
        if verbose:
            logger.loud_log("Deleting temporary files.")
        else:
            logger.silent_log("Deleting temporary files.")

        try:
            import shutil
            shutil.rmtree(outputs.tmp_dir)
            logger.silent_log(f"Temporary files deleted.")
        except FileNotFoundError:
            print(f"Directory '{outputs.tmp_dir}' does not exist.")
        except PermissionError:
            print(f"Permission denied while trying to delete '{outputs.tmp_dir}'.")
        except Exception as e:
            print(f"An error occurred: {e}")

    # Bundle results
    if bundle:
        archive_path = bundle_results(output_dir, prefix)
        if verbose:
            logger.loud_log(f"Results bundled into: {archive_path}")
        else:
            logger.silent_log(f"Results bundled into: {archive_path}")

    end_time = logger.stop_timer(start_time, verbose)

    logger.loud_log(f"Total Runtime: {end_time}")

    logger.loud_log("RdRpCATCH completed successfully.")


    return outputs.extended_rdrpcatch_output

if __name__ == "__main__":
    main()
