"""
Wrapper for the RdRpCATCH package.

"""
from .rdrpcatch_scripts import utils
from .rdrpcatch_scripts import paths
from .rdrpcatch_scripts import run_pyhmmer
from .rdrpcatch_scripts import fetch_dbs
from .rdrpcatch_scripts import format_pyhmmer_out
import os
from pathlib import Path
from .rdrpcatch_scripts import run_seqkit
from .rdrpcatch_scripts import plot
import pandas as pd
import warnings
# from .rdrpcatch_scripts import gui
from .rdrpcatch_scripts import mmseqs_tax
from rich.console import Console

def main():
    pass

def run_download(destination_dir):

    console = Console()

    fetch_dbs.db_downloader(Path(destination_dir)).download_db()
    fetch_dbs.db_downloader(Path(destination_dir)).extract_db()
    fetch_dbs.db_downloader(Path(destination_dir)).del_tar()
    console.log(f"RdRp databases downloaded and extracted successfully. Databases are stored in:\n"
                f"{destination_dir}/DBs ")
    


# def run_gui():
#
#     gui_runner = gui.colabscanner_gui()
#     gui_runner.run()


def run_scan(input_file, output_dir, db_options, db_dir, seq_type, verbose, e,incdomE,domE,incE,z, cpus, length_thr, gen_code):


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

    ## Check sequence length in .fasta files, if >100000, pyHMMER breaks
    if seq_type == 'nuc':
        utils.fasta_checker(input_file).check_seq_length(300000)
    if seq_type == 'prot':
        utils.fasta_checker(input_file).check_seq_length(100000)




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

    if verbose:
        logger.loud_log("Fetching mmseqs databases.")
    else:
        logger.silent_log("Fetching mmseqs databases.")
    mmseqs_db_path = fetch_dbs.db_fetcher(db_dir).fetch_mmseqs_db_path("mmseqs_refseq_riboviria_20250211")

    if verbose:
        logger.loud_log(f"mmseqs database fetched from: {mmseqs_db_path}")
    else:
        logger.silent_log(f"mmseqs database fetched from: {mmseqs_db_path}")

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


        run_seqkit.seqkit(input_file,  outputs.seqkit_seq_output_path, log_file, threads=cpus).run_seqkit_seq(length_thr)
        if verbose:
            logger.loud_log(f"Filtered sequence written to: { outputs.seqkit_seq_output_path}")
        else:
            logger.silent_log(f"Filtered sequence written to: { outputs.seqkit_seq_output_path}")

        ## Translate nucleotide sequences to protein sequences with seqkit
        if verbose:
            logger.loud_log("Translating nucleotide sequences to protein sequences.")
        else:
            logger.silent_log("Translating nucleotide sequences to protein sequences.")

        if not os.path.exists(outputs.seqkit_translate_output_dir):
            outputs.seqkit_translate_output_dir.mkdir(parents=True)

        run_seqkit.seqkit( outputs.seqkit_seq_output_path, outputs.seqkit_translate_output_path, log_file, threads=cpus).run_seqkit_translate(gen_code, 6)

        if verbose:
            logger.loud_log(f"Translated sequence written to: {outputs.seqkit_translate_output_path}")
        else:
            logger.silent_log(f"Translated sequence written to: {outputs.seqkit_translate_output_path}")

        for db_name,db_path in zip(db_name_list, db_path_list):

            if verbose:
                logger.loud_log(f"HMM output path: {outputs.hmm_output_path(db_name)}")
            else:
                logger.silent_log(f"HMM output path: {outputs.hmm_output_path(db_name)}")

            start_hmmsearch_time = logger.start_timer()
            run_pyhmmer.pyhmmsearch(outputs.hmm_output_path(db_name), outputs.seqkit_translate_output_path, db_path, cpus, e, incdomE, domE, incE,
                                              z).run_pyhmmsearch()
            end_hmmsearch_time = logger.stop_timer(verbose)
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

            format_pyhmmer_out.hmmsearch_format_helpers(outputs.formatted_hmm_output_path(db_name), seq_type).highest_bitscore_hits(
                outputs.best_hit_path(db_name))
            if verbose:
                logger.loud_log(f"Highest Bitscore hits written to: {outputs.best_hit_path(db_name)}")
            else:
                logger.silent_log(f"Highest Bitscore hits written to: {outputs.best_hit_path(db_name)}")

            set_dict[db_name] = format_pyhmmer_out.hmmsearch_format_helpers(outputs.formatted_hmm_output_path(db_name),
                                                                            seq_type).hmm_to_contig_set()
            translated_set_dict[db_name] = format_pyhmmer_out.hmmsearch_format_helpers(outputs.formatted_hmm_output_path(db_name),
                                                                                       'prot').hmm_to_contig_set()

            # Convert to pandas dataframe, add db_name column and append to df_list
            df = pd.read_csv(outputs.best_hit_path(db_name), sep='\t')
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


        combined_df.to_csv(outputs.combined_tsv_path, sep="\t", index=False)

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

        # Write a fasta file with all the contigs
        if not os.path.exists(outputs.fasta_output_dir):
            outputs.fasta_output_dir.mkdir(parents=True)

        utils.fasta(input_file).write_fasta(utils.fasta(input_file).extract_contigs(combined_set), outputs.fasta_nuc_out_path)
        if not os.path.exists(outputs.gff_output_dir):
            outputs.gff_output_dir.mkdir(parents=True)
        format_pyhmmer_out.hmmsearch_output_writter().write_hmmsearch_hits(outputs.combined_tsv_path, seq_type, outputs.rdrpcatch_output, outputs.gff_output_path)
        rdrp_coords_list = format_pyhmmer_out.hmmsearch_output_writter().get_rdrp_coords(outputs.rdrpcatch_output)
        utils.fasta(outputs.seqkit_translate_output_path).write_fasta_coords(rdrp_coords_list,outputs.fasta_prot_out_path, seq_type)

        if verbose:
            logger.loud_log(f"Contigs written to: {outputs.fasta_nuc_out_path}")
            logger.loud_log(f"Translated contigs written to: {outputs.fasta_prot_out_path}")
        else:
            logger.silent_log(f"Contigs written to: {outputs.fasta_nuc_out_path}")
            logger.silent_log(f"Translated contigs written to: {outputs.fasta_prot_out_path}")

        if not os.path.exists(outputs.mmseqs_tax_output_dir):
            outputs.mmseqs_tax_output_dir.mkdir(parents=True)

        if verbose:
            logger.loud_log("Running mmseqs easy-taxonomy for taxonomic annotation.")
        else:
            logger.silent_log("Running mmseqs easy-taxonomy for taxonomic annotation.")

        mmseqs_tax.mmseqs(outputs.fasta_prot_out_path, mmseqs_db_path, outputs.mmseqs_tax_output_prefix,
                          outputs.mmseqs_tax_output_dir, 7, cpus, outputs.mmseqs_tax_log_path).run_mmseqs_easy_tax_lca()

        if verbose:
            logger.loud_log("Running mmseqs easy-search for taxonomic annotation.")
        else:
            logger.silent_log("Running mmseqs easy-search for taxonomic annotation.")

        if not os.path.exists(outputs.mmseqs_e_search_output_dir):
            outputs.mmseqs_e_search_output_dir.mkdir(parents=True)


        mmseqs_tax.mmseqs(outputs.fasta_prot_out_path, mmseqs_db_path, outputs.mmseqs_e_search_output_dir,
                          outputs.mmseqs_e_search_output_path, 7, cpus, outputs.mmseqs_e_search_log_path).run_mmseqs_e_search()

        utils.mmseqs_parser(outputs.mmseqs_tax_output_lca_path, outputs.mmseqs_e_search_output_path).tax_to_rdrpcatch(
            outputs.rdrpcatch_output, outputs.extended_rdrpcatch_output, seq_type)


    elif seq_type == 'prot':

        if verbose:
            logger.loud_log("Protein sequence detected.")
        else:
            logger.silent_log("Protein sequence detected.")

        set_dict = {}
        df_list = []

        for db_name,db_path in zip (db_name_list, db_path_list):

            if verbose:
                logger.loud_log(f"HMM output path: {outputs.hmm_output_path(db_name)}")
            else:
                logger.silent_log(f"HMM output path: {outputs.hmm_output_path(db_name)}")
            start_hmmsearch_time = logger.start_timer()
            hmm_out = run_pyhmmer.pyhmmsearch(outputs.hmm_output_path(db_name), input_file, db_path, cpus, e, incdomE, domE, incE, z).run_pyhmmsearch()
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

            format_pyhmmer_out.hmmsearch_formatter(hmm_out, outputs.formatted_hmm_output_path(db_name), seq_type)
            if verbose:
                logger.loud_log(f"Formatted Pyhmmer output written to: {outputs.formatted_hmm_output_path(db_name)}")
            else:
                logger.silent_log(f"Formatted Pyhmmer output written to: {outputs.formatted_hmm_output_path(db_name)}")

            # Extract Highest Bitscore hits from the formatted hmm output

            if not os.path.exists(outputs.best_hit_dir):
                outputs.best_hit_dir.mkdir(parents=True)

            format_pyhmmer_out.hmmsearch_format_helpers(outputs.formatted_hmm_output_path(db_name),seq_type).highest_bitscore_hits(outputs.best_hit_path(db_name))

            if verbose:
                logger.loud_log(f"Highest Bitscore hits written to: {outputs.best_hit_path(db_name)}")
            else:
                logger.silent_log(f"Highest Bitscore hits written to: {outputs.best_hit_path(db_name)}")

            set_dict[db_name] = format_pyhmmer_out.hmmsearch_format_helpers(outputs.formatted_hmm_output_path(db_name),seq_type).hmm_to_contig_set()

            # Convert to pandas dataframe, add db_name column and append to df_list
            df = pd.read_csv(outputs.best_hit_path(db_name), sep='\t')
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

        combined_df.to_csv(outputs.combined_tsv_path, sep="\t", index=False)

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
        if not os.path.exists(outputs.fasta_output_dir):
            outputs.fasta_output_dir.mkdir(parents=True)

        utils.fasta(input_file).write_fasta(utils.fasta(input_file).extract_contigs(combined_set), outputs.fasta_prot_out_path)

        if verbose:
            logger.loud_log(f"Contigs written to: {outputs.fasta_prot_out_path}")
        else:
            logger.silent_log(f"Contigs written to: {outputs.fasta_prot_out_path}")

        if not os.path.exists(outputs.gff_output_dir):
            outputs.gff_output_dir.mkdir(parents=True)

        format_pyhmmer_out.hmmsearch_output_writter().write_hmmsearch_hits(outputs.combined_tsv_path, seq_type, outputs.rdrpcatch_output,outputs.gff_output_path)
        rdrp_coords_list = format_pyhmmer_out.hmmsearch_output_writter().get_rdrp_coords(outputs.rdrpcatch_output)
        utils.fasta(input_file).write_fasta_coords(rdrp_coords_list,outputs.fasta_prot_out_path, seq_type)

        if verbose:
            logger.loud_log(f"RdRpCATCH output file written to: {outputs.fasta_prot_out_path}")
        else:
            logger.silent_log(f"RdRpCATCH output file written to: {outputs.fasta_prot_out_path}")

        if not os.path.exists(outputs.mmseqs_tax_output_dir):
            outputs.mmseqs_tax_output_dir.mkdir(parents=True)

        if verbose:
            logger.loud_log("Running mmseqs easy-taxonomy for taxonomic annotation.")
        else:
            logger.silent_log("Running mmseqs easy-taxonomy for taxonomic annotation.")


        mmseqs_tax.mmseqs(outputs.fasta_prot_out_path, mmseqs_db_path, outputs.mmseqs_tax_output_prefix,
                          outputs.mmseqs_tax_output_dir, 7, cpus, outputs.mmseqs_tax_log_path).run_mmseqs_easy_tax_lca()

        if not os.path.exists(outputs.mmseqs_e_search_output_dir):
            outputs.mmseqs_e_search_output_dir.mkdir(parents=True)

        if verbose:
            logger.loud_log("Running mmseqs easy-search for taxonomic annotation.")
        else:
            logger.silent_log("Running mmseqs easy-search for taxonomic annotation.")

        mmseqs_tax.mmseqs(outputs.fasta_prot_out_path, mmseqs_db_path, outputs.mmseqs_e_search_output_dir,
                          outputs.mmseqs_e_search_output_path, 7, cpus, outputs.mmseqs_e_search_log_path).run_mmseqs_e_search()

        utils.mmseqs_parser(outputs.mmseqs_tax_output_lca_path, outputs.mmseqs_e_search_output_path).tax_to_rdrpcatch(
            outputs.rdrpcatch_output, outputs.extended_rdrpcatch_output, seq_type)


if __name__ == "__main__":
    main()
