import polars as pl
import re
from pathlib import Path



def calculate_true_coverage(starts: list, ends: list) -> int:
    """Optimized coverage calculation using interval merging"""
    if not starts:
        return 0

    intervals = sorted(zip(starts, ends))
    merged = []
    current_start, current_end = intervals[0]

    for start, end in intervals[1:]:
        if start <= current_end + 1:  # Handle adjacent ranges
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end

    merged.append((current_start, current_end))
    return sum(end - start + 1 for start, end in merged)


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
        self.hmm_output_file = hmm_raw
        hmm_custom = str(hmm_raw.with_suffix('.custom.tsv'))

        # Parse and process the data using Polars DataFrame operations
        data_df = pl.read_csv(hmm_custom, separator='\t')
        data_df = self.calculate_norm_bitscore_profile(data_df)
        data_df = self.calculate_norm_bitscore_contig(data_df)
        data_df = self.calculate_coverage_stats(data_df)

        
        if seq_type == 'prot':
            self.export_processed_file_aa(data_df, hmm_processed)
        elif seq_type == 'nuc':
            self.export_processed_file_dna(data_df, hmm_processed)


    def calculate_norm_bitscore_profile(self, data_df):
        """
        Calculates the normalized bitscore for each profile.

        :param data_df: Dictionary containing the parsed data.
        :type data: dict
        :return: Dictionary containing the parsed data with normalized bitscores.
        :rtype: dataframe
        """
        data_df = (data_df.with_columns([
                        # Normalized bitscores
                        (pl.col('score') / pl.col('qlen')).alias('norm_bitscore_profile')]))
        return data_df

    def calculate_norm_bitscore_contig(self, data_df):
        """
        Calculates the normalized bitscore for each contig.

        :param data_df: Dictionary containing the parsed data.
        :type data: dict
        :return: Dictionary containing the parsed data with normalized bitscores.
        :rtype: dataframe
        """
        data_df = (data_df.with_columns([
                        # Normalized bitscores
                        (pl.col('score') / pl.col('tlen')).alias('norm_bitscore_contig')]))
        return data_df

    def calculate_coverage_stats(self, data_df):
        """
        Calculates the coverage statistics for each profile.

        :param data: Dictionary containing the parsed data.
        :type data: dict
        :return: Dictionary containing the parsed data with coverage statistics.
        :rtype: dict
        """

        df = data_df.with_columns(
            pl.col("env_from").cast(pl.Int64),
            pl.col("env_to").cast(pl.Int64),
            pl.col("hmm_from").cast(pl.Int64),
            pl.col("hmm_to").cast(pl.Int64),
            pl.col("ali_from").cast(pl.Int64),
            pl.col("ali_to").cast(pl.Int64)
        )

        stats_df = (
            df
            .with_row_index("row_id")
            .join(
                df.group_by(["t_name", "q_name"])
                .agg(
                    pl.col("env_from").alias("starts"),
                    pl.col("env_to").alias("ends"),
                    pl.col("hmm_from").alias("hmm_starts"),
                    pl.col("hmm_to").alias("hmm_ends"),
                    pl.col("ali_from").alias("ali_starts"),
                    pl.col("ali_to").alias("ali_ends"),
                    pl.col("tlen").first().alias("tlen"),
                    pl.col("qlen").first().alias("qlen"),
                    pl.col("score").first().alias("score"),
                    pl.col("env_from").min().alias("RdRp_start"),
                    pl.col("env_to").max().alias("RdRp_end"),
                    pl.len().alias("row_count")
                )
                .with_columns(
                    contig_coverage=pl.when(pl.col("row_count") == 1)
                    .then(pl.col("ends").list.first() - pl.col("starts").list.first() + 1)
                    .otherwise(
                        pl.struct(["starts", "ends"])
                        .map_elements(lambda x: calculate_true_coverage(x["starts"], x["ends"]),return_dtype=pl.Int64)
                    ),
                    profile_coverage=pl.when(pl.col("row_count") == 1)
                    .then(pl.col("hmm_ends").list.first() - pl.col("hmm_starts").list.first() + 1)
                    .otherwise(
                        pl.struct(["hmm_starts", "hmm_ends"])
                        .map_elements(lambda x: calculate_true_coverage(x["hmm_starts"], x["hmm_ends"]),return_dtype=pl.Int64)
                    ),
                    aligned_coverage=pl.when(pl.col("row_count") == 1)
                    .then(pl.col("ali_ends").list.first() - pl.col("ali_starts").list.first() + 1)
                    .otherwise(
                        pl.struct(["ali_starts", "ali_ends"])
                        .map_elements(lambda x: calculate_true_coverage(x["ali_starts"], x["ali_ends"]),return_dtype=pl.Int64)
                    )
                )
                .with_columns(
                    contig_coverage=(pl.col("contig_coverage") / pl.col("tlen")).alias("contig_coverage"),
                    profile_coverage=(pl.col("profile_coverage") / pl.col("qlen")).alias("profile_coverage"),
                    ID_score=(pl.col("score") / pl.col("aligned_coverage")).alias("ID_score")
                )
                .select(
                    ["t_name", "q_name", "contig_coverage", "profile_coverage", "ID_score", "RdRp_start", "RdRp_end"]),
                on=["t_name", "q_name"]
            )
            .sort("row_id")
            .drop("row_id")
        )
        # Group by contig and profile name, keep the first occurrence of all columns
        stats_df = (
            stats_df
            .group_by(["t_name", "q_name"])
            .agg(
                pl.col("*").first()  # Keep the first occurrence of all columns
            )
            .sort(["t_name", "q_name"])
        )

        return stats_df


    def export_processed_file_aa(self, data_df, outfile):
        """
        Exports the processed hmmsearch output file for protein sequences.

        :param data_df: Polars DataFrame containing the parsed data.
        :type data_df: pl.DataFrame
        :param outfile: Path to the output file.
        :type outfile: str
        :return: None
        """
        # Select and rename columns for output
        output_df = data_df.select([
            pl.col('t_name').alias('Contig_name'),
            pl.lit("-").alias('Translated_contig_name (frame)'),
            pl.col('tlen').alias('Sequence_length(AA)'),
            pl.col('q_name').alias('Profile_name'),
            pl.col('qlen').alias('Profile_length'),
            pl.col('E-value'),
            pl.col('score'),
            # pl.col("acc").alias("hmm_accuracy"),
            pl.col('norm_bitscore_profile'),
            pl.col('norm_bitscore_contig'),
            pl.col('ID_score'),
            pl.col('RdRp_start').alias('RdRp_from(AA)'),
            pl.col('RdRp_end').alias('RdRp_to(AA)'),
            pl.col('profile_coverage'),
            pl.col('contig_coverage')
        ])
        
        output_df.write_csv(outfile, separator="\t")

    def export_processed_file_dna(self, data_df, outfile):
        """
        Exports the processed hmmsearch output file for DNA sequences.

        :param data_df: Polars DataFrame containing the parsed data.
        :type data_df: pl.DataFrame
        :param outfile: Path to the output file.
        :type outfile: str
        :return: None
        """
        # Extract contig name and frame from translated sequence name
        output_df = (data_df
            .with_columns([
                pl.col('t_name').str.extract(r'(.+)_frame=[-]?\d').alias('Contig_name'),
                pl.col('t_name').alias('Translated_contig_name (frame)')
            ])
            .select([
                pl.col('Contig_name'),
                pl.col('Translated_contig_name (frame)'),
                pl.col('tlen').alias('Sequence_length(AA)'),
                pl.col('q_name').alias('Profile_name'),
                pl.col('qlen').alias('Profile_length'),
                pl.col('E-value'),
                pl.col('score'),
                # pl.col("acc").alias("hmm_accuracy"),
                pl.col('norm_bitscore_profile'),
                pl.col('norm_bitscore_contig'),
                pl.col('ID_score'),
                pl.col('RdRp_start').alias('RdRp_from(AA)'),
                pl.col('RdRp_end').alias('RdRp_to(AA)'),
                pl.col('profile_coverage'),
                pl.col('contig_coverage')
            ]))
        
        output_df.write_csv(outfile, separator="\t")

class hmmsearch_format_helpers:

    def __init__(self, hmm_outfn, seq_type, logger=None):
        self.hmm_outfn = hmm_outfn
        self.seq_type = seq_type
        self.logger = logger

    def hmm_to_contig_set(self):
        """
        Returns a set of all contig names in the data.

        :return: Set of contig names.
        :rtype: set
        """
        df = pl.read_csv(self.hmm_outfn, separator='\t')
        result = set(df['Contig_name'].unique())
        if self.logger:
            self.logger.silent_log(f"Found {len(result)} unique contigs")
        return result

    def highest_bitscore_hits(self, filtered_file):
        """
        Filters the hmmsearch output file based on the highest bitscore for each contig.

        :param filtered_file: Path to the filtered output file.
        :type filtered_file: str
        :return: None
        """
        df = pl.read_csv(self.hmm_outfn, separator='\t')
        if self.logger:
            self.logger.silent_log(f"Processing {len(df)} hits for highest bitscore")
        
        # Get total hits per contig
        hit_counts = df.group_by('Contig_name').agg(
            pl.count().alias('Total_positive_profiles')
        )
        
        # Get best hits by score
        best_hits = df.join(hit_counts, on='Contig_name').sort('score', descending=True).group_by('Contig_name').first()
        
        if self.logger:
            self.logger.silent_log(f"Found {len(best_hits)} best hits")
        
        best_hits.write_csv(filtered_file, separator='\t')

    def highest_norm_bit_prof_hits(self, filtered_file):
        """
        Filters the hmmsearch output file based on the highest normalized bitscore for each contig.

        :param filtered_file: Path to the filtered output file.
        :type filtered_file: str
        :return: None
        """
        df = pl.read_csv(self.hmm_outfn, separator='\t')
        if self.logger:
            self.logger.silent_log(f"Processing {len(df)} hits for highest normalized bitscore")
        
        # Get best hits by normalized bitscore
        best_hits = df.sort('norm_bitscore_profile', descending=True).group_by('Contig_name').first()
        
        if self.logger:
            self.logger.silent_log(f"Found {len(best_hits)} best hits")
        
        best_hits.write_csv(filtered_file, separator='\t')

    def lowest_evalue_hits(self, filtered_file):
        """
        Filters the hmmsearch output file based on the lowest E-value for each contig.

        :param filtered_file: Path to the filtered output file.
        :type filtered_file: str
        :return: None
        """
        df = pl.read_csv(self.hmm_outfn, separator='\t')
        if self.logger:
            self.logger.silent_log(f"Processing {len(df)} hits for lowest E-value")
        
        # Get best hits by lowest E-value
        best_hits = df.sort('E-value').group_by('Contig_name').first()
        
        if self.logger:
            self.logger.silent_log(f"Found {len(best_hits)} best hits")
        
        best_hits.write_csv(filtered_file, separator='\t')

    def extract_col(self, index):
        """
        Extracts a column from the hmmsearch output file based on index.

        :param index: Index of the column to extract.
        :type index: int
        :return: List of values from the specified column.
        :rtype: list
        """
        df = pl.read_csv(self.hmm_outfn, separator='\t')
        return df.select(df.columns[index]).to_series().to_list()

class hmmsearch_output_writter:

    def __init__(self, logger=None):
        """
        Constructor for the hmmsearch_output_writter class.
        
        :param logger: Logger instance for output
        :type logger: utils.Logger
        """
        self.logger = logger

    def write_hmmsearch_hits(self, hmmsearch_out_file, seq_type, rdrpcatch_out, gff_out):
        """
        Writes the hmmsearch hits to a GFF file.

        :param hmmsearch_out_file: Path to the hmmsearch output file.
        :type hmmsearch_out_file: str
        :param seq_type: Type of sequence (prot or nuc).
        :type seq_type: str
        :param rdrpcatch_out: Path to the RdRpCATCH output file.
        :type rdrpcatch_out: str
        :param gff_out: Path to the GFF output file.
        :type gff_out: str
        :return: None
        """
        from .utils import write_combined_results_to_gff, convert_record_to_gff3_record
        
        df = pl.read_csv(hmmsearch_out_file, separator='\t')
        
        # Write the RdRpCATCH output file first
        df.write_csv(rdrpcatch_out, separator='\t')
        
        # Create GFF format with attributes as a struct
        write_combined_results_to_gff(gff_out, df)
        # print(df.columns)
        # gff_df = df.with_columns([
        #     pl.col('Contig_name'),
        #     pl.col('db_name').alias('source'),
        #     pl.lit('protein_match').alias('type'),
        #     pl.col('RdRp_from(AA)'),
        #     pl.col('RdRp_to(AA)'),
        #     pl.col('score'),
        #     pl.lit('+').alias('strand'),
        #     pl.lit('.').alias('phase')])
        # print(gff_df)
        # print(gff_df.columns)
        # with open(gff_out, 'w') as out_handle:
        #     out_handle.write('##gff-version 3\n')
        #     for row in gff_df.iter_rows(named=True):
        #         # print(row)
        #         # print(row['Contig_name'])
        #         gff_line = "\t".join(
        #             [row['Contig_name'],
        #              row['source'],
        #              row['type'],
        #              row['RdRp_from(AA)'],
        #              row['RdRp_to(AA)'],
        #              row['score'],
        #              row['strand'],
        #              row['phase'],
        #              row['attributes']])
        #         out_handle.write(f"{gff_line}\n")
        # gff_df = gff_df.with_columns([
        #     pl.struct([
        #         pl.col('Contig_name'),
        #         pl.col('Profile_name'),
        #         pl.col('E-value').cast(pl.Utf8),
        #         pl.col('score').cast(pl.Utf8),
        #         pl.col('profile_coverage').cast(pl.Utf8),
        #         pl.col('contig_coverage').cast(pl.Utf8),
        #         pl.col('ID_score').cast(pl.Utf8)
        #     ]).map_elements(lambda x: f"ID=RdRp_{x[0]};Profile={x[1]};E-value={x[2]};score={x[3]};profile_coverage={x[4]};contig_coverage={x[5]};ID_score={x[6]}").alias('attributes')
        # ])
        
        # # Write GFF file
        # with open(gff_out, 'w') as out_handle:
        #     out_handle.write('##gff-version 3\n')
        #     gff_df.write_csv(out_handle, separator='\t', has_header=False)

    def get_rdrp_coords(self, rdrpcatch_out):
        """
        Gets the RdRp coordinates from the RdRpCATCH output file.

        :param rdrpcatch_out: Path to the RdRpCATCH output file.
        :type rdrpcatch_out: str
        :return: List of tuples containing contig name and RdRp coordinates.
        :rtype: list
        """
        # Convert the path to use combined.tsv instead of rdrpcatch_output.tsv
        combined_file = str(Path(rdrpcatch_out).parent / Path(rdrpcatch_out).stem.replace('_rdrpcatch_output', '_combined.tsv'))
        if self.logger:
            self.logger.silent_log(f"Reading coordinates from {combined_file}")
        
        df = pl.read_csv(combined_file, separator='\t')
        if self.logger:
            self.logger.silent_log(f"Found {len(df)} rows in combined file")
            self.logger.silent_log(f"Column names: {df.columns}")
        
        coords = df.select([
            'Translated_contig_name (frame)',
            'RdRp_from(AA)',
            'RdRp_to(AA)'
        ]).rows()
        
        if self.logger:
            self.logger.silent_log(f"Extracted {len(coords)} coordinate sets")
            self.logger.silent_log(f"First few coordinates: {coords[:3]}")
        return coords

class hmmsearch_combiner:

    def __init__(self, hmmsearch_files, combined_file, logger=None):
        """
        Constructor for the hmmsearch_combiner class.

        :param hmmsearch_files: List of paths to the hmmsearch output files.
        :type hmmsearch_files: list
        :param combined_file: Path to the combined output file.
        :type combined_file: str
        :param logger: Logger instance for output
        :type logger: utils.Logger
        """
        self.hmmsearch_files = hmmsearch_files
        self.combined_file = combined_file
        self.logger = logger
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
        # Read and process each file
        processed_dfs = []
        if self.logger:
            self.logger.silent_log(f"Processing {len(hmmsearch_files)} hmmsearch output files")
        
        for f in hmmsearch_files:
            if self.logger:
                self.logger.silent_log(f"Processing file: {f}")
            
            df = pl.read_csv(f, separator='\t')
            # Extract database name from filename
            db_name = Path(f).stem.split('_hmm_output')[0].split('_')[-1]
            
            # Add database name
            df = df.with_columns([
                pl.lit(db_name).alias('db_name')
            ])
            
            # Get total hits per contig
            hit_counts = df.groupby('Contig_name').agg(
                pl.count().alias('Total_positive_profiles')
            )
            df = df.join(hit_counts, on='Contig_name')
            
            if self.logger:
                self.logger.silent_log(f"Found {len(df)} hits for database {db_name}")
            
            processed_dfs.append(df)
        
        # Combine all processed DataFrames
        combined_df = pl.concat(processed_dfs)
        if self.logger:
            self.logger.silent_log(f"Combined {len(processed_dfs)} dataframes with total {len(combined_df)} rows")
        
        # Write combined DataFrame to file
        combined_df.write_csv(combined_file, separator='\t')
        if self.logger:
            self.logger.silent_log(f"Written combined results to: {combined_file}")
        
        return combined_file

