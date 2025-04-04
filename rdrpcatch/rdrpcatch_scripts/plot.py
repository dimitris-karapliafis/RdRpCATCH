import warnings
# Filter numpy warnings before any imports that might trigger them
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="numpy")
warnings.filterwarnings("ignore", message=".*subnormal.*")

class Plotter:

    def __init__(self,  upset_outdir, tsv_outdir, prefix):
        self.upset_outdir = upset_outdir
        self.tsv_outdir = tsv_outdir
        self.prefix = prefix

    def upset_plotter(self, analysis_dict):
        ''' Create an upset plot for the analysis results for a given e-value threshold

        :param analysis_dict:
        :param general_outdir:
        :param eval:
        :return:

        '''
        from matplotlib import pyplot as plt
        import upsetplot
        import os

        upset_data = upsetplot.from_contents(analysis_dict)
        # write upset data to a tsv file
        upset_data.to_csv(os.path.join(self.tsv_outdir, f"{self.prefix}_upset_data.tsv"), sep="\t")
        upsetplot.UpSet(upset_data, subset_size="count", show_counts=True, sort_by='cardinality').plot()
        plt.savefig(os.path.join(self.upset_outdir, f"{self.prefix}_upset_plot.png"), bbox_inches='tight', dpi=300)
        plt.close()

    # def plot_evalue(self, combined_df):
    #
    #     sns.set(style="whitegrid")
    #     plt.figure(figsize=(10, 6))
    #     ax = sns.boxplot(x='db_name', y='E-value', data=combined_df, showfliers=False)
    #     plt.title(f"E-value distribution", fontweight='bold')
    #     plt.savefig(os.path.join(self.upset_outdir, f"{self.prefix}_evalue_plot.png"), bbox_inches='tight', dpi=300)
    #     plt.close()

    def plot_evalue(self, combined_df):
        import polars as pl
        import altair as alt
        import os
        
        # Ensure the E-value column contains only positive numbers and convert to log scale
        df = combined_df.filter(pl.col('E-value') > 0).with_columns([
            pl.col('E-value').log10().alias('log10_evalue')
        ])
        
        chart = alt.Chart(df).mark_boxplot().encode(
            x=alt.X('db_name:N', title='Database'),
            y=alt.Y('log10_evalue:Q', title='log10(E-value)'),
            color='db_name:N'
        ).properties(
            title='E-value Distribution',
            width=600,
            height=400
        )
        
        chart.save(os.path.join(self.upset_outdir, f"{self.prefix}_evalue_plot.html"))

    def plot_score(self, combined_df):
        import altair as alt
        import os
        
        chart = alt.Chart(combined_df).mark_boxplot().encode(
            x=alt.X('db_name:N', title='Database'),
            y=alt.Y('score:Q', title='Score'),
            color='db_name:N'
        ).properties(
            title='Score Distribution',
            width=600,
            height=400
        )
        
        chart.save(os.path.join(self.upset_outdir, f"{self.prefix}_score_plot.html"))

    def plot_norm_bitscore_profile(self, combined_df):
        import altair as alt
        import os
        
        chart = alt.Chart(combined_df).mark_boxplot().encode(
            x=alt.X('db_name:N', title='Database'),
            y=alt.Y('norm_bitscore_profile:Q', title='Normalized Bitscore (Profile)'),
            color='db_name:N'
        ).properties(
            title='Normalized Bitscore Distribution (Profile)',
            width=600,
            height=400
        )
        
        chart.save(os.path.join(self.upset_outdir, f"{self.prefix}_norm_bitscore_plot_profile.html"))

    def plot_norm_bitscore_contig(self, combined_df):
        import altair as alt
        import os
        
        chart = alt.Chart(combined_df).mark_boxplot().encode(
            x=alt.X('db_name:N', title='Database'),
            y=alt.Y('norm_bitscore_contig:Q', title='Normalized Bitscore (Contig)'),
            color='db_name:N'
        ).properties(
            title='Normalized Bitscore Distribution (Contig)',
            width=600,
            height=400
        )
        
        chart.save(os.path.join(self.upset_outdir, f"{self.prefix}_norm_bitscore_contig_plot.html"))

    def plot_ID_score(self, combined_df):
        import altair as alt
        import os
        
        chart = alt.Chart(combined_df).mark_boxplot().encode(
            x=alt.X('db_name:N', title='Database'),
            y=alt.Y('ID_score:Q', title='Identity Score'),
            color='db_name:N'
        ).properties(
            title='Identity Score Distribution',
            width=600,
            height=400
        )
        
        chart.save(os.path.join(self.upset_outdir, f"{self.prefix}_ID_score_plot.html"))

    def plot_profile_coverage(self, combined_df):
        import altair as alt
        import os
        
        chart = alt.Chart(combined_df).mark_boxplot().encode(
            x=alt.X('db_name:N', title='Database'),
            y=alt.Y('profile_coverage:Q', title='Profile Coverage'),
            color='db_name:N'
        ).properties(
            title='Profile Coverage Distribution',
            width=600,
            height=400
        )
        
        chart.save(os.path.join(self.upset_outdir, f"{self.prefix}_profile_coverage_plot.html"))

    def plot_contig_coverage(self, combined_df):
        import altair as alt
        import os
        
        chart = alt.Chart(combined_df).mark_boxplot().encode(
            x=alt.X('db_name:N', title='Database'),
            y=alt.Y('contig_coverage:Q', title='Contig Coverage'),
            color='db_name:N'
        ).properties(
            title='Contig Coverage Distribution',
            width=600,
            height=400
        )
        
        chart.save(os.path.join(self.upset_outdir, f"{self.prefix}_contig_coverage_plot.html"))






