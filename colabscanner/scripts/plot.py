import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import upsetplot
import os


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
        # Ensure the E-value column contains only positive numbers
        combined_df['E-value'] = combined_df['E-value'].clip(lower=1e-1000)

        sns.set(style="whitegrid")
        plt.figure(figsize=(10, 6))

        ax = sns.boxplot(x='db_name', y='E-value', data=combined_df, showfliers=False)
        ax.set_yscale('log')
        ymin, ymax = ax.get_ylim()
        ax.set_ylim([max(ymin, 1e-1000), ymax])

        ax.grid(axis='y', which='major', linestyle='--', alpha=0.7)

        plt.title(f"E-value distribution", fontweight='bold')

        # Customize y-axis ticks and labels
        # yticks = [1e-200, 1e-150, 1e-100, 1e-50, 1e-20, 1e-10, 1e-5, 1e-1, 1]
        # ax.set_yticks(yticks)
        # ax.set_yticklabels([f"{y:.2e}" for y in yticks])

        plt.savefig(os.path.join(self.upset_outdir, f"{self.prefix}_evalue_plot_logscale.png"), bbox_inches='tight',
                    dpi=300)
        plt.close()


    def plot_score(self, combined_df):

        sns.set(style="whitegrid")
        plt.figure(figsize=(10, 6))
        ax = sns.boxplot(x='db_name', y='score', data=combined_df)
        plt.title(f"Bitscore distribution", fontweight='bold')
        plt.savefig(os.path.join(self.upset_outdir, f"{self.prefix}_score_plot.png"), bbox_inches='tight', dpi=300)
        plt.close()

    def plot_norm_bitscore_profile(self, combined_df):

        sns.set(style="whitegrid")
        plt.figure(figsize=(10, 6))
        ax = sns.boxplot(x='db_name', y='norm_bitscore_profile', data=combined_df)
        plt.title(f"Normalized bitscore distribution", fontweight='bold')
        plt.savefig(os.path.join(self.upset_outdir, f"{self.prefix}_norm_bitscore_plot.png"), bbox_inches='tight', dpi=300)
        plt.close()


    def plot_norm_bitscore_contig(self, combined_df):

        sns.set(style="whitegrid")
        plt.figure(figsize=(10, 6))
        ax = sns.boxplot(x='db_name', y='norm_bitscore_contig', data=combined_df)
        plt.title(f"Normalized bitscore distribution", fontweight='bold')
        plt.savefig(os.path.join(self.upset_outdir, f"{self.prefix}_norm_bitscore_plot.png"), bbox_inches='tight', dpi=300)
        plt.close()

    def plot_ID_score(self, combined_df):

        sns.set(style="whitegrid")
        plt.figure(figsize=(10, 6))
        ax = sns.boxplot(x='db_name', y='ID_score', data=combined_df)
        plt.title(f"Identity score distribution", fontweight='bold')
        plt.savefig(os.path.join(self.upset_outdir, f"{self.prefix}_ID_score_plot.png"), bbox_inches='tight', dpi=300)
        plt.close()

    def plot_profile_coverage(self, combined_df):

        sns.set(style="whitegrid")
        plt.figure(figsize=(10, 6))
        ax = sns.boxplot(x='db_name', y='profile_coverage', data=combined_df)
        plt.title(f"Profile coverage distribution", fontweight='bold')
        plt.savefig(os.path.join(self.upset_outdir, f"{self.prefix}_profile_coverage_plot.png"), bbox_inches='tight', dpi=300)
        plt.close()

    def plot_contig_coverage(self, combined_df):

        sns.set(style="whitegrid")
        plt.figure(figsize=(10, 6))
        ax = sns.boxplot(x='db_name', y='contig_coverage', data=combined_df)
        plt.title(f"Contig coverage distribution", fontweight='bold')
        plt.savefig(os.path.join(self.upset_outdir, f"{self.prefix}_contig_coverage_plot.png"), bbox_inches='tight', dpi=300)
        plt.close()






