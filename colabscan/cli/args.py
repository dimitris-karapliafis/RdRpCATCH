import rich_click as click
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.syntax import Syntax
from pathlib import Path
from ..colabscanner_wrapper import run_scan, run_download

console = Console()
def parse_comma_separated_options(ctx, param, value):
    if not value:
        return ['all']

    allowed_choices = ['RVMT', 'NeoRdRp', 'NeoRdRp.2.1', 'TSA_Olendraite_fam', 'TSA_Olendraite_gen', 'RDRP-scan',
                       'Lucaprot', 'all']
    lower_choices = [choice.lower() for choice in allowed_choices]
    options = value.split(',')
    lower_options = [option.lower() for option in options]

    for option in options:
        if option.lower() not in lower_choices:
            raise click.BadParameter(f"Invalid choice: '{option}' (choose from {', '.join(allowed_choices)})")

    return lower_options


@click.group()
def cli():
    """ColabScan: A package for scanning sequences for RNA virus RNA dependent RNA Polymerases."""
    pass


@cli.command("scan", help="Scan sequences for RdRps.")
@click.option("-i", "--input",
              help="Path to the input FASTA file.",
              type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path), required=True)
@click.option("-o", "--output",
              help="Path to the output directory.",
              type=click.Path(exists=False, file_okay=False, writable=True, path_type=Path), required=True)
@click.option("-db_dir", "--db_dir",
              help="Path to the directory containing colabscan databases.",
              type=click.Path(exists=True, dir_okay=True, readable=True, path_type=Path),required=True)
@click.option("-dbs", "--db_options",
              callback=parse_comma_separated_options,
              default="all",
              help="Comma-separated list of databases to search against. Valid options: RVMT, NeoRdRp, NeoRdRp.2.1,"
                   " TSA_Olendraite_fam, TSA_Olendraite_gen, RDRP-scan,Lucaprot, all")
@click.option("-seq_type", "--seq_type",
              type=click.STRING,
              default=None,
              help="Type of sequence to search against: (prot,nuc) Default: unknown")
@click.option("-v", "--verbose",
              is_flag=True,
              help="Print verbose output.")
@click.option('-e', '--evalue',
              type=click.FLOAT,
              default=1e-5,
              help="E-value threshold for HMMsearch. (default: 1e-5)")
@click.option('-incE', '--incevalue',
              type=click.FLOAT,
              default=1e-5,
              help="Inclusion E-value threshold for HMMsearch. (default: 1e-5)")
@click.option('-domE', '--domevalue',
              type=click.FLOAT,
              default=1e-5,
              help="Domain E-value threshold for HMMsearch. (default: 1e-5)")
@click.option('-incdomE', '--incdomevalue',
              type=click.FLOAT,
              default=1e-5,
              help="Inclusion domain E-value threshold for HMMsearch. (default: 1e-5)")
@click.option('-z', '--zvalue',
              type=click.INT,
              default=1000000,
              help="Number of sequences to search against. (default: 1000000)")
@click.option('-cpus', '--cpus',
              type=click.INT,
              default=1,
              help="Number of CPUs to use for HMMsearch. (default: 1)")
@click.option('-length_thr', '--length_thr',
              type=click.INT,
              default=400,
              help="Minimum length threshold for seqkit seq. (default: 400)")
@click.option('-gen_code', '--gen_code',
              type=click.INT,
              default=1,
              help="Genetic code to use for translation. (default: 1)")
@click.pass_context
def scan(ctx, input, output, db_options, db_dir, seq_type, verbose, evalue,
         incevalue, domevalue, incdomevalue, zvalue, cpus, length_thr, gen_code):
    """Scan sequences for RdRps."""


    # Create a rich table for displaying parameters
    table = Table(title="Scan Parameters")
    table.add_column("Parameter", style="cyan")
    table.add_column("Value", style="green")

    table.add_row("Input File", str(input))
    table.add_row("Output Directory", str(output))
    table.add_row("Databases", ", ".join(db_options))
    table.add_row("Database Directory", str(db_dir))
    table.add_row("Sequence Type", seq_type or "unknown")
    table.add_row("Verbose Mode", "ON" if verbose else "OFF")
    table.add_row("E-value", str(evalue))
    table.add_row("Inclusion E-value", str(incevalue))
    table.add_row("Domain E-value", str(domevalue))
    table.add_row("Inclusion Domain E-value", str(incdomevalue))
    table.add_row("Z-value", str(zvalue))
    table.add_row("CPUs", str(cpus))
    table.add_row("Length Threshold (if nucleotide input, < contigs are filtered out", str(length_thr))
    table.add_row("Genetic Code", str(gen_code))


    console.print(Panel(table, title="Scan Configuration"))

    run_scan(
        input_file=input,
        output_dir=output,
        db_options=db_options,
        db_dir=db_dir,
        seq_type=seq_type,
        verbose=verbose,
        e=evalue,
        incE=incevalue,
        domE=domevalue,
        incdomE=incdomevalue,
        z=zvalue,
        cpus=cpus,
        length_thr=length_thr,
        gen_code=gen_code
    )


@cli.command("download", help="Download ColabScan databases.")
@click.option('-dest', '--destination_dir',
              help="Path to the directory to download HMM databases.",
              default=None,
              type=click.Path(exists=False, dir_okay=True, writable=True, path_type=Path), required=True)
@click.pass_context
def download(ctx, destination_dir):
    """Download ColabScan databases"""

    console.print(Panel(f"Downloading databases to: {destination_dir}", title="Download Status"))
    run_download(destination_dir=destination_dir)


# @cli.command("gui", help="Launch the GUI.")
# @click.pass_context
# def gui(ctx):
#     """Launch the GUI."""
#
#     console.print(Panel("Starting ColabScan GUI...", title="GUI Launch"))
#     run_gui()


if __name__ == '__main__':
    cli(obj={})

