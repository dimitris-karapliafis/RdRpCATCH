import warnings
# Filter numpy warnings before any imports that might trigger them
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="numpy")
warnings.filterwarnings("ignore", message=".*subnormal.*")

import rich_click as click
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.syntax import Syntax
from rich.progress import Progress, BarColumn, TextColumn, DownloadColumn, TimeRemainingColumn
from pathlib import Path
import datetime
from ..rdrpcatch_wrapper import run_scan
from ..rdrpcatch_scripts.fetch_dbs import ZenodoDownloader, db_fetcher
import os
import shutil
import requests

console = Console()

## FUNCTIONS
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


def format_size(bytes_size: int) -> str:
    """Convert bytes to human-readable format without external dependencies"""
    units = ["B", "KB", "MB", "GB", "TB"]
    unit_idx = 0
    size = float(bytes_size)

    while size >= 1024 and unit_idx < len(units) - 1:
        size /= 1024
        unit_idx += 1

    return f"{size:.2f} {units[unit_idx]}"



## CLI ENTRY POINT

@click.group()
def cli():
    """RdRpCATCH - RNA-dependent RNA polymerase Collaborative Analysis Tool with Collections of pHMMs"""
    pass

@cli.command("scan", help="Scan sequences for RdRps.")
@click.option("-i", "--input",
              help="Path to the input FASTA file.",
              type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path), required=True)
@click.option("-o", "--output",
              help="Path to the output directory.",
              type=click.Path(exists=False, file_okay=False, writable=True, path_type=Path), required=True)
@click.option("-db_dir", "--db_dir",
              help="Path to the directory containing RdRpCATCH databases.",
              type=click.Path(exists=True, dir_okay=True, readable=True, path_type=Path),required=True)
@click.option("-dbs", "--db_options",
              callback=parse_comma_separated_options,
              default="all",
              help="Comma-separated list of databases to search against. Valid options: RVMT, NeoRdRp, NeoRdRp.2.1,"
                   " TSA_Olendraite_fam, TSA_Olendraite_gen, RDRP-scan,Lucaprot, all")
@click.option("--custom-dbs",
              help="Path to directory containing custom MSAs/pHMM files to use as additional databases",
              type=click.Path(exists=True, path_type=Path))
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
              help='Genetic code to use for translation. (default: 1) Possible genetic codes (supported by seqkit translate) :      1: The Standard Code       \n'
                     '2: The Vertebrate Mitochondrial Code      \n'
                     '3: The Yeast Mitochondrial Code       \n'
                     '4: The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code        \n'
                     '5: The Invertebrate Mitochondrial Code        \n'
                     '6: The Ciliate, Dasycladacean and Hexamita Nuclear Code       \n'
                     '9: The Echinoderm and Flatworm Mitochondrial Code     \n'
                    '10: The Euplotid Nuclear Code      \n'
                    '11: The Bacterial, Archaeal and Plant Plastid Code     \n'
                    '12: The Alternative Yeast Nuclear Code     \n'
                    '13: The Ascidian Mitochondrial Code        \n'
                    '14: The Alternative Flatworm Mitochondrial Code        \n'
                    '16: Chlorophycean Mitochondrial Code       \n'
                    '21: Trematode Mitochondrial Code       \n'
                    '22: Scenedesmus obliquus Mitochondrial Code        \n'
                    '23: Thraustochytrium Mitochondrial Code        \n'
                    '24: Pterobranchia Mitochondrial Code       \n'
                    '25: Candidate Division SR1 and Gracilibacteria Code        \n'
                    '26: Pachysolen tannophilus Nuclear Code        \n'
                    '27: Karyorelict Nuclear        \n'
                    '28: Condylostoma Nuclear       \n'
                    '29: Mesodinium Nuclear     \n'
                    '30: Peritrich Nuclear      \n'
                    '31: Blastocrithidia Nuclear        \n')
@click.option('-bundle', '--bundle',
              is_flag=True,
              default=False,
              help="Bundle the output files into a single archive. (default: False)")
@click.option('-keep_tmp', '--keep_tmp',
              is_flag=True,
              default=False,
              help="Keep temporary files (Expert users) (default: False)")
@click.pass_context
def scan(ctx, input, output, db_options, db_dir, custom_dbs, seq_type, verbose, evalue,
         incevalue, domevalue, incdomevalue, zvalue, cpus, length_thr, gen_code, bundle, keep_tmp):
    """Scan sequences for RdRps."""

    # Create a rich table for displaying parameters
    table = Table(title="Scan Parameters")
    table.add_column("Parameter", style="cyan")
    table.add_column("Value", style="green")

    table.add_row("Input File", str(input))
    table.add_row("Output Directory", str(output))
    table.add_row("Databases", ", ".join(db_options))
    table.add_row("Database Directory", str(db_dir))
    if custom_dbs:
        table.add_row("Custom Databases", str(custom_dbs))
    table.add_row("Sequence Type", seq_type or "unknown")
    table.add_row("Verbose Mode", "ON" if verbose else "OFF")
    table.add_row("E-value", str(evalue))
    table.add_row("Inclusion E-value", str(incevalue))
    table.add_row("Domain E-value", str(domevalue))
    table.add_row("Inclusion Domain E-value", str(incdomevalue))
    table.add_row("Z-value", str(zvalue))
    table.add_row("CPUs", str(cpus))
    table.add_row("Length Threshold", str(length_thr))
    table.add_row("Genetic Code", str(gen_code))
    table.add_row("Bundle Output", "ON" if bundle else "OFF")
    table.add_row("Save Temporary Files", "ON" if keep_tmp else "OFF")

    console.print(Panel(table, title="Scan Configuration"))

    # Add custom databases if provided
    if custom_dbs:
        db = db_fetcher(db_dir)
        if os.path.isfile(custom_dbs):
            db.add_custom_db(custom_dbs)
        else:
            for item in os.listdir(custom_dbs):
                item_path = os.path.join(custom_dbs, item)
                if os.path.isfile(item_path) and item_path.endswith(('.hmm', '.h3m', '.msa', '.sto', '.fasta', '.fa')):
                    db.add_custom_db(item_path)
                elif os.path.isdir(item_path):
                    db.add_custom_db(item_path, item)

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
        gen_code=gen_code,
        bundle=bundle,
        keep_tmp=keep_tmp
    )

# @cli.command("download", help="Download RdRpCATCH databases.")
# @click.option("--destination_dir", "-dest",
#               help="Path to the directory to download HMM databases.",
#               type=click.Path(exists=False, file_okay=False, writable=True, path_type=Path), required=True)
# @click.option("--check-updates", "-u",
#               is_flag=True,
#               help="Check for database updates")
# @click.pass_context
# def download(ctx, destination_dir, check_updates):
#     """Download RdRpCATCH databases."""
#
#     # if check_updates:
#     #     db = db_fetcher(destination_dir)
#     #     version_info = db.check_db_updates()
#     #     if version_info:
#     #         console.print("Current database versions:")
#     #         for db_name, info in version_info.items():
#     #             console.print(f"- {db_name}: {info}")
#     #     else:
#     #         console.print("No version information available")
#     #     return
#
#     run_download(destination_dir)
#
# # @cli.command("gui", help="Launch the GUI.")
# # @click.pass_context
# # def gui(ctx):
# #     """Launch the GUI."""
# #
# #     console.print(Panel("Starting ColabScan GUI...", title="GUI Launch"))
# #     run_gui()



@cli.command("download", help="Download &  update RdRpCATCH databases. If databases are already installed in the "
                              "specified directory,"
                              " it will check for updates and download the latest version if available.")
@click.option("--destination_dir", "-dest",
              help="Path to directory to download databases",
              type=click.Path(path_type=Path, file_okay=False, writable=True),
              required=True)
@click.option("--concept-doi", default="10.5281/zenodo.14358348",
              help="Zenodo Concept DOI for database repository")
def download(destination_dir: Path, concept_doi: str):
    """Handle database download/update workflow"""
    downloader = ZenodoDownloader(concept_doi, destination_dir)

    try:

        current_version = downloader.get_current_version()
        if downloader.lock_file.exists():
            console.print("[red]× Another download is already in progress[/red]")
            raise click.Abort()

        if downloader.needs_update() or not current_version:
            downloader.lock_file.touch(exist_ok=False)
            with Progress(
                    TextColumn("[progress.description]{task.description}"),
                    BarColumn(),
                    TextColumn("{task.completed:.2f}/{task.total:.2f} MB"),
                    TimeRemainingColumn(),
                    transient=True
            ) as progress:
                # Setup main download task
                main_task = progress.add_task("[cyan]Database Manager", total=4)

                # Phase 1: Metadata fetching
                progress.update(main_task, description="Fetching Zenodo metadata...")
                metadata = downloader._fetch_latest_metadata()
                progress.advance(main_task)

                # Phase 2: Prepare download
                progress.update(main_task, description="Analyzing package...")
                tarball_info = downloader._get_tarball_info()
                file_size_mb = tarball_info["size"] / (1024 * 1024)
                progress.advance(main_task)

                # Phase 3: Download with progress
                progress.update(main_task,
                                description="Downloading RdRpCATCH databases...",
                                total=file_size_mb)

                if not downloader.temp_dir.exists():
                    downloader.temp_dir.mkdir(parents=True, exist_ok=True)

                temp_tar = downloader.temp_dir / "download.tmp"

                with requests.get(tarball_info["url"], stream=True) as response:
                    response.raise_for_status()
                    with open(temp_tar, "wb") as f:
                        downloaded = 0
                        for chunk in response.iter_content(chunk_size=8192):
                            f.write(chunk)
                            downloaded += len(chunk)
                            progress.update(main_task, advance=len(chunk) / (1024 * 1024))

                # Phase 4: Verification & installation
                progress.update(main_task, description="Verifying checksum...")
                if not downloader._verify_checksum(temp_tar, tarball_info["checksum"]):
                    raise ValueError("Checksum verification failed")

                progress.update(main_task, description="Installing databases...")
                downloader.extract_and_verify(temp_tar)
                version_info = downloader.get_latest_version_info()
                downloader.atomic_write_version(version_info)
                progress.advance(main_task)

            # Success message
            size_str = format_size(tarball_info["size"])
            console.print(
                f"\n[bold green]✓ Successfully downloaded version {version_info['record_id']}[/bold green]",
                f"Release date: {version_info['created']}",
                f"Size: {size_str}",
                sep="\n"
            )

        else:
            installed_date = current_version["downloaded"]
            console.print(
                f"[green]✓ Databases are current[/green]",
                f"Version ID: {current_version['record_id']}",
                f"Installed: {installed_date}",
                sep="\n"
            )
    except FileExistsError:
        console.print("[red]× Another download is already in progress![/red]")
        console.print(f"Lock file exists: {downloader.lock_file}")
        raise click.Abort()

    except Exception as e:
        console.print(f"\n[red]× Download failed: {str(e)}[/red]")
        if downloader.temp_dir.exists():
            shutil.rmtree(downloader.temp_dir)
        raise click.Abort()

    finally:
        # Cleanup operations
        if downloader.lock_file.exists():
            downloader.lock_file.unlink()
        if downloader.temp_dir.exists():
            shutil.rmtree(downloader.temp_dir)


if __name__ == '__main__':
    cli(obj={})

