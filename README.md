# RdRpCATCH
## RNA-dependent RNA polymerase Collaborative Analysis Tool with Collections of pHMMs



RdRpCATCH is collaborative effort to combine various publicly available RNA virus RNA-dependent RNA polymerase pHMM databases in one tool
to facilitate their detection  in (meta-)transcriptomics data.


RdRpCATCH  is written in Python and uses the pyHMMER3
library to perform pHMM searches.  In addition, the tool scans each sequence (aa or nt) in the input file with the selected databases and provides the best hit (hit with the highest bitscore across all databases) as output.
In addition, RdRpCATCH provides information about the number of profiles
that were positive for each sequence across all pHMM databases, and taxonomic information based on the MMseqs2 easy-taxonomy and search modules against a custom RefSeq Riboviria database.

### Version 0.0.7 -> 0.0.8 Changelog
- Added support for custom pHMM databases. See the [Setting up custom pHMM databases](#setting-up-custom-phmm-databases) section for more information.
- All specified flags use '-' instead of '_' (e.g. `--db-dir` instead of `--db_dir`). 
- Fixed issue with specifying the Lucaprot_HMM and Zayed_HMM databases in the `--db-options` argument.
- Command `rdrpcatch download` renamed as `rdrpcatch databases` for clarity, as it now supports adding custom pHMM
databases to the RdRpCATCH databases. This is facilitated by the `--add-custom-db` argument.
- Added none option to the `--db-options` argument to search only against custom databases.



** The tool has been modified to use [rolypoly](https://code.jgi.doe.gov/UNeri/rolypoly) code/approaches **



![rdrpcatch_flowchart_v0.png](images%2Frdrpcatch_illustration.png)

### Supported databases
- NeoRdRp <sup>1</sup> : 1182 pHMMs 
- NeoRdRp2 <sup>2</sup>: 19394 pHMMs  
- RVMT <sup>3</sup>: 710 pHMMs  
- RdRp-Scan <sup>4</sup> : 68 pHMMs
- TSA_Olendraite_fam <sup>5</sup>: 77 pHMMs 
- TSA_Olendraite_gen <sup>6</sup> : 341 pHMMs
- LucaProt_HMM<sup>7 </sup> : 754 pHMMs
- Zayed_HMM<sup>8 </sup> : 2489 pHMMs

1. Sakaguchi, S. et al. (2022) 'NeoRdRp: A comprehensive dataset for identifying RNA-dependent RNA polymerases of various RNA viruses from metatranscriptomic data', *Microbes and Environments*, 37(3). [doi:10.1264/jsme2.me22001](https://doi.org/10.1264/jsme2.me22001)
2. Sakaguchi, S., Nakano, T. and Nakagawa, S. (2024) 'Neordrp2 with improved seed data, annotations, and scoring', *Frontiers in Virology*, 4. [doi:10.3389/fviro.2024.1378695](https://doi.org/10.3389/fviro.2024.1378695)
3. Neri, U. et al. (2022) 'Expansion of the global RNA virome reveals diverse clades of bacteriophages', *Cell*, 185(21). [doi:10.1016/j.cell.2022.08.023](https://doi.org/10.1016/j.cell.2022.08.023)
4. Charon, J. et al. (2022) 'RDRP-Scan: A bioinformatic resource to identify and annotate divergent RNA viruses in metagenomic sequence data', *Virus Evolution*, 8(2). [doi:10.1093/ve/veac082](https://doi.org/10.1093/ve/veac082)
5. Olendraite, I., Brown, K. and Firth, A.E. (2023) 'Identification of RNA virus–derived rdrp sequences in publicly available transcriptomic data sets', *Molecular Biology and Evolution*, 40(4). [doi:10.1093/molbev/msad060](https://doi.org/10.1093/molbev/msad060)
6. Olendraite, I. (2021) 'Mining diverse and novel RNA viruses in transcriptomic datasets', Apollo. Available at: [https://www.repository.cam.ac.uk/items/1fabebd2-429b-45c9-b6eb-41d27d0a90c2](https://www.repository.cam.ac.uk/items/1fabebd2-429b-45c9-b6eb-41d27d0a90c2)
7. Hou, X. and He, Y. et al. (2024) 'Using artificial intelligence to document the hidden RNA virosphere', *Cell*, 187(24). [doi:10.1016/j.cell.2024.09.027](https://doi.org/10.1016/j.cell.2024.09.027)
8. Zayed, A. A., et al. (2022)  'Cryptic and abundant marine viruses at the evolutionary origins of Earth’s RNA virome.' *Science*, 376(6589), 156–162. [doi:10.1126/science.abm5847](https://doi.org/10.1126/science.abm5847)



## Installation


#### Prerequisites
For the installation process, conda is required. If you don't have conda installed, you can find instructions on how to
https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html  
Mamba is a faster alternative to conda. If you have it installed, you can use it instead of conda.  

#### Installation steps

The package is available as a bioconda package. You can install it using the following command:

```bash
conda env create rdrpcatch -c bioconda rdrpcatch
```

Alternatively, you can install RdRpCATCH from python package index (PyPI) using pip. This requires the installation of the dependencies
manually. The dependencies are:
- mmseqs2
- seqkit

The dependencies can be installed using conda or mamba. Follow these steps:

Create a new conda environment and install the dependencies:
```bash
conda env create -n rdrpcatch python=3.12
conda activate rdrpcatch
conda install -c bioconda mmseqs2==17.b804f seqkit==2.10.0
```
Install the tool from pip:
```bash
pip install rdrpcatch
```

Activate the environment and download the RdRpCATCH databases:

```bash 
conda activate rdrpcatch
rdrpcatch databases --destination_dir path/to/store/databases
```

* Note 1: The databases are large files and may take some time to download (~ 3 GB).
* Note 2: The databases are stored in the specified directory, and the path is required to run RdRpCATCH.
* Note 3: If you encounter an SSL error while downloading, please try again. The error seems to appear sporadically during testing, and a simple re-initiation of the downloading process seems to fix it. 




## Usage
RdRpCATCH can be used as a CLI tool as follows:

```bash 
# make sure the conda environment is activated
# conda activate rdrpcatch

# scan the input fasta file with the selected databases
rdrpcatch scan -i path/to/input.fasta -o path/to/output_dir -db-dir path/to/database
```

## Input description
The input file can be one or more nucleotide or protein sequences in multi-fasta format. 
The output directory is where the results will be stored. We recommend specifying the type of the sequence in the command line,
An optional argument `--seq-type` (nuc or prot) can be used to specify if the input fasta file sequences are nucleotide or amino acid.


## Setting up custom pHMM databases
It is possible to use custom pHMM databases with RdRpCATCH. As a prerequisite, you need to install the RdRpCATCH 
databases using the `rdrpcatch databases` command as described above, to a directory of your choice.

The custom databases should be formatted as follows:

- First create a directory and give it a descriptive name, e.g. `my_custom_rdrp_database`. Important: The name should not contain  comma `,` characters.
- Inside the directory put your custom pHMM HMMER pressed database. You can use the `hmmpress` command of HMMER to create the pressed database from your custom HMM file. This creates a set of files with the same name as the original HMM file, but with different extensions (e.g. `.h3f`, `.h3i`, `.h3m`, `.h3p`).  The directory should contain all these files. Please refer to the HMMER manual for more information on how to create a pressed database from an HMM file. (http://eddylab.org/software/hmmer/Userguide.pdf) 
- Next you can add the directory to the custom databases that are readable by RdRpCATCH. This can be done by using the rdrpcatch databases command as follows:

```bash
rdrpcatch databases --add-custom-db path/to/my_custom_rdrp_database --destination-dir path/that/contains/rdrpcatch/databases
```

- This will add the custom database to the list of databases that can be used with RdRpCATCH.
- The custom database can then be used with the `rdrpcatch scan` command by specifying the `--custom-dbs` argument as follows:
- 
```bash
rdrpcatch scan -i path/to/input.fasta -o path/to/output_dir -db-dir path/to/database --custom-dbs custom_database_name
```

- The `custom_database_name` should be the name of the directory that contains the custom pHMM files, without the path.
- For example, if the custom database is stored in `path/to/my_custom_rdrp_database`, you would use `--custom-dbs my_custom_rdrp_database` in the command line.
- You can add multiple custom databases by installing them in the same way and specifying them  by separating them with commas, e.g. `--custom-dbs my_custom_rdrp_database,another_custom_database`.
- The custom databases can be used in combination with the pre-compiled databases provided by RdRpCATCH. To do this, you can specify the `--db_options` argument with the names of the pre-compiled databases you want to use, and specify the custom databases with the `--custom-dbs` argument.
- For example, if you want to use the NeoRdRp and RVMT databases along with your custom database, you would use the following command:

```bash
rdrpcatch scan -i path/to/input.fasta -o path/to/output_dir -db-dir path/to/database --db-options NeoRdRp,RVMT --custom-dbs my_custom_rdrp_database
```

- Note: By default, RdRpCATCH will search against all pre-compiled databases if no `--db_options` argument is specified. If you want to use only the custom databases, you can specify `--db_options none` to avoid searching against the pre-compiled databases.




## Commands
The following two commands are available in RdRpCATCH:  
* [`rdrpcatch scan`](#rdrpcatch-scan)  
* [`rdrpcatch databases`](#rdrpcatch-download)

### rdrpcatch databases:
Command to download pre-compiled databases from Zenodo and to set up custom databases. If the databases are already downloaded in the specified directory
, the command will check for updates and download the latest version if available.

| Argument | Short Flag | Type | Description                                                 |
|----------|------------|------|-------------------------------------------------------------|
| `--destination_dir` | `-dest` | PATH | Path to the directory to download HMM databases. [required] |
| `--concept-doi` | `` | TEXT | Zenodo Concept DOI for database repository                  |
| `--help` | `` |  | Show help message and exit                                  |
| `--add-custom-db` | `` | PATH | Path to the directory containing custom pHMM files to add to the RdRpCATCH databases. |

### rdrpcatch scan:
Search a given input using selected RdRp databases.  

| Argument         | Short Flag    | Type | Description                                                                                                                                                                                  |
|------------------|---------------|------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--input`        | `-i`          | FILE | Path to the input FASTA file. [required]                                                                                                                                                     |
| `--output`       | `-o`          | DIRECTORY | Path to the output directory. [required]                                                                                                                                                     |
| `--db-dir`       | `-db-dir`     | PATH | Path to the directory containing RdRpCATCH databases. [required]                                                                                                                             |
| `--db-options`   | `-dbs`        | TEXT | Comma-separated list of pre-installed databases to search against. Valid options: RVMT, NeoRdRp, NeoRdRp.2.1, TSA_Olendraite_fam, TSA_Olendraite_gen, RDRP-scan, Lucaprot_HMM,Zayed_HMM, all |
| `--custom-dbs`   |               | PATH | Comma-separated list of custom databases to search against. Valid options: names of the directories that the custom databases are stored in.                                                 |
| `--seq-type`     | `-seq-type`   | TEXT | Type of sequence to search against: (prot,nuc) Default: unknown                                                                                                                              |
| `--verbose`      | `-v`          | FLAG | Print verbose output.                                                                                                                                                                        |
| `--evalue`       | `-e`          | FLOAT | E-value threshold for HMMsearch. (default: 1e-5)                                                                                                                                             |
| `--incevalue`    | `-incE`       | FLOAT | Inclusion E-value threshold for HMMsearch. (default: 1e-5)                                                                                                                                   |
| `--domevalue`    | `-domE`       | FLOAT | Domain E-value threshold for HMMsearch. (default: 1e-5)                                                                                                                                      |
| `--incdomevalue` | `-incdomE`    | FLOAT | Inclusion domain E-value threshold for HMMsearch. (default: 1e-5)                                                                                                                            |
| `--zvalue`       | `-z`          | INTEGER | Number of sequences to search against. (default: 1000000)                                                                                                                                    |
| `--cpus`         | `-cpus`       | INTEGER | Number of CPUs to use for HMMsearch. (default: 1)                                                                                                                                            |
| `--length-thr`   | `-length-thr` | INTEGER | Minimum length threshold for seqkit seq. (default: 400)                                                                                                                                      |
| `--gen-code`     | `-gen-code`   | INTEGER | Genetic code to use for translation. (default: 1)                                                                                                                                            |
| `--bundle`       | `-bundle`     |  | Bundle the output files into a single archive. (default: False)                                                                                                                              |
| `--keep-tmp`     | `-keep-tmp`   |  | Keep the temporary files generated during the analysis. (default: False)                                                                                                                     |



#### Output files  
rdrpcatch scan will create a folder with the following structure:

| Output | Description                                                                  |
|--------|------------------------------------------------------------------------------|
| `{prefix}_rdrpcatch_output_annotated.tsv` | A tab-separated file containing the results of the RdRpCATCH analysis.       |
| `{prefix}_rdrpcatch_fasta` | A directory containing the sequences that were identified as RdRp sequences. |
| `{prefix}_rdrpcatch_plots` | A directory containing the plots generated during the analysis.              |
| `{prefix}_gff_files` | A directory containing the GFF files generated during the analysis. (For now only based on protein sequences) |
| `tmp` | A directory containing temporary files generated during the analysis. (Only available if the -keep_tmp flag is used )|

#### Output table fields
A summary of the results is stored in the `{prefix}_rdrpcatch_output_annotated.tsv` file, which contains the following fields:
| Field | Description                                                                                                         |
|-------|---------------------------------------------------------------------------------------------------------------------|
| `Contig_name` | The name of the contig.                                                                                             |
| `Translated_contig_name (frame)` | The name of the translated contig and the frame of the RdRp sequence.                                               |
| `Sequence_length(AA)` | The length of the RdRp sequence in amino acids.                                                                     |
| `Total_databases_that_the_contig_was_detected(No_of_Profiles)` | The name of databases and the number of profiles that the RdRp sequence was detected by.                            |
| `Best_hit_Database` | The database with the best hit.                                                                                     |
| `Best_hit_profile_name` | The name of the profile with the best hit.                                                                          |
| `Best_hit_profile_length` | The length of the profile with the best hit.                                                                        |
| `Best_hit_e-value` | The e-value of the best hit.                                                                                        |
| `Best_hit_bitscore` | The bitscore of the best hit.                                                                                       |
| `RdRp_from(AA)` | The start position of the RdRp sequence, in relation to the amino acid sequence.                                    |
| `RdRp_to(AA)` | The end position of the RdRp sequence, in relation to the amino acid sequence.                                      |
| `Best_hit_profile_coverage` | The fraction of the profile that was covered by the RdRp sequence.                                                  |
| `Best_hit_contig_coverage` | The fraction of the contig that was covered by the RdRp sequence. (Based on aminoacid sequence)                     |
| `MMseqs_Taxonomy_2bLCA` | The taxonomy of the RdRp sequence based on MMseqs2 easy-taxonomy module against a custom RefSeq Riboviria database. |
| `MMseqs_TopHit_accession` | The accession of the top hit in the RefSeq Riboviria database.                                                      |
| `MMseqs_TopHit_fident` | The fraction of identical matches of the top hit in the RefSeq Riboviria database.                                  |
| `MMseqs_TopHit_alnlen` | The alignment length of the top hit in the RefSeq Riboviria database.                                               |
| `MMseqs_TopHit_eval` | The e-value of the top hit in the RefSeq Riboviria database.                                                        |
| `MMseqs_TopHit_bitscore` | The bitscore of the top hit in the RefSeq Riboviria database.                                                       |
| `MMseqs_TopHit_qcov` | The query coverage of the top hit in the RefSeq Riboviria database.                                                 |
| `MMseqs_TopHit_lineage` | The lineage of the top hit in the RefSeq Riboviria database.                                                        |

## Citations
Manuscript still in preparation. If you use RdRpCATCH, please cite this GitHub repository 
A precompiled version of the used databases is available at Zenodo DOI: [10.5281/zenodo.14358348](https://doi.org/10.5281/zenodo.14358348).  
If you use RdRpCATCH, please cite the [underlying third party databases](#supported-databases) :

## Acknowledgements
RdRpCATCH is a collaborative effort and we would like to thank all the authors and developers of the underlying databases. 

## Contact
Dimitris Karapliafis (dimitris.karapliafis@wur.nl), potentially via slack/teams or an issue in the main repo.

##TODO:
- [ ] loud logging is linking to the utils.py file, not the actual line of code causing the error.
- [ ] drop `db_dir` argument and use global/environment/config variable that is set after running the `download` command


## Contributing
TBD up to Dimitris and Anne

## Licence
[MIT](LICENSE)
