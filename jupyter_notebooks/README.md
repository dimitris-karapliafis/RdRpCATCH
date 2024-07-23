# ColabScan
Making it easier for non-bioinfos to search their data for RdRps.  
   WORK IN PROGRESS!!!     

   Note! The output of is not intended to be a conclusive or definitive "last word" about anything, only the hmmsearch/search-tool output (match coordinates, alignment statistics etc) providing no metadata apart apart from the match/profile name and info about which public resource the hit came from - therefor it's the user responsibilty to go to (and maybe cite?) the original resource to look for the metadata they provided.

### Notebooks  
**[Preprocess.ipynb](Preprocess.ipynb)** - Downloads data (various RNA virus resources NeoRdRp, RVMT, RdRp-Scan, ... more to come?), bundles it, formats as HMM DBs and dumps somewhere (this repo for now, but ideally we'd have a periodical/on-change deposit in a continous Zenodo repo for archivability/reproducibilty and data versioning.  
**[ColabScan.ipynb](ColabScan.ipynb)** - the main `hmmsearch` workflow. Download the output of the preprocess notebook, gets input contigs from user, translates them to end-to-end, six frame, amino acid seqs, and runs that against the HMM DBs.


### TO DO:  
1. Write a guide on how to add more DBs
2. Add additional search tools, deduplication & lookup tables, etc. Initally probably will reuse code from RVMT (e.g. [searcher](https://github.com/UriNeri/RVMT/blob/main/Discovery_pipeline/RdRP_searchs/V3_runner_Iter_RdRp_Search_psihmmseqs2.sh) and [profiler](https://github.com/UriNeri/RVMT/blob/main/Domains_Annotation/Profiler_motifs.sh)).
3. Add heuristic steps to pre-filter the input:   
"Discard contigs with rRNA/tRNA (infernal/cmscan?) or mapped to silva 16S (bbmap) ?"  "Use standard ORF prediction instead of all poossible six-frame translations?"  
"Upload host genome/transcriptome to use for filtering *likely* non-viral contigs)?"
4. Add "Advanced" options (better control on flags for whatever search tool used)
5. Write a user guide /tutorial. 


### Resources currently in use
NeoRdRp  - [GitHub](https://github.com/shoichisakaguchi/NeoRdRp), [Paper](https://doi.org/10.1264/jsme2.ME22001)  
RVMT - [Zenodo](https://zenodo.org/record/7368133), [GitHub](https://github.com/UriNeri/RVMT), [Paper](https://doi.org/10.1016/j.cell.2022.08.023).  
RdRp-Scan - [GitHub](https://github.com/JustineCharon/RdRp-scan) [Paper](https://doi.org/10.1093/ve/veac082)    
   	  ⤷ (which IIRC incorporated PALMdb, [GitHub](https://github.com/rcedgar/palmdb), [Paper](https://doi.org/10.7717/peerj.14055)).  
TSA_Olendraite - [Paper](https://doi.org/10.1093/molbev/msad060), [Data](https://drive.google.com/drive/folders/1liPyP9Qt_qh0Y2MBvBPZQS6Jrh9X0gyZ?usp=drive_link), [Thesis](https://www.repository.cam.ac.uk/items/1fabebd2-429b-45c9-b6eb-41d27d0a90c2)
      
      
### Caveats:  
1. "More is not always better":  
*a*. It is very likely that diversifying the subject set of a search can yield more hits that could have been missed otherwise. However, more HMMs doesn't (necceraily) mean more diveristy. There can be a "dimisnhing returns" thing going on here.   
*b*. Regardless of how anyone defines true/false negatives/positives, aggregation of results is oblivious to which hits are what - so keep in mind, you may also have more false positives.
2. No guarentee that matches are ground trouth level of real actual viruses!  
 Everyone tries to reduce the number of false positives/negatives in their tool/DB/etc, but there is no (*for now!!!*)  no consensous on what is a true/false positive/negative, and different people have a different idea of what cutoffs/evidence are to be considered "stringent" and what are to be considered "promiscus". I make no claims that this "aggregation" approch would only give you true positives, even if you do use highly stringent alignment statistic cutoffs - among else because as anything in life this is a "garbge in garbge out" thing - consider having a perfect match between your input sequence to e.g. what someone thought was an RdRp but is actually a reverse transcriptase (false positive) or an nrEVE (no consenous are they false positives) - You'd need to consider the hits critacly based on your own set of definitions and maybe other lines of evidence (e.g. was the input extracted from total RNA or from a ds-RNA enrichment? are there any other viral genes in close proximity? can you find high-identity matches between your sequence and potential hosts or DNA entites?). Right now, this ColabScan idea makes no attempt to call the shots on such matters.
3. From a computation/total-runtime prespective this doesn't make any sense - but that's totally fine!    
This is intended for people who (rightfully) do not have the time to go over each and every new DB, learn how to write and run code for each, or get familiar with the various quirks of homology searches at large. Every time you intitate the ColabScan notebook it installs all the dependencies (e.g. hmmer and seqkit) from scratch (using resources) and fetching the DBs (using network bandwidth). We are reling on Google good graces assuming this wouldn't use much resources so we can get this search done via their free tier. If you plan on using this multiple times or for large assemblies, I suggest downloading and trying to run it via  [jupyter-lab](https://jupyter.org/install) in a seprate conde/venv/something. We might add a seperate code/notebook for running this locally, if there is any want for this.

### Contribute
Please do! pull-requests are 100% welcome.   
Coders - I'm a self addmitted bad programer, refactors welcome.  
Potential users - testers needed!!! this is designed for you so we need your feedback!  
Tool makers - add your tool.  
DB makers - add your DB.


### Closing remarks
I've decided to make this upon hearing from more and more wet-lab researchers that they keep missing the viruses in their samples if they are only using the tool that is the easiest to use "off the shelf". Giving way to a single tool/approach to monopolise a research niche is (IMO) dangerous. It can make other tools/approaches results' to be considered by default as less reliable just because of their lower adoption rate (and not because of actual scientific content grounds). The final nail in the coffin for making this was that I was standing in front of a poster presented by an extremely talented experimentalist (that is also familiar with bioinfo!) who sequenced her dsRNA extractions and shared that the off-the-shelf tool didn’t identify half or so of the sequences she can *experimentally demonstrate* are de-facto viruses, while running an HMM based method was able to help (even though only one such alternative approach was tested..)




