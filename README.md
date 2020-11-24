Pervasive duplication of tumor suppressors in Afrotherians during the evolution of large bodies and reduced cancer risk
=====
<a content="https://www.biorxiv.org/content/10.1101/2020.09.10.291906v1" href="https://www.biorxiv.org/content/10.1101/2020.09.10.291906v1" rel="me noopener noreferrer" style="vertical-align:center;"><img alt="bioRxiv" src="https://www.biorxiv.org/sites/default/files/bioRxiv_article.jpg"><a/>

[![DOI](https://zenodo.org/badge/315503791.svg)](https://zenodo.org/badge/latestdoi/315503791)

This is the home of the reproducible manuscript for our publication. Here, you will find the instructions for how to properly set up your system to reproduce our results.

Authors
-----
* **Juan M Vazquez**
<a itemprop="sameAs" content="https://orcid.org/0000-0001-8341-2390" href="https://orcid.org/0000-0001-8341-2390" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a>[![alt text][1.2]][mainmanmanny] [:globe_with_meridians:](https://vazquez.bio)  
*Department of Integrative Biology, University of California - Berkeley*
* **Vincent J Lynch** <a itemprop="sameAs" content="https://orcid.org/0000-0001-5311-3824" href="https://orcid.org/0000-0001-5311-3824" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a>[![alt text][1.2]][devoevomed] [:globe_with_meridians:](https://arts-sciences.buffalo.edu/biological-sciences/faculty/faculty-directory/vincent-lynch.html) [:e-mail:](mailto:vjlynch@buffalo.edu)  
*Department of Biology, SUNY Buffalo*

Abstract
----

>The risk of developing cancer is correlated with body size and lifespan within species. Between species, however, there is no correlation between cancer and either body size or lifespan, indicating that large, long-lived species have evolved enhanced cancer protection mechanisms. Elephants and their relatives (Proboscideans) are a particularly interesting lineage for the exploration of mechanisms underlying the evolution of augmented cancer resistance because they evolved large bodies recently and are closely related to smaller bodied species (Afrotherians). Here, we explore the contribution of gene duplication to body size and cancer risk in Afrotherians. Unexpectedly, we found that tumor suppressor duplication was pervasive in Afrotherian genomes, rather than restricted to Proboscideans. Proboscideans, however, have unique duplicates in pathways that may underlie some aspects of their remarkable anti-cancer cell biology. These data suggest that duplication of tumor suppressor genes facilitated the evolution of increased body size by compensating for increased cancer risk.

System Requirements
---

***Local Computer:*** The `RecSearch` and _de novo_ transcritome assembly steps are the most resource-intensive parts of this pipeline. `RecSearch` was run using workstation equipped with two 20-core processors @ 2.4 GHz and 128 GB RAM; however, `RecSearch` has a memory saver option that lowers the memory requirement to 8-16 GB, depending on the genome. The `snakemake` rule for RecSearch can be edited to use the `memory_saver_level = 2` in order to run this script on common equipment.

***Computing Clusters & Cloud Computing***: Individual parts of this pipeline were tested on the Midway2 Computing Cluster at the University of Chicago, which uses a SLURM job scheduler. Courtesy of `snakemake`, the entire pipeline should be compatible with various remote computing workflows, with caveats discussed in "**Usage**."

***Operating System***: This pipeline was developed on Linux, with additional testing on Mac computers. I would strongly urge Windows uses to set up [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10) before running this - and really any - computational analysis.

Setup
----

1. Install [`conda`](https://docs.conda.io/en/latest/miniconda.html), [`bioconda`](https://bioconda.github.io/), and [`snakemake`](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) on your computer.
1. Install `gfServer` and `gfClient` from UCSC: [http://hgdownload.soe.ucsc.edu/admin/exe/](http://hgdownload.soe.ucsc.edu/admin/exe/)
1. Clone this repo:  
`git clone https://github.com/docmanny/smRecSearch.git`
1. Create and install the conda environment needed for this pipeline:  
`conda env create --name RecSearch --file envs/conda_env.yaml`
1. Activate the environment using `conda activate RecSearch`.

Data needed prior to use
----

### Genomes (FASTA)
The core program underlying this pipeline, [`RecSearch`](https://github.com/docmanny/RecSearch), is flexible in its use of genome types and search algorithms. To reproduce this paper, you will need the genomes in the following table. Be sure to save them to the `data/genomes` folder.

|     Species                           |     Common Name                 |     Highest Quality Genome      | Citation/Link                                                                 |
|---------------------------------------|---------------------------------|---------------------------------|--------------------------------------------------------------------------|
|     Choloepus hoffmanni               |     Hoffmans two-toed sloth     | choHof-C_hoffmanni-2.0.1_HiC    | [DNAZoo](https://www.dnazoo.org/assemblies/Choloepus_hoffmanni)          |
|     Chrysochloris asiatica            |     Cape golden mole            | chrAsi1                         | [NCBI: chrAsi1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000296735.1/)  |
|     Dasypus novemcinctus              |     Nine-banded armadillo       | dasNov3                         | [NCBI: dasNov3](https://www.ncbi.nlm.nih.gov/assembly/GCA_000208655.2/)  |
|     Echinops telfairi                 |     Lesser Hedgehog Tenrec      | echTel2                         | [NCBI: echTel2](https://www.ncbi.nlm.nih.gov/assembly/GCA_000313985.1/)  |
|     Elephantulus edwardii             |     Cape elephant shrew         | eleEdw1                         | [NCBI: eleEdw1](https://www.ncbi.nlm.nih.gov/assembly/GCA_000299155.1/)  |
|     Elephas maximus                   |     Asian elephant              | eleMaxD                         | [Palkopoulou _et al._ 2018]     |
|     Loxodonta africana                |     African savanna elephant    |     loxAfr4                     | [ftp://ftp.broadinstitute.org/pub/assemblies/mammals/elephant/loxAfr4]() |
|     Loxodonta cyclotis                |     African forest elephant     |     loxCycF                     | [Palkopoulou _et al._ 2018]     |
|     Mammut americanum                 |     American mastodon           |     mamAmeI                     | [Palkopoulou _et al._ 2018]     |
|     Mammuthus columbi                 |     Columbian mammoth           |     mamColU                     | [Palkopoulou _et al._ 2018]     |
|     Mammuthus primigenius             |     Woolly mammoth              |     mamPriV                     | [Palkopoulou _et al._ 2015](https://doi.org/10.1016/j.cub.2015.04.007)   |
|     Orycteropus afer                  |     Aardvark                    |     oryAfe2                     | [DNAZoo](https://www.dnazoo.org/assemblies/Orycteropus_afer)             |
|     Palaeoloxodon antiquus            |     Straight tusked elephant    |     palAntN                     | [Palkopoulou _et al._ 2018]     |
|     Procavia capensis                 |     Rock hyrax                  |     proCap-Pcap_2.0_HiC         | [DNAZoo](https://www.dnazoo.org/assemblies/Procavia_capensis)            |
|     Trichechus manatus latirostris    |     Manatee                     |     triManLat2                  | [DNAZoo](https://www.dnazoo.org/assemblies/Trichechus_manatus)           |
| Homo sapiens  | Human | hg38 | [NCBI: hg38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39)|

### Port Table Configuration for BLAT/gfServer/gfClient
To use BLAT and gfServer/gfClient, you must specify a port for each gfServer. Additionally, various parts of this pipeline require interconverting between species names and genome assemblies. To facilitate this, a file has been included named `portTable.csv` in the `data` folder. Various genomes and species have already been included in this file along with suggested ports. If on your system, these ports are used by other processes, they can be changed without issue.

### Query Files
You should save all query sequences in `data/input`. Note that currently the pipeline assumes a FASTA input sequence with the extension ".fa", however, this can be easily changed in `rules/RecSearch.smk`. For convenience, the `AvA.fa` file containing the master list of sequences is included in `data/input`.

### SraRunTable

Our pipeline will generate _de novo_ transcriptomes for target genomes using SRA identifiers and the HISAT2-StringTie pipeline. A table of SRAs for _Loxodonta africana_, _Trichechus manatus_, and _Dasypus novemcinctus_ is included in `data/SraRunTable/SraRunTable.csv`.

### Other lines of evidence for Reciprocal Best-Hits (optional)
In addition to performing RBH Searches, this pipeline can intersect the hits with other lines of evidence to validate the results, and return a list of "evidenced" hits. To do so, simply download the other evidence as either a BED or a GFF file into either `data/BED` or `data/GFF`, respectively.


## Usage

I highly suggest reading both the [`snakemake` tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html) and looking through the different snakefiles in the `rules/` folder to familiarize yourself with the rules.

While it is possible to run the entire pipeline from start to end using `snakemake --use-conda publication/manuscript.pdf`, I would *strongly* recommend the following execution order:

1. Generate the modified `hg38` for reciprocal best hit searches using your query file of interest:  
`snakemake --use-conda output/recBlastDBPrep/hg38_maskRep_noVarChr_fragWithGenes.2bit`
1. Next, confirm that `RecSearch` runs correctly for one genome:  
`snakemake --use-conda -npr output/loxAfr4/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/RBB/loxAfr4_RecBlastOutput.bed.rbb`
1. Generate the transcriptomic evidence for functional duplicates:
`snakemake --use-conda -npr data/BED/loxAfr4-finalGuide.bed`

Once these steps are troubleshooted, it is possible to repeat the RecSearch and evidence steps individually for each genome; or continue to the next step and allow `snakemake` to generate them automatically.

4. Generate the table of genes per genome with copy number and ECNC scores (plus required files):  
`snakemake --use-conda -npr output/geneCopyTable/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/atlantogenata-GeneCopyTable_{RBB,evidenced}_filtered_long.csv`
5. Generate the maximum-likelihood tree for ancestral gene copy numbers:  
`snakemake --use-conda -npr  output/iqtree/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/maxLikelihood_model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/atlantogenata_RBB_filtered_dyn.iqtree`
6. Generate the final manuscript:  
`snakemake --use-conda -npr publication/manuscript.pdf`

Issues
------
**If you run into any issues with this pipeline**, or see that there are incomplete rules, please submit an Issue using a [reproducible example](https://stackoverflow.com/help/minimal-reproducible-example); including all error codes, log files, and any leads or investigative work that you did will go a long way towards making this manuscript and workflow even better, faster! While I can't guarantee that this will work on any computer, I can at least attest that it works well on Linux and Mac systems, in addition to SLURM clusters.  


[Palkopoulou _et al._ 2018]: https://doi.org/10.1073/pnas.1720554115
[1.2]: http://i.imgur.com/wWzX9uB.png
[mainmanmanny]: https://twitter.com/TheMainManManny
[devoevomed]: https://twitter.com/devoevomed?lang=en
