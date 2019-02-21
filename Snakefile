import pandas as pd


########################################################################################################################
######## Config
########################################################################################################################

configfile: "config.yaml"
localrules: all



pt = pd.read_csv(
              config["port_table"], 
              header=None, 
              names=["Translated_BLAT_Port", "Untranslated_BLAT_Port", "Genome", "Name", "Species"], 
              na_filter=True
          ).dropna().set_index("Species", drop=False)

## Output directories
dir_genome = "data/genome"
dir_2bit = "data/2bit"
dir_flags = "flags"
dir_idx = "output/hisat2-build"
dir_bam = "output/hisat2"
dir_sra = "data/SraRunTable"
dir_gff3 = "output/StringTie/deNovo"
dir_stringtie = "output/StringTie/final"

## Log directories
log_idx = "logs/hisat2-build"
log_bam = "logs/hisat2"
log_bai = "logs/samtools-index"
log_gff3 = "logs/StringTie/deNovo"
log_stringtie = "logs/StringTie/final"
log_flags = "logs/flags"


genomes = pt["Genome"]
species = pt.index

genome_regex = "{genome, ((([A-Za-z]{1,3})([A-Za-z]{3})+)((\d{1,2}m?|[A-Z])(\.\d+)?(\.softmask)?)?(\d\.[a-z]{3}\.[a-z]{3}\.\d{8})?(-[A-Z][A-Za-z0-9]*([_.][A-Za-z]+){0,2})?(_\d{2}[A-Z][a-z]{2}\d{4}_[A-Za-z0-9]{0,6})?([-_][Vv]?\d{1,2}(\.\d{1,2}){0,2}k?)?(_HiC)?)|(mm10)|(rn6)|(hg19)|(hg38(_maskRep_noVarChr_fragWithGenes)?)|(dm6)|(fr2)}"

SRA_regex = "{sra, [SED]RR\d+}"

########################################################################################################################
######## Rules
########################################################################################################################

rule all:
    input:
        "final thing"

rule test:
    input:
        genomes=expand("/usr/db/{g}", g=pt["Genome"])
    shell:
        "echo {genomes}"


########################################################################################################################
######## Imports
########################################################################################################################

# Rules for assembling transcriptome
include: "rules/RNA-seq.smk"
# Rules for RecBlat
include: "rules/RecBLAT.smk"
# Rules for intersecting the RecBlat output
include: "rules/intersections.smk"
