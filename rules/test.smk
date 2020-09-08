import pandas as pd


########################################################################################################################
######## Config
########################################################################################################################

def get_SRA_list(sra_file):
    """Reads a .tsv SRA file and returns a dataframe indexed by SRA Run accessions"""
    return pd.read_table(sra_file).set_index("Run", drop=False)


def get_SRA_acc(species, dir_sra, suffix="_SraRunTable.tsv"):
    """Returns a list of SRA Run Accession numbers from a SraRunTable"""
    try:
        species = species.replace(" ", "_")
        sra_table = species + suffix
        sra_file = dir_sra.rstrip("/") + "/" + sra_table
        df = get_SRA_list(sra_file)
        return list(df.index)
    except FileNotFoundError:
        return ""
    except AttributeError:
        return ""


def generate_genome_SRA_list(genome_df, dir_sra, suffix="_SraRunTable.tsv"):
    """Generates a set of species-SRA strings from a list of species"""
    species_sra="{species}_{sra}"
    final_list = []
    for species, genome in zip(genome_df.index, genome_df):
        sra_list = get_SRA_acc(species, dir_sra, suffix)
        if sra_list:
            final_list.extend(map(lambda sra: (genome, sra), sra_list))
        else:
            continue
    return [species_sra.format(species=i[0], sra=i[1]) for i in final_list]


########################################################################################################################
######## Config
########################################################################################################################

## Output directories
dir_genome = "data/genome"
dir_idx = "output/hisat2-build"
dir_bam = "output/hisat2"
dir_sra = "data/SraRunTable"
dir_stringtie = "output/StringTie"

## Log directories
log_idx = "log/hisat2-build/"
log_bam = "log/hisat2/"
log_bai = "log/samtools-index"
log_stringtie = "log/StringTie"

genomes = pt["Genome"]
species = pt.index

########################################################################################################################
######## Targets
########################################################################################################################
#ALL_BAM = expand("{path}/{sample}/{sample}.bam", path = BAM_dir, sample = ALL_SAMPLES)
#ALL_EXP = expand("{path}/{sample}/{sample}.exp", path = BAM_dir, sample = ALL_SAMPLES)
#ALL_BAMSTAT = expand("{path}/{sample}/{sample}.sam_stat", path = BAM_dir, sample = ALL_SAMPLES)
#ALL_TIN = expand("{path}/{sample}/{sample}.summary.txt", path = BAM_dir, sample = ALL_SAMPLES)
#ALL_COVERAGE = expand("{path}/{sample}/{sample}.geneBodyCoverage.r", path = BAM_dir, sample = ALL_SAMPLES)
#ALL_HTSEQCOUNT = expand(config["BAM_dir"] + "/{sample}/{sample}.htseq_count", sample = ALL_SAMPLES)

#rule all:
#     input: ALL_BAM + ALL_EXP + ALL_BAMSTAT + ALL_TIN + ALL_COVERAGE
#    input: ALL_BAM + ALL_EXP + ALL_BAMSTAT + ALL_TIN + ALL_COVERAGE + ALL_LOG + ALL_GTF + ALL_HTSEQCOUNT

localrules: all


########################################################################################################################
######## Rules
########################################################################################################################

rule all:
    input:
        expand("{path}/{genome_sra}.bai",
               path = dir_bam,
               genome_sra = generate_genome_SRA_list(genomes, dir_sra=dir_sra)
               )

rule hisat2_index:
    input:
        expand("{path}/{{genome}}.fa", path=dir_genome)
    output:
        expand("{path}/{{genome}}", path=dir_idx)
#    params:
#        hi=expand("{path}/{g}", path=dir_idx, g=pt["Genome"])
    threads: 10
    log:
        expand("{path}/{{genome}}", path=log_idx)
    shell:
        "hisat2-build -p {threads} {input} {output}"


rule hisat2:
    input:
        idx = expand("{path}/{{genome}}", path=dir_idx)
    output:
        bam = expand("{path}/{{genome}}_{{sra}}.bam", path=dir_bam)
    log:
        expand("{path}/{{genome}}_{{sra}}.log", path = log_bam)
    params:
        summary = expand("{path}/{{genome}}_{{sra}}.summary", path = log_bam)
    threads: 10
    shell:
        "hisat2 -p {threads} "
        "-x {input.idx} "
        "--dta "
        "--summary-file {params.summary} "
        "--met-stderr "
        "--sra-acc {wildcards.sra} "
        "2> {log} | "
        "samtools view -b -u | "
        "samtools sort - "
        "--threads {threads} "
        "-o {output}"

rule samtools_index:
    input:
        bam = expand("{path}/{{genome}}_{{sra}}.bam", path=dir_bam)
    output:
        bai = expand("{path}/{{genome}}_{{sra}}.bai", path=dir_bam)
    log: expand("{path}/{{genome}}_{{sra}}.log", path=log_bai)
    threads: 10
    shell:
        "samtools index {input.bam}"
