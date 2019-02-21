import pandas as pd


########################################################################################################################
######## Config
########################################################################################################################

genome_regex = "{genome, ((([A-Za-z]{1,3})([A-Za-z]{3})+)((\d{1,2}m?|[A-Z])(\.\d+)?(\.softmask)?)?(\d\.[a-z]{3}\.[a-z]{3}\.\d{8})?(-[A-Z][A-Za-z0-9]*([_.][A-Za-z]+){0,2})?(_\d{2}[A-Z][a-z]{2}\d{4}_[A-Za-z0-9]{0,6})?([-_][Vv]?\d{1,2}(\.\d{1,2}){0,2}k?)?(_HiC)?)|(mm10)|(rn6)|(hg19)|(hg38(_maskRep_noVarChr_fragWithGenes)?)|(dm6)|(fr2)}"

SRA_regex = "{sra, [SED]RR\d+}"

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
    genome_sra="{genome}_{sra}"
    final_list = []
    for species, genome in zip(genome_df.index, genome_df):
        sra_list = get_SRA_acc(species, dir_sra, suffix)
        if sra_list:
            final_list.extend(map(lambda sra: (genome, sra), sra_list))
        else:
            continue
    return [genome_sra.format(genome=i[0], sra=i[1]) for i in final_list]

def sra_from_wildcards(wildcards, genome_df = pt, suffix="_SraRunTable.tsv"):
    return get_SRA_acc(genome_df.loc[genome_df['Genome'] == wildcards.genome].index[0], dir_sra, suffix)

def stmerge_inputs(wildcards, genome_df = pt, suffix="_SraRunTable.tsv"):
    sras = sra_from_wildcards(wildcards, genome_df = pt, suffix="_SraRunTable.tsv")
    return ["{path}/{genome}_{sra}.gff3".format(path = dir_gff3,
                                                genome=wildcards.genome,
                                                sra=s)
                for s in sras
           ]

########################################################################################################################
######## Rules
########################################################################################################################

rule all_alignments:
    input:
        expand("data/BED/{genome_sra}-final.gff3",
               genome_sra = generate_genome_SRA_list(genomes, dir_sra=dir_sra)
               )

rule hisat2_index:
    input:
        "{path}/{genome}{extension}".format(genome = genome_regex, path=dir_genome, extension = ".fa")
    output:
        protected("{path}/{genome}{extension}".format(genome = genome_regex, path=dir_idx, extension = ""))
#    params:
#        hi=expand("{path}/{g}", path=dir_idx, g=pt["Genome"])
    threads: 10
    log: "{path}/{genome}{extension}".format(genome = genome_regex, path=log_idx, extension = ".idx.log")
    shell:
        "hisat2-build -p {threads} {input} {output}"

rule hisat2:
    input:
        idx = "{path}/{genome}{extension}".format(genome = genome_regex, path=dir_idx, extension = "")
    output:
        bam = temp("{path}/{genome}_{sra}{extension}".format(genome = genome_regex, sra = SRA_regex, path=dir_bam, extension = ".bam"))
        #bam = temp(expand("{path}/{{genome}}_{{sra}}.bam", path=dir_bam))
    log:
        "{path}/{genome}_{sra}{extension}".format(genome = genome_regex, sra = SRA_regex, path=log_bam, extension = ".bam.log")
    params:
        summary = "{path}/{genome}_{sra}{extension}".format(genome = genome_regex, sra = SRA_regex, path=log_bam, extension = ".bam.summary")
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
        bam = "{path}/{genome}_{sra}{extension}".format(genome = genome_regex, sra = SRA_regex, path=dir_bam, extension = ".bam")
    output:
        bai = temp("{path}/{genome}_{sra}{extension}".format(genome = genome_regex, sra = SRA_regex, path=dir_bam, extension = ".bai"))
    log: "{path}/{genome}_{sra}{extension}".format(genome = genome_regex, sra = SRA_regex, path=log_bam, extension = ".bai.log")
    threads: 10
    shell:
        "samtools index {input.bam}"

rule stringtie_denovo:
    input:
        bam = "{path}/{genome}_{sra}{extension}".format(genome = genome_regex, sra = SRA_regex, path=dir_bam, extension = ".bam")
    output:
        gff3 = temp("{path}/{genome}_{sra}{extension}".format(genome = genome_regex, sra = SRA_regex, path=dir_gff3, extension = ".gff3"))
    log: "{path}/{genome}_{sra}{extension}".format(genome = genome_regex, sra = SRA_regex, path=log_gff3, extension = ".gff3.log")
    params:
    threads: 10
    shell:
        "stringtie {input.bam} "
        "-p {threads} "
        "-l {wildcards.genome}_{wildcards.sra} "
        "-v "
        "-o {output.gff3} "
        "2> {log}"

rule stringtie_merge:
    input: stmerge_inputs
    output: "{path}/{genome}{extension}".format(genome = genome_regex, path=dir_stringtie, extension = "-guide.gff3")
    log: "{path}/{genome}{extension}".format(genome = genome_regex, path=log_stringtie, extension = "-guide.gff3.log")
    threads: 10
    shell:
        "stringtie --merge "
        "{input}"
        "-l {wildcards.genome} "
        "-v "
        "-o {output} "
        "2> {log}"

rule stringtie_final:
    input:
        bam = "{path}/{genome}_{sra}{extension}".format(genome = genome_regex, sra = SRA_regex, path=dir_bam, extension = ".bam"),
        guide = "{path}/{genome}{extension}".format(genome = genome_regex, path=dir_stringtie, extension = "-guide.gff3")
    output:
        gff3 = protected("{path}/{genome}_{sra}{extension}".format(genome = genome_regex, sra = SRA_regex, path=dir_stringtie, extension = "-final.gff3"))
    log: "{path}/{genome}_{sra}{extension}".format(genome = genome_regex, sra = SRA_regex, path=log_stringtie, extension = "-final.gff3.log")
    threads: 10
    shell:
        "stringtie {input.bam} "
        "-p {threads} "
        "-G {input.guide} "
        "-l {wildcards.genome} "
        "-vBAC "
        "-o {output.gff3} "
        "2> {log}"

rule GFF3ToBed:
    input:
        "{path}/{{genome}}_{{sra}}{extension}".format(path=dir_stringtie, extension = "-final.gff3")
    output:
        "data/BED/{genome}_{sra}.bed".format(genome = genome_regex, sra = SRA_regex)
    shell:
        "gff2bed < {input} > {output}"