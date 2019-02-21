import pandas as pd

def get_port_for_genome(df, genome, translated = True):
    if translated:
        return int(df.loc[df["Genome"]==genome].Translated_BLAT_Port)
    else:
        return int(df.loc[df["Genome"]==genome].Untranslated_BLAT_Port)

def get_species_port_translated(df):
    return ["{{path}}/translated-{0}-{1}".format(i, j) for i, j in zip(df.Translated_BLAT_Port, df.Genome)]

def get_species_port_untranslated(df):
    return ["{{path}}/untranslated-{0}-{1}".format(i, j) for i, j in zip(df.Untranslated_BLAT_Port, df.Genome)]

def get_species_port_closed(df):
    l = []
    l += ["{{path}}/translated-{0}-{1}-closed".format(j, i) for i, j in zip(df.Genome, df.Translated_BLAT_Port)]
    l += ["{{path}}/untranslated-{0}-{1}-closed".format(j, i) for i, j in zip(df.Genome, df.Untranslated_BLAT_Port)]
    return l

def get_port(f):
    with open(f) as infile:
        return int(infile.read())

def filtered_beds(wildcards):
    d = "output/{genome}/AvA-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/filtered/{genome}_intersected_{{type}}.bed".format(
           genome = wildcards.genome, pc_score = wildcards.pc_score, pc_ident=wildcards.pc_ident, pc_qspan = wildcards.pc_qspan, rgenome = wildcards.rgenome)
    t = glob_wildcards(d)
    return expand(d, type=t.type)

# Parameters

# Output Directories


# Log Directories

rule RBHB:
    input:
        expand("output/{genome}/AvA-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/filtered/{genome}_evidenced.bed",
               genome="loxAfr3", pc_score = 0.1, pc_ident = 0.8, pc_qspan = 0.5, rgenome = "hg38_maskRep_noVarChr_fragWithGenes")

rule start_translated_gfServer:
    input:
        twoBitFile="{path}/{{genome}}{extension}".format(path=dir_2bit, extension = ".2bit")
    output:
        "{path}/translated-{genome}".format(path=dir_flags, genome = genome_regex)
    log:
        "{path}/translated-{genome}.log".format(path=log_flags, genome = genome_regex)
    shell:
        "./code/start_Translated_gfServer.sh {wildcards.genome} > {output} 2>{log}"

rule start_untranslated_gfServer:
    input:
        twoBitFile="{path}/{{genome}}{extension}".format(path=dir_2bit, extension = ".2bit")
    output:
        "{path}/untranslated-{genome}".format(path=dir_flags, genome = genome_regex)
    log:
        "{path}/untranslated-{genome}.log".format(path=log_flags, genome = genome_regex)
    shell:
        "./code/start_Untranslated_gfServer.sh {wildcards.genome} > {output} 2>{log}"

rule RecBlat:
    input:
        qfile="data/input/CheAbi_NEIL1_Proteins.txt", #"data/input/hg_UP000005640_querylist.fasta",
        twoBitFile="{path}/{{genome}}{extension}".format(path=dir_2bit, extension = ".2bit"),
        tportfile = "{path}/translated-{{genome}}".format(path=dir_flags),
        utportfile = "{path}/untranslated-{{rgenome}}".format(path=dir_flags),
        annoTable = "output/recBlastDBPrep/hg38_geneAcc_hashTable.tsv"
    output:
        ancient(protected("output/{genome}/AvA-pcScore{{pc_score}}_pcIdent{{pc_ident}}_pcQuerySpan{{pc_qspan}}_reverse-{{rgenome}}/output/{{genome}}_RecBlastOutput.bed".format(genome = genome_regex, rgenome = genome_regex.replace("genome", "rgenome"))))
    params:
        tport = lambda wildcards, df=pt: get_port_for_genome(df, wildcards.genome, translated=True),
        utport = lambda wildcards, df=pt: get_port_for_genome(df, wildcards.rgenome, translated=False)
    threads: 10
    shell:
        "./code/rbb.py "
        "--query-file {input.qfile} "
        "--forward-port {params.tport} "
        "--reverse-port {params.utport} "
        "--forward-species {wildcards.genome} "
        "--forward-twobit {input.twoBitFile} "
        "--reverse-species 'Homo sapiens' "
        "--reverse-twobit 'hg38_maskRep_noVarChr_fragWithGenes.2bit' "
        "--annotation_lookup_tsv {input.annoTable} "
        "--perc-score {wildcards.pc_score} "
        "--perc-identity {wildcards.pc_ident} "
        "--perc-query-span {wildcards.pc_qspan} "
        "--max-processes {threads} "
        #"gfServer stop localhost {params.tport} && "
        #"rm {input.tportfile}"

rule reciprocalBestHits:
    input:
        "output/{genome}/AvA-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/output/{genome}_RecBlastOutput.bed"
    output:
        "output/{genome}/AvA-pcScore{{pc_score}}_pcIdent{{pc_ident}}_pcQuerySpan{{pc_qspan}}_reverse-{{rgenome}}/RBB/{{genome}}_RecBlastOutput.bed.rbb".format(genome = genome_regex, rgenome = genome_regex.replace("genome", "rgenome"))
    script:
        "../code/rbhb_from_bed.py"

rule filterRBHBbyBED:
    input:
        a="output/{genome}/AvA-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/RBB/{genome}_RecBlastOutput.bed.rbb",
        b="data/BED/{genome}_{type}.bed"
    output:
        "output/{genome}/AvA-pcScore{{pc_score}}_pcIdent{{pc_ident}}_pcQuerySpan{{pc_qspan}}_reverse-{{rgenome}}/filtered/{{genome}}_intersected_{{type}}.bed".format(genome = genome_regex, rgenome = genome_regex.replace("genome", "rgenome"))
    log:
        "logs/{genome}/AvA-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/filtered/{genome}_intersected_{type}.log"
    shell:
        "bedtools intersect -u -a {input.a} -b {input.b} | bedtools sort > {output} 2>{log}"

rule evidencedRBHB:
    input:
        filtered_beds
    output:
        "output/{genome}/AvA-pcScore{{pc_score}}_pcIdent{{pc_ident}}_pcQuerySpan{{pc_qspan}}_reverse-{{rgenome}}/filtered/{{genome}}_evidenced.bed".format(genome = genome_regex, rgenome = genome_regex.replace("genome", "rgenome"))
    #log:
    #    "logs/{genome}/AvA-pcScore{pc_score}_pcIdent{pc_ident}_pcQuerySpan{pc_qspan}_reverse-{rgenome}/filtered/{genome}_intersected_{type}.log"
    shell:
        #"cat {input} | sort | uniq | bedtools sort > {output}"
        "echo {input}"