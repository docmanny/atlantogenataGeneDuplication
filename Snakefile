import pandas as pd


########################################################################################################################
######## Config
########################################################################################################################

configfile: "config.yaml"
localrules: all, linkBed, linkGeneLists, linkCopyTable, linkFinalBed, WebGestalt_oneShot, WebGestalt_intersection, WebGestalt_duplicates, WebGestalt_ML


wildcard_constraints:
    genome="((([A-Za-z]{1,3})([A-Za-z]{3})+)((\d{1,2}m?|[A-Z])(\.\d+)?(\.softmask)?)?(\d\.[a-z]{3}\.[a-z]{3}\.\d{8})?(-[A-Z][A-Za-z0-9]*([_.][A-Za-z]+){0,2})?(_\d{2}[A-Z][a-z]{2}\d{4}_[A-Za-z0-9]{0,6})?([-_][Vv]?\d{1,2}(\.\d{1,2}){0,2}k?)?(_HiC)?)|(mm10)|(rn6)|(hg19)|(hg38(_maskRep_noVarChr_fragWithGenes)?)|(dm6)|(fr2)",
    rgenome="((([A-Za-z]{1,3})([A-Za-z]{3})+)((\d{1,2}m?|[A-Z])(\.\d+)?(\.softmask)?)?(\d\.[a-z]{3}\.[a-z]{3}\.\d{8})?(-[A-Z][A-Za-z0-9]*([_.][A-Za-z]+){0,2})?(_\d{2}[A-Z][a-z]{2}\d{4}_[A-Za-z0-9]{0,6})?([-_][Vv]?\d{1,2}(\.\d{1,2}){0,2}k?)?(_HiC)?)|(mm10)|(rn6)|(hg19)|(hg38(_maskRep_noVarChr_fragWithGenes)?)|(dm6)|(fr2)"


pt = pd.read_csv(
              config["port_table"], 
              header=None, 
              names=["Translated_BLAT_Port", "Untranslated_BLAT_Port", "Genome", "Name", "Species"], 
              na_filter=True,
              comment='#'
          ).dropna().set_index("Species", drop=False)

## Output directories
dir_genome = "data/genome"
dir_2bit = "data/2bit"
dir_flags = "flags"
dir_idx = "output/hisat2-build"
dir_bam = "output/hisat2"
dir_sra = "data/SraRunTable"
dir_gff3 = "output/StringTie/deNovo"
dir_stringtieMerge = "output/StringTie/merge"
dir_stringtie = "output/StringTie/final"
dir_stringtieMergeFinal = "output/StringTie/finalMerge"

## Log directories
log_idx = "logs/hisat2-build"
log_bam = "logs/hisat2"
log_bai = "logs/samtools-index"
log_gff3 = "logs/StringTie/deNovo"
log_stringtie = "logs/StringTie/final"
log_stringtieMergeFinal = "logs/StringTie/finalMerge"
log_flags = "logs/flags"


genomes = pt["Genome"]
species = pt.index

genome_regex = "{genome, ((([A-Za-z]{1,3})([A-Za-z]{3})+)((\d{1,2}m?|[A-Z])(\.\d+)?(\.softmask)?)?(\d\.[a-z]{3}\.[a-z]{3}\.\d{8})?(-[A-Z][A-Za-z0-9]*([_.][A-Za-z]+){0,2})?(_\d{2}[A-Z][a-z]{2}\d{4}_[A-Za-z0-9]{0,6})?([-_][Vv]?\d{1,2}(\.\d{1,2}){0,2}k?)?(_HiC)?)|(mm10)|(rn6)|(hg19)|(hg38(_maskRep_noVarChr_fragWithGenes)?)|(dm6)|(fr2)}"

SRA_regex = "{sra, [SED]RR\d+}"

wildcard_constraints:
    genome="((([A-Za-z]{1,3})([A-Za-z]{3})+)((\d{1,2}m?|[A-Z])(\.\d+)?(\.softmask)?)?(\d\.[a-z]{3}\.[a-z]{3}\.\d{8})?(-[A-Z][A-Za-z0-9]*([_.][A-Za-z]+){0,2})?(_\d{2}[A-Z][a-z]{2}\d{4}_[A-Za-z0-9]{0,6})?([-_][Vv]?\d{1,2}(\.\d{1,2}){0,2}k?)?(_HiC)?)|(mm10)|(rn6)|(hg19)|(hg38(_maskRep_noVarChr_fragWithGenes)?)|(dm6)|(fr2)",
    genomeA="((([A-Za-z]{1,3})([A-Za-z]{3})+)((\d{1,2}m?|[A-Z])(\.\d+)?(\.softmask)?)?(\d\.[a-z]{3}\.[a-z]{3}\.\d{8})?(-[A-Z][A-Za-z0-9]*([_.][A-Za-z]+){0,2})?(_\d{2}[A-Z][a-z]{2}\d{4}_[A-Za-z0-9]{0,6})?([-_][Vv]?\d{1,2}(\.\d{1,2}){0,2}k?)?(_HiC)?)|(mm10)|(rn6)|(hg19)|(hg38(_maskRep_noVarChr_fragWithGenes)?)|(dm6)|(fr2)",
    genomeB="((([A-Za-z]{1,3})([A-Za-z]{3})+)((\d{1,2}m?|[A-Z])(\.\d+)?(\.softmask)?)?(\d\.[a-z]{3}\.[a-z]{3}\.\d{8})?(-[A-Z][A-Za-z0-9]*([_.][A-Za-z]+){0,2})?(_\d{2}[A-Z][a-z]{2}\d{4}_[A-Za-z0-9]{0,6})?([-_][Vv]?\d{1,2}(\.\d{1,2}){0,2}k?)?(_HiC)?)|(mm10)|(rn6)|(hg19)|(hg38(_maskRep_noVarChr_fragWithGenes)?)|(dm6)|(fr2)",
    sra="[SED]RR\d+",
    hit_type = "(RBB)|(evidenced)",
    enrichmentDB = "pathway_KEGG|pathway_Panther|pathway_Reactome|pathway_Wikipathway|pathway_Wikipathway_cancer"
    

########################################################################################################################
######## Rules
########################################################################################################################

#rule all:
#    input:
#        "final thing"

rule test:
    input:
        genomes=expand("{dir}/{g}", g=pt["Genome"], dir=dir_genome)
    shell:
        "echo {genomes}"


########################################################################################################################
######## Imports
########################################################################################################################

# Rules for assembling transcriptome
include: "rules/RNA-seq.smk"
# Rules for RecBlat
include: "rules/RecBLAT.smk"
# Rules for using HAL and progressiveCactus
include: "rules/halTools.smk"
# Rules for intersecting the RecBlat output
include: "rules/intersections.smk"
# Rules for other things
include: "rules/extras.smk"
# Rules for WebGestalt
include: "rules/webgestalt.smk"


# rule report:
#     input:
#         RBB="output/{genome}/{query}-pcScore{pcS}_pcIdent{pcI}_pcQuerySpan{pcQS}_reverse-{rgenome}/RBB/{genome}_RecBlastOutput.bed.rbb",
#         Evidenced="output/{genome}/{query}-pcScore{pcS}_pcIdent{pcI}_pcQuerySpan{pcQS}_reverse-{rgenome}/filtered/{genome}_evidenced.bed",
#         ECNC="output/{genome}/{query}-pcScore{pcS}_pcIdent{pcI}_pcQuerySpan{pcQS}_reverse-{rgenome}/ECNC/{genome}_RecBlastOutput.bed.rbb.ecnc"
#     output: "docs/{genome}-{query}-pcScore{pcS}_pcIdent{pcI}_pcQuerySpan{pcQS}_reverse-{rgenome}-Report.html"
#     script: "analysis/test.Rmd"
#     
# 
# rule report_test:
#     input:
#         RBB="output/loxAfr3/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/RBB/loxAfr3_RecBlastOutput.bed.rbb",
#         Evidenced="output/loxAfr3/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/filtered/loxAfr3_evidenced.bed",
#         ECNC="output/loxAfr3/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/ECNC/loxAfr3_RecBlastOutput.bed.rbb.ecnc"
#     output:
#         "docs/test-report.html"
#     script:
#         "analysis/test2.Rmd"

rule knit_paper:
    input:
        RMarkdown = "paper_PLOS/paper_PLOS_draft.Rmd",
        UP000005640 = "data/input/UP000005640.withheader.tsv",
        gene_list = "output/recBlastDBPrep/AvA_geneList.txt",
        gene_list_filtered = "output/recBlastDBPrep/AvA_geneList_filtered.txt",
        genomeTable = "output/other/genomeTable.csv",
        stabletraits_progress = "data/stableTraits/eutheria.progress",
        time_tree = "data/stableTraits/eutheria.tree"
    output: "paper_PLOS/paper_PLOS_draft.html"
    script:
        "code/renderpaper.R"
        
