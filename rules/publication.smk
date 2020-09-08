rule ORATable:
    input:
        directory("output/ORA/maxLikelihood_model_MK+FQ+I+G4-dataType_MORPH-asrMin_0.8/AvA-pcScore0.1_pcIdent0.8_pcQuerySpan0.5_reverse-hg38_maskRep_noVarChr_fragWithGenes/atlantogenata_RBB_filtered_dyn/pathway_Reactome-top100/")
    output:
         "../output/pubFiles/ORA_table_pub.csv"
    # conda: "../envs/R-geneLists.yaml"
    script: "../code/pub_ORATable.R"


rule bodySizeTree:
    input:
        bsize = "output/AncSizeRec/AnnotatedTreeForFigure.tree",
        nexus = "data/stableTraits/eutheria.nexus",
        ancstates = "data/stableTraits/eutheria.ancstates",
        brlens = "data/stableTraits/eutheria.brlens",
        genomeTable = "output/other/genomeTable.csv"
    output:
        outTree_full = "output/pubFiles/bodysize-Eutheria.csv",
        outTree_Atlantogenata = "output/pubFiles/bodysize-Atlantogenata.csv"
    script: "../code/pub_bodySizeTree.R"


rule lifespanPGLS_cancerRisk:
    input:
        anage = "data/AnAge/anage_build14.txt",
        Atlantogenata = "output/pubFiles/bodysize-Atlantogenata.csv"
    output:
        cancerSuccep = "output/pubFiles/cancerRisk-Atlantogenata.csv",
        latex = "output/pubFiles/lifespan-PGLS.tex",
        html = "output/pubFiles/lifespan-PGLS.html"
    script: "../code/pub_lifespanPGLS_cancerRisk.R"









# rule knit_paper:
#     input:
#         RMarkdown = "paper_PLOS/paper_PLOS_draft.Rmd",
#         UP000005640 = "data/input/UP000005640.withheader.tsv",
#         gene_list = "output/recBlastDBPrep/AvA_geneList.txt",
#         gene_list_filtered = "output/recBlastDBPrep/AvA_geneList_filtered.txt",
#         genomeTable = "output/other/genomeTable.csv",
#         stabletraits_progress = "data/stableTraits/eutheria.progress",
#         time_tree = "data/stableTraits/eutheria.tree"
#     output: "paper_PLOS/paper_PLOS_draft.html"
#     script:
#         "code/renderpaper.R"
