#!/usr/bin/env python3

import click
from pathlib import Path
from RecBlast.RecBlast import RecSearch

@click.command()
@click.option("-q", "--query-file", type=click.Path(exists=True))
@click.option("--query-file-type", type=str, default="fasta")
@click.option("-p", "--max-processes", type=int, default=40)
@click.option("-fp", "--forward-port")
@click.option("-rp", "--reverse-port")
@click.option("-fs", "--forward-species", type=str)
@click.option("-ft", "--forward-twobit", type=click.Path(exists=False))
@click.option("-rs", "--reverse-species", type=str)
@click.option("-rt", "--reverse-twobit", type=click.Path(exists=False))
@click.option("-ps", "--perc-score", type=str)
@click.option("-pi", "--perc-identity", type=str)
@click.option("-pq", "--perc-query-span", type=str)
@click.option("--annotation_lookup_tsv", type=str)
@click.option("--output-root", type=str, default="./output")
def __main__(query_file, forward_port, forward_species, forward_twobit,
             reverse_port, reverse_species, reverse_twobit, query_file_type="fasta",
             perc_score=0.1, perc_identity=0.5, perc_query_span=0.5, max_processes=40,
             annotation_lookup_tsv="./output/recBlastDBPrep/hg38_geneAcc_hashTable.tsv",
             output_root = "./output"):
    perc_score = float(perc_score)
    perc_identity = float(perc_identity)
    perc_query_span = float(perc_query_span)

    forward_twobit = Path(forward_twobit)
    reverse_twobit = Path(reverse_twobit)

    #print(forward_twobit, reverse_twobit, output_root, perc_identity, perc_score, perc_query_span, query_file, sep="\n")
    output_location = Path(output_root, forward_twobit.stem)
    print(output_location)

    recblast = RecSearch(target_species=forward_species, query_species=reverse_species,
                         forward_search_type="tblat", reverse_search_type="blat",
                         sequence_source="twobit", verbose=2)
    recblast.max_processes = max_processes
    recblast.set_queries(query_file,
                         infile_type=query_file_type)
    #print(recblast.records)
    recblast.forward_search_settings['database_port'] = {forward_species: forward_port}
    recblast.forward_search_settings['database'] = {forward_species: str(forward_twobit.name)}
    recblast.forward_search_settings['database_path'] = str(forward_twobit.parent)
    recblast.forward_search_criteria = dict(perc_score=perc_score,
                                            perc_ident=perc_identity,
                                            perc_query_span=perc_query_span)
    recblast.sequence_source_settings['database'] = {forward_species: str(forward_twobit.name)}
    recblast.sequence_source_settings['database_path'] = str(forward_twobit.parent)
    recblast.memory_saver_level = 1
    recblast.reverse_search_settings['database'] = {reverse_species: str(reverse_twobit.name)}
    recblast.reverse_search_settings['database_path'] = str(reverse_twobit.parent)
    recblast.reverse_search_settings['database_port'] = {reverse_species: reverse_port}
    recblast.set_translation_annotation_parameters(method="table", key_value_order=False,
                                                   tsv_location=annotation_lookup_tsv)
    recblast(run_name="{0}-pcScore{1}_pcIdent{2}_pcQuerySpan{3}_reverse-{4}".format(Path(query_file).stem, perc_score, perc_identity, perc_query_span, reverse_twobit.stem),
             output_type="bed-complete",
             output_location=output_location)
    

if __name__ == "__main__":
    __main__()

exit()
