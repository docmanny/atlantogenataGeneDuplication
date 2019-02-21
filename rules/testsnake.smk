import pandas as pd
from snakemake.utils import validate


def get_SRA_list(sra_file):
    """Reads a .tsv SRA file and returns a dataframe indexed by SRA Run accessions"""
    return pd.read_table(sra_file).set_index("Run", drop=False)


def get_SRA_acc(sra_table):
    """Returns a list of SRA Run Accession numbers from a SraRunTable"""
    try:
        df = get_SRA_list(sra_table)
        return list(df.index)
    except FileNotFoundError:
        return ""
    except AttributeError:
        return ""


def generate_species_SRA_list(species_list, delimiter="."):
    """Generates a set of species-SRA strings from a list of species"""
    species_sra="{species}.{sra}"
    final_list = []
    for species in species_list:
        sra_list = get_SRA_acc(species)
        if sra_list:
            final_list.extend(map(lambda sra: (species, sra), sra_list))
        else:
            continue
    return [species_sra.format(species=i[0], sra=i[1]) for i in final_list]

## Config file just has the line "port_table: './data/portTable.csv'"
#configfile: "config.yaml"
## I'll just put that here instead for this example
port_table="./data/portTable.csv"

samples = pd.read_csv(
              config["port_table"], 
              header=None, 
              names=["Translated_BLAT_Port", "Untranslated_BLAT_Port", "twoBit", "Common_Name", "Species"], 
              na_filter=True
          ).dropna()
samples["Species"] = samples["Species"].str.replace(' ', '_')
samples = samples.set_index("Species", drop=False)

        
rule test:
    output:
        generate_species_SRA_list(set(samples.index))
    shell:
        "touch {output}"