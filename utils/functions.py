from typing import List
def query(d:List[dict], k:str, v:str) -> dict:
    """Return the first dictionary in a list of dictionaries where the value of key k matches v."""
    return [x for x in d if x[k] == v][0]

def query_all(d:List[dict], k:str, v:str, k_out:str) -> List[str]:
    """Return a list of values from key k_out in a list of dictionaries where the value of key k matches v."""
    return [x[k_out] for x in d if x[k] == v]

def _genome_assembly(wildcards):
    """Get genome assembly file path by sample name"""
    species = query(
        config["samples"] + \
        config["samples-dna"] + \
        config["samples-lcm"], 
        "name", 
        wildcards.sample
    )["species"]
    return query(config["references"], "species", species)["assembly"]
