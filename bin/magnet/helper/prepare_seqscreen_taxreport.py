import os
from collections import defaultdict
import pandas as pd
import argparse
import pathlib
import warnings
warnings.filterwarnings("ignore")

from ete3 import NCBITaxa

def get_species_taxid(taxid, ncbi_taxa_db, valid_kingdom, ret_subspecies=False):
    lineage = ncbi_taxa_db.get_lineage(taxid)
    if bool(set(lineage) & valid_kingdom):
        taxid2rank_dict = ncbi_taxa_db.get_rank(lineage)
        species_rank_visited = False
        for lineage_taxid in lineage:
            if species_rank_visited and ret_subspecies:
                return lineage_taxid
            if taxid2rank_dict[lineage_taxid] == 'species':
                species_rank_visited = True
                if not ret_subspecies or lineage_taxid == lineage[-1]:
                    return lineage_taxid
    return None

def filter_seqscreen_taxonomy(seqscreen_output, ncbi_taxa_db, call_subspecies, valid_kingdom={2, 4751, 2157, 10239}):
    classification_result_df = pd.read_csv(os.path.join(seqscreen_output, 
                                                        'taxonomic_identification', 
                                                        'taxonomic_assignment', 
                                                        'taxonomic_results.txt'), 
                                           sep='\t', usecols=[0,1,2])
    total_read_count, _ = classification_result_df.shape
    
    taxid_count_dict = defaultdict(int)
    taxid_species_lookup = dict()
    error_count = 0
    for taxid_list in classification_result_df['taxid'].values:
        for taxid in taxid_list.split(','):
            try:
                taxid = int(taxid)
                try:
                    species_taxid = taxid_species_lookup[taxid]
                except KeyError:
                    species_taxid = get_species_taxid(taxid, ncbi_taxa_db, valid_kingdom, call_subspecies)
                    taxid_species_lookup[taxid] = species_taxid

                if species_taxid is not None:
                    taxid_count_dict[species_taxid] += 1
            except ValueError:
                error_count += 1

    taxid_queries = []
    for taxid in taxid_count_dict:
        taxid_queries.append({"taxid": taxid, "abundance":taxid_count_dict[taxid]/total_read_count})
            
    return pd.DataFrame(taxid_queries)

def main():
    parser = argparse.ArgumentParser(description="Universal Taxonomic Classification Verifier.")

    parser.add_argument("-c", "--classification", type=pathlib.Path, required=True, help="Path to the output of Taxonomic Classification Report.")
    parser.add_argument("-s", "--seqscreen", type=pathlib.Path, required=True, help="Path to the Seqscreen output directory.")
    parser.add_argument("--subspecies", action='store_true', required=False, help="Verify taxonomic classification at subspecies rank. Default:[off]")
    parser.set_defaults(subspecies=False)
    
    args = parser.parse_args()
    output_csv = args.classification
    seqscreen_output = args.seqscreen
    call_subspecies = args.subspecies
    
    ncbi_taxa_db = NCBITaxa()
    classification_df = filter_seqscreen_taxonomy(seqscreen_output, ncbi_taxa_db, call_subspecies)
    classification_df.to_csv(output_csv, index=False)
    
if __name__ == "__main__":
    main()
