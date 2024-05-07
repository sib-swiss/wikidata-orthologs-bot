import os
import time
from dataclasses import dataclass

import pandas as pd
import requests
from tqdm import tqdm
from wikidataintegrator import wdi_core, wdi_helpers

# Not parallelized, really slow


@dataclass
class WdProp:
    instance_of = "P31"
    ensembl_gene_id = "P594"
    uniprot_id = "P352"
    encodes = "P688"
    ortholog = "P684"
    found_in_taxon = "P703"


def is_oma_url_valid(oma_url: str) -> bool:
    """Check if the OMA URL for a given UniProt ID responds."""
    try:
        response = requests.head(oma_url, timeout=5)
        return response.status_code == 200
    except requests.RequestException:
        return False


if __name__ == "__main__":
    start_time = time.time()
    ensembl_id_mapper: dict = wdi_helpers.id_mapper(WdProp.ensembl_gene_id)
    print(f"🔢 {len(ensembl_id_mapper)} Ensembl IDs found in Wikidata")

    data_folder = "data"
    genes_not_found = set()
    error_loading_files = set()
    more_than_1_encodes = set()
    missing_encodes_infos = set()
    oma_url_not_found = set()
    oma_url_found = set()
    wd_items = {}
    valid_orthos = []

    for filename in tqdm(os.listdir(data_folder), desc="Processing OMA orthologs files", unit="file"):
        if filename.endswith(".csv"):
            # NOTE: to filter for human orthologs add `and "9606" in filename`
            # print(f"📂 Processing {filename}")
            file_path = os.path.join(data_folder, filename)
            try:
                df = pd.read_csv(file_path)
            except Exception as _e:
                error_loading_files.add(filename)
                continue

            for _i, row in df.iterrows():
                gene1 = str(row["gene1"])
                gene2 = str(row["gene2"])
                # TODO: when the ID is only integers, it's NCBIGene ID, not ensembl
                if gene1 not in ensembl_id_mapper:
                    genes_not_found.add(gene1)
                    continue
                if gene2 not in ensembl_id_mapper:
                    genes_not_found.add(gene2)
                    continue

                # print(f"wd id! {ensembl_id_mapper[gene1]}")
                if gene1 not in wd_items:
                    try:
                        prot1 = wdi_core.WDItemEngine(wd_item_id=ensembl_id_mapper[gene1]).get_wd_json_representation()[
                            "claims"
                        ][WdProp.encodes]
                    except Exception as _e:
                        missing_encodes_infos.add(gene1)
                        continue
                    wd_items[gene1] = prot1
                else:
                    prot1 = wd_items[gene1]

                if gene2 not in wd_items:
                    try:
                        prot2 = wdi_core.WDItemEngine(wd_item_id=ensembl_id_mapper[gene2]).get_wd_json_representation()[
                            "claims"
                        ][WdProp.encodes]
                    except Exception as _e:
                        missing_encodes_infos.add(gene2)
                        continue
                    wd_items[gene2] = prot2
                else:
                    prot2 = wd_items[gene2]

                if len(prot1) > 1:
                    more_than_1_encodes.add(gene1)
                    # print(f"⚠️ {gene1} has more than 1 encodes property")
                    # continue
                if len(prot2) > 1:
                    more_than_1_encodes.add(gene2)
                    # continue

                try:
                    prot1_wdid = prot1[0]["mainsnak"]["datavalue"]["value"]["id"]
                    prot1_uniprot = prot1[0]["references"][0]["snaks"][WdProp.uniprot_id][0]["datavalue"]["value"]
                except Exception as _e:
                    missing_encodes_infos.add(gene1)
                    continue

                try:
                    prot2_wdid = prot2[0]["mainsnak"]["datavalue"]["value"]["id"]
                    prot2_uniprot = prot2[0]["references"][0]["snaks"][WdProp.uniprot_id][0]["datavalue"]["value"]
                except Exception as _e:
                    missing_encodes_infos.add(gene2)
                    continue

                oma_url1 = f"https://omabrowser.org/oma/vps/{prot1_uniprot}/"
                oma_url2 = f"https://omabrowser.org/oma/vps/{prot2_uniprot}/"

                if oma_url1 not in oma_url_found:
                    if not is_oma_url_valid(oma_url1):
                        oma_url_not_found.add(oma_url1)
                        continue
                    else:
                        oma_url_found.add(oma_url1)

                if oma_url2 not in oma_url_found:
                    if not is_oma_url_valid(oma_url2):
                        oma_url_not_found.add(oma_url2)
                        continue
                    else:
                        oma_url_found.add(oma_url2)

                print("🔗", gene1, gene2, oma_url1, oma_url2)
                # It takes about 0.5s to 2s for each pair. 1 file has 4500 pairs... So 1 hour
                valid_orthos.append(
                    {
                        "gene1": gene1,
                        "gene2": gene2,
                        "gene1_wdid": ensembl_id_mapper[gene1],
                        "gene2_wdid": ensembl_id_mapper[gene2],
                        "prot1_wdid": prot1_wdid,
                        "prot2_wdid": prot2_wdid,
                        "prot1_uniprot": prot1_uniprot,
                        "prot2_uniprot": prot2_uniprot,
                        "oma_url1": oma_url1,
                        "oma_url2": oma_url2,
                    }
                )

                # TODO: get UniProt ID for these genes `encodes`
                # Generate OMA URL: https://omabrowser.org/oma/vps/A0A3B3BVC0/
                # Some missing match in OMA browser: P62805 (Q22676641)

                # print(prot1_json)
                # print("CLAIMING", prot1_wdid, prot1_uniprot)
                # print(prot1.get_distinct_value_props(property_constraint_pid=WdProp.encodes))
                # encodes_value = prot1.get_single_value(WdProp.encodes)
                # if encodes_value:
                #     print(f"🔍 Property {encodes_prop} exists for item {ensembl_id_mapper[gene1]} with value {encodes_value}")
                # else:
                #     print(f"❌ Property {encodes_prop} does not exist for item {ensembl_id_mapper[gene1]}")

    print(f"⚠️ {len(genes_not_found)} genes not found in Wikidata")
    print(f"⚠️ {len(error_loading_files)} files could not be loaded")
    print(f"⚠️ {len(more_than_1_encodes)} genes that encodes for more than 1 proteins found in Wikidata")
    print(f"⚠️ {len(missing_encodes_infos)} genes with missing encodes infos in Wikidata")
    print(f"⚠️ {len(oma_url_not_found)} OMA URL not found")
    print(f"✅ {len(valid_orthos)} orthologs pairs found with proper OMA URL")

    valid_orthos_df = pd.DataFrame(valid_orthos)
    # csv_output_path = os.path.join(data_folder, "valid_ortholog_pairs.csv")
    valid_orthos_df.to_csv("valid_ortholog_pairs.csv", index=False)

    print(f"⏱️ Total time taken: {(time.time() - start_time)/60:.2f} minutes")

    # Get UniProt ID from Ensemble ID:
    # https://www.ebi.ac.uk/proteins/api/proteins/Ensembl:ENSP00000351276?offset=0&size=100&format=json
    # Many possible sequences for 1 Ensembl ID
    # if target_id.lower().startswith("ensembl:"):
    #     target_id = target_id[len("ensembl:") :]
    #     res = requests.get(
    #         f"https://www.ebi.ac.uk/proteins/api/proteins/Ensembl:{target_id}?offset=0&size=100&format=json",
    #         timeout=TIMEOUT,
    #     ).json()
    #     res[0]["accession"]
    #     return res[0]["sequence"]["sequence"], res[0]["protein"]["recommendedName"]["fullName"]["value"]
