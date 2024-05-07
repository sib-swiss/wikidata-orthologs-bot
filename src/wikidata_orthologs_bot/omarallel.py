import os
import time
from collections import defaultdict
from dataclasses import dataclass

import pandas as pd
import requests
from pandarallel import pandarallel
from tqdm import tqdm
from wikidataintegrator import wdi_core, wdi_helpers


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


def aggregate_results(results, valid_orthos, combined_errors) -> tuple:
    for result in results:
        if result["data"]:
            valid_orthos.append(result["data"])
        for error_type, error_set in result["errors"].items():
            combined_errors[error_type].update(error_set)
    return valid_orthos, combined_errors


def process_row(row, ensembl_mapper):
    data = None
    errors = defaultdict(set)
    gene1 = str(row["gene1"])
    gene2 = str(row["gene2"])

    # NOTE: when the ID is only integers, it's NCBIGene ID, not ensembl
    if gene1 not in ensembl_mapper:
        errors["genes_not_found"].add(gene1)
        return {"data": data, "errors": errors}
    if gene2 not in ensembl_mapper:
        errors["genes_not_found"].add(gene2)
        return {"data": data, "errors": errors}

    try:
        prot1 = wdi_core.WDItemEngine(wd_item_id=ensembl_mapper[gene1]).get_wd_json_representation()["claims"][
            WdProp.encodes
        ]
    except Exception as _e:
        errors["missing_encodes_infos"].add(gene1)
        return {"data": data, "errors": errors}

    try:
        prot2 = wdi_core.WDItemEngine(wd_item_id=ensembl_mapper[gene2]).get_wd_json_representation()["claims"][
            WdProp.encodes
        ]
    except Exception as _e:
        errors["missing_encodes_infos"].add(gene2)
        return {"data": data, "errors": errors}

    if len(prot1) > 1:
        errors["more_than_1_encodes"].add(gene1)
        # return {"data": data, "errors": errors}
    if len(prot2) > 1:
        errors["more_than_1_encodes"].add(gene2)
        # return {"data": data, "errors": errors}

    try:
        prot1_wdid = prot1[0]["mainsnak"]["datavalue"]["value"]["id"]
        prot1_uniprot = prot1[0]["references"][0]["snaks"][WdProp.uniprot_id][0]["datavalue"]["value"]
    except Exception as _e:
        errors["missing_encodes_infos"].add(gene1)
        return {"data": data, "errors": errors}

    try:
        prot2_wdid = prot2[0]["mainsnak"]["datavalue"]["value"]["id"]
        prot2_uniprot = prot2[0]["references"][0]["snaks"][WdProp.uniprot_id][0]["datavalue"]["value"]
    except Exception as _e:
        errors["missing_encodes_infos"].add(gene2)
        return {"data": data, "errors": errors}

    oma_url1 = f"https://omabrowser.org/oma/vps/{prot1_uniprot}/"
    oma_url2 = f"https://omabrowser.org/oma/vps/{prot2_uniprot}/"

    if not is_oma_url_valid(oma_url1):
        errors["oma_url_not_found"].add(oma_url1)
        return {"data": data, "errors": errors}

    if not is_oma_url_valid(oma_url2):
        errors["oma_url_not_found"].add(oma_url2)
        return {"data": data, "errors": errors}

    return {
        "data": {
            "gene1": gene1,
            "gene2": gene2,
            "gene1_wdid": ensembl_mapper[gene1],
            "gene2_wdid": ensembl_mapper[gene2],
            "prot1_wdid": prot1_wdid,
            "prot2_wdid": prot2_wdid,
            "prot1_uniprot": prot1_uniprot,
            "prot2_uniprot": prot2_uniprot,
            "oma_url1": oma_url1,
            "oma_url2": oma_url2,
        },
        "errors": {},
    }


if __name__ == "__main__":
    start_time = time.time()
    pandarallel.initialize(progress_bar=True)
    # pandarallel.initialize(progress_bar=True, nb_workers=32)

    ensembl_mapper: dict = wdi_helpers.id_mapper(WdProp.ensembl_gene_id)
    print(f"🔢 {len(ensembl_mapper)} Ensembl IDs found in Wikidata")

    data_folder = "data"
    error_loading_files = set()
    valid_orthos = []
    combined_errors = defaultdict(set)

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

            results = df.parallel_apply(process_row, axis=1, args=(ensembl_mapper,))
            valid_orthos, combined_errors = aggregate_results(results, valid_orthos, combined_errors)

    print(f"⚠️ {len(error_loading_files)} files could not be loaded")
    for error_type, errors in combined_errors.items():
        print(f"{error_type}: {len(errors)} errors")
    print(f"✅ {len(valid_orthos)} orthologs pairs found with proper OMA URL")

    valid_orthos_df = pd.DataFrame(valid_orthos)
    valid_orthos_df.to_csv("valid_ortholog_pairs.csv", index=False)

    print(f"⏱️ Total time taken: {(time.time() - start_time)/60:.2f} minutes")

    # NOTE: Get UniProt ID from Ensemble ID, e.g. https://www.ebi.ac.uk/proteins/api/proteins/Ensembl:ENSP00000351276?offset=0&size=100&format=json
    # Many possible sequences for 1 Ensembl ID
    # if target_id.lower().startswith("ensembl:"):
    #     target_id = target_id[len("ensembl:") :]
    #     res = requests.get(
    #         f"https://www.ebi.ac.uk/proteins/api/proteins/Ensembl:{target_id}?offset=0&size=100&format=json",
    #         timeout=TIMEOUT,
    #     ).json()
    #     res[0]["accession"]
    #     return res[0]["sequence"]["sequence"], res[0]["protein"]["recommendedName"]["fullName"]["value"]
