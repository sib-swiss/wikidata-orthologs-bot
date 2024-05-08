import argparse
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
    ncbi_gene_id = "P351"
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


def aggregate_results(results, valid_orthos: list[dict[str, str]], combined_errors: dict[str, set]) -> tuple:
    """Agregate the results and errors of the parallel processing of the ortholog pairs."""
    for result in results:
        if result["data"]:
            valid_orthos.append(result["data"])
        for error_type, error_set in result["errors"].items():
            combined_errors[error_type].update(error_set)
    return valid_orthos, combined_errors


def process_row(row, id_mapper: dict[str, str]):
    """Process a row of the ortholog pairs DataFrame."""
    errors = defaultdict(set)
    gene1 = str(row["gene1"])
    gene2 = str(row["gene2"])

    try:
        prot1 = wdi_core.WDItemEngine(wd_item_id=id_mapper[gene1]).get_wd_json_representation()["claims"][
            WdProp.encodes
        ]
    except Exception as _e:
        errors["missing_encodes_infos"].add(gene1)
        errors["pairs_lost_due_to_missing_prot_info"].add(gene1 + gene2)
        return {"data": None, "errors": errors}

    try:
        prot2 = wdi_core.WDItemEngine(wd_item_id=id_mapper[gene2]).get_wd_json_representation()["claims"][
            WdProp.encodes
        ]
    except Exception as _e:
        errors["missing_encodes_infos"].add(gene2)
        errors["pairs_lost_due_to_missing_prot_info"].add(gene1 + gene2)
        return {"data": None, "errors": errors}

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
        return {"data": None, "errors": errors}

    try:
        prot2_wdid = prot2[0]["mainsnak"]["datavalue"]["value"]["id"]
        prot2_uniprot = prot2[0]["references"][0]["snaks"][WdProp.uniprot_id][0]["datavalue"]["value"]
    except Exception as _e:
        errors["missing_encodes_infos"].add(gene2)
        return {"data": None, "errors": errors}

    oma_url1 = f"https://omabrowser.org/oma/vps/{prot1_uniprot}/"
    oma_url2 = f"https://omabrowser.org/oma/vps/{prot2_uniprot}/"

    if not is_oma_url_valid(oma_url1):
        errors["oma_url_not_found"].add(oma_url1)
        return {"data": None, "errors": errors}

    if not is_oma_url_valid(oma_url2):
        errors["oma_url_not_found"].add(oma_url2)
        return {"data": None, "errors": errors}

    return {
        "data": {
            "gene1": gene1,
            "gene2": gene2,
            "gene1_wdid": id_mapper[gene1],
            "gene2_wdid": id_mapper[gene2],
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
    parser = argparse.ArgumentParser(description="Process OMA ortholog pairs")
    parser.add_argument(
        "nb_workers", type=int, nargs="?", default=None, help="Number of workers to use for parallel processing"
    )
    args = parser.parse_args()

    if args.nb_workers:
        pandarallel.initialize(nb_workers=args.nb_workers)
    else:
        pandarallel.initialize()

    ensembl_mapper: dict = wdi_helpers.id_mapper(WdProp.ensembl_gene_id)
    print(f"🔢 {len(ensembl_mapper)} Ensembl IDs found in Wikidata")
    ncbigene_mapper: dict = wdi_helpers.id_mapper(WdProp.ncbi_gene_id)
    print(f"🔢 {len(ncbigene_mapper)} NCBI gene IDs found in Wikidata")
    id_mapper = {**ncbigene_mapper, **ensembl_mapper}
    print(f"🔢 {len(id_mapper)} Ensembl and NCBI gene IDs found in Wikidata")
    del ensembl_mapper
    del ncbigene_mapper

    data_folder = "data"
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
                combined_errors["error_loading_files"].add(filename)
                continue

            # Find missing genes before filtering
            missing_gene1 = set(df["gene1"].astype(str)) - set(id_mapper.keys())
            missing_gene2 = set(df["gene2"].astype(str)) - set(id_mapper.keys())
            combined_errors["gene_not_found"].update(missing_gene1)
            combined_errors["gene_not_found"].update(missing_gene2)

            # Remove all rows for which the gene1 or gene2 are not in the id_mapper
            df = df[(df["gene1"].astype(str).isin(id_mapper)) & (df["gene2"].astype(str).isin(id_mapper))]
            # If no rows remain after filtering, skip processing
            if df.empty:
                continue

            # NOTE: we could speed up the process by making a set out of the gene1 and gene2 columns
            # and getting the wikidata IDs + prot IDs for these sets before running the row processing

            results = df.parallel_apply(process_row, axis=1, args=(id_mapper,))
            valid_orthos, combined_errors = aggregate_results(results, valid_orthos, combined_errors)

    for error_type, errors in combined_errors.items():
        print(f"⚠️ {error_type}: {len(errors)}")
    print(f"✅ {len(valid_orthos)} orthologs pairs found with proper OMA URL")

    valid_orthos_df = pd.DataFrame(valid_orthos)
    valid_orthos_df.to_csv("valid_ortholog_pairs.csv", index=False)

    print(f"⏱️ Total time taken: {(time.time() - start_time)/60:.2f} minutes")
