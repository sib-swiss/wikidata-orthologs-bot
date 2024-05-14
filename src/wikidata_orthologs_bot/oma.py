import argparse
import os
import re
import time
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime

import pandas as pd
import requests
from dotenv import load_dotenv
from tqdm import tqdm
from wikidataintegrator import wdi_core, wdi_helpers, wdi_login

# Load environment variables from .env file
load_dotenv(".env")

# Inspired by https://github.com/SuLab/scheduled-bots/blob/main/scheduled_bots/geneprotein/OrthologBot.py

# Test with:
# TRP53: https://www.wikidata.org/wiki/Q14818098
# Interesting line for TRP53 in orthologs_9606-10090.csv:
# ENSG00000141510,ENSMUSG00000059552,Euarchontoglires,314146

# Insuline: https://www.wikidata.org/wiki/Q21163221
# ENSG00000254647,ENSMUSG00000000215,Euarchontoglires,314146

# Wikidata integrator notebook example: https://public-paws.wmcloud.org/46883698/example%20ema%20bot.ipynb

# This script currently adds a reference to the OMA browser URL for each ortholog pair.
# It does not check if the ortholog pair already exists in Wikidata, so it might create duplicates if ran multiple times


@dataclass
class WdProp:
    instance_of = "P31"
    ensembl_gene_id = "P594"
    entrez_gene_id = "P351"  # aka. NCBI gene ID
    ncbi_taxonomy_id = "P685"
    uniprot_id = "P352"
    encodes = "P688"
    ortholog = "P684"
    found_in_taxon = "P703"
    reference_url = "P854"
    stated_in = "P248"


prop = WdProp()

OMA_WD_ITEM_ID = "Q7104801"
GLOBAL_REF_MODE = "STRICT_KEEP_APPEND"
EDIT_SUMMARY = "Add orthologs from OMA (Orthologous MAtrix) database"
DATA_DIR = "data"
OMA_DIR = os.path.join(DATA_DIR, "oma")


def download_oma_files():
    zip_file_path = os.path.join(DATA_DIR, "OMA_orthologs.zip")
    os.makedirs(OMA_DIR, exist_ok=True)
    if not os.listdir(OMA_DIR):
        print("Downloading OMA files")
        if not os.path.isfile(zip_file_path):
            resp = requests.get("https://www.bgee.org/ftp/current/homologous_genes/OMA_orthologs.zip", timeout=600)
            with open(zip_file_path, "wb") as file:
                file.write(resp.content)
        with zipfile.ZipFile(zip_file_path, "r") as zip_ref:
            zip_ref.extractall(OMA_DIR)
    else:
        print("Files already exist in ./data/oma/, skipping download")


def is_oma_url_valid(oma_url: str) -> bool:
    """Check if the OMA URL for a given UniProt ID responds."""
    try:
        response = requests.head(oma_url, timeout=5)
        return response.status_code == 200
    except requests.RequestException:
        return False


def main():
    start_time = time.time()
    parser = argparse.ArgumentParser(description="Process OMA ortholog pairs")
    parser.add_argument(
        "--write", action="store_true", default=False, help="Enable write mode which requires login and password"
    )
    args = parser.parse_args()
    write = args.write

    if "WDUSER" in os.environ and "WDPASS" in os.environ:
        wd_user = os.environ["WDUSER"]
        wd_login = wdi_login.WDLogin(user=wd_user, pwd=os.environ["WDPASS"])
    else:
        raise ValueError("WDUSER and WDPASS must be specified in `.env` file or as environment variables")

    wdi_core.WDItemEngine.setup_logging(
        log_dir="./data/logs",
        logger_name="WD_logger",
        log_name=f"OmaOrthologBot-{datetime.now().strftime('%Y%m%d_%H:%M')}.log",
    )
    download_oma_files()

    # Map gene IDs to Wikidata IDs
    # id_mapper = {**wdi_helpers.id_mapper(prop.entrez_gene_id), **wdi_helpers.id_mapper(prop.ensembl_gene_id)}
    geneid_to_wdid: dict[str, str] = wdi_helpers.id_mapper(prop.ensembl_gene_id)
    print(f"🔢 {len(geneid_to_wdid)} Ensembl IDs found in Wikidata")
    # Map NCBI taxon ID to wikidata ID
    taxon_to_wdid = wdi_helpers.id_mapper(prop.ncbi_taxonomy_id)

    # NOTE: Mapping Wikidata gene IDs to Wikidata protein IDs. Does not work much because multiple proteins can be encoded by a gene
    # gene_encodes_prot_mapper = {
    #     v: k.replace("http://www.wikidata.org/entity/", "")
    #     for k, v in wdi_helpers.id_mapper(prop.encodes).items()
    # }

    # first_five_items = {k: id_mapper[k] for k in list(id_mapper)[:5]}
    # print(first_five_items)

    valid_orthos = []
    oma_url_found = set()
    wd_items = {}
    errors = defaultdict(set)

    pattern = re.compile(r"orthologs_(\d+)-(\d+)\.csv")
    for filename in tqdm(os.listdir(OMA_DIR), desc="Processing OMA orthologs files", unit="file"):
        if filename.endswith(".csv"):
            # NOTE: to filter for human orthologs add `and "9606" in filename`
            # print(f"📂 Processing {filename}")
            file_path = os.path.join(OMA_DIR, filename)
            try:
                df = pd.read_csv(file_path)
            except Exception as _e:
                errors["error_loading_files"].add(filename)
                continue

            match_taxon = pattern.search(filename)
            if not match_taxon:
                print(f"❌ Taxon not found in filename: {filename}")
                continue

            # Extract taxons from filename
            taxon1, taxon2 = match_taxon.groups()
            if taxon1 not in taxon_to_wdid:
                errors["taxon_not_found"].add(taxon1)
                continue
            if taxon2 not in taxon_to_wdid:
                errors["taxon_not_found"].add(taxon2)
                continue

            # Find missing genes before filtering
            missing_gene1 = set(df["gene1"].astype(str)) - set(geneid_to_wdid.keys())
            missing_gene2 = set(df["gene2"].astype(str)) - set(geneid_to_wdid.keys())
            errors["gene_not_found"].update(missing_gene1)
            errors["gene_not_found"].update(missing_gene2)

            # Remove all rows for which the gene1 or gene2 are not in the id_mapper
            df = df[(df["gene1"].astype(str).isin(geneid_to_wdid)) & (df["gene2"].astype(str).isin(geneid_to_wdid))]
            # If no rows remain after filtering, skip processing
            if df.empty:
                continue

            for _i, row in df.iterrows():
                gene1 = str(row["gene1"])
                gene2 = str(row["gene2"])

                if gene1 not in wd_items:
                    try:
                        prot1 = wdi_core.WDItemEngine(wd_item_id=geneid_to_wdid[gene1]).get_wd_json_representation()[
                            "claims"
                        ][prop.encodes]
                    except Exception as _e:
                        errors["missing_encodes_infos"].add(gene1)
                        errors["pairs_lost_due_to_missing_prot_info"].add(gene1 + gene2)
                        continue
                    wd_items[gene1] = prot1
                else:
                    prot1 = wd_items[gene1]

                if gene2 not in wd_items:
                    try:
                        prot2 = wdi_core.WDItemEngine(wd_item_id=geneid_to_wdid[gene2]).get_wd_json_representation()[
                            "claims"
                        ][prop.encodes]
                    except Exception as _e:
                        errors["missing_encodes_infos"].add(gene2)
                        errors["pairs_lost_due_to_missing_prot_info"].add(gene1 + gene2)
                        continue
                    wd_items[gene2] = prot2
                else:
                    prot2 = wd_items[gene2]

                # TODO: add for multiple proteins? Currently we just add for the 1st when many
                if len(prot1) > 1:
                    errors["more_than_1_encodes"].add(gene1)
                if len(prot2) > 1:
                    errors["more_than_1_encodes"].add(gene2)

                try:
                    prot1_wdid = prot1[0]["mainsnak"]["datavalue"]["value"]["id"]
                    prot1_uniprot = prot1[0]["references"][0]["snaks"][prop.uniprot_id][0]["datavalue"]["value"]
                except Exception as _e:
                    errors["missing_encodes_infos"].add(gene1)
                    continue

                try:
                    prot2_wdid = prot2[0]["mainsnak"]["datavalue"]["value"]["id"]
                    prot2_uniprot = prot2[0]["references"][0]["snaks"][prop.uniprot_id][0]["datavalue"]["value"]
                except Exception as _e:
                    errors["missing_encodes_infos"].add(gene2)
                    continue

                oma_url1 = f"https://omabrowser.org/oma/vps/{prot1_uniprot}/"
                oma_url2 = f"https://omabrowser.org/oma/vps/{prot2_uniprot}/"

                if oma_url1 not in oma_url_found:
                    if not is_oma_url_valid(oma_url1):
                        errors["oma_url_not_found"].add(oma_url1)
                        continue
                    else:
                        oma_url_found.add(oma_url1)

                if oma_url2 not in oma_url_found:
                    if not is_oma_url_valid(oma_url2):
                        errors["oma_url_not_found"].add(oma_url2)
                        continue
                    else:
                        oma_url_found.add(oma_url2)

                # Generate references using OMA wikidata item + reference URL to OMA

                # item1 = wdi_core.WDItemID(value=prot2_wdid, prop_nr=prop.ortholog, references=[(prop.stated_in, oma_wd_item), (prop.reference_url, oma_url1)])
                # Docs: https://github.com/SuLab/WikidataIntegrator/blob/main/wikidataintegrator/wdi_core.py#L300
                # https://github.com/SuLab/WikidataIntegrator/blob/main/wikidataintegrator/wdi_helpers/__init__.py#L59
                item1 = wdi_core.WDItemEngine(
                    wd_item_id=geneid_to_wdid[gene1],
                    data=[
                        wdi_core.WDItemID(
                            geneid_to_wdid[gene2],
                            prop.ortholog,
                            references=[
                                [
                                    wdi_core.WDItemID(OMA_WD_ITEM_ID, prop.stated_in, is_reference=True),
                                    wdi_core.WDExternalID(oma_url1, prop.reference_url, is_reference=True),
                                ]
                            ],
                            qualifiers=[
                                wdi_core.WDItemID(taxon_to_wdid[taxon2], prop.found_in_taxon, is_qualifier=True)
                            ],
                        )
                    ],
                    append_value=[prop.ortholog],
                    global_ref_mode=GLOBAL_REF_MODE,
                    # THIS OVERWRITES, we would need to write a custom ref handler to keep the existing references and only add the OMA browser reference when not present already
                    # global_ref_mode="CUSTOM",
                    # ref_handler=ref_handlers.update_retrieved_if_new,
                    # TODO: the fast_run parameter could avoid writing items that already exist. But we want to update existing items with new references, so not really helpful here
                    # fast_run_use_refs=True,
                    # fast_run=True,
                    # fast_run_base_filter={prop.ensembl_gene_id: "", prop.found_in_taxon: taxon_to_wdid[taxon1]},
                    # core_props=core_props
                )
                wdi_helpers.try_write(
                    item1,
                    gene1,
                    prop.ensembl_gene_id,
                    edit_summary=EDIT_SUMMARY,
                    login=wd_login,
                    write=write,
                )

                # Add the reverse ortholog
                item2 = wdi_core.WDItemEngine(
                    wd_item_id=geneid_to_wdid[gene2],
                    data=[
                        wdi_core.WDItemID(
                            geneid_to_wdid[gene1],
                            prop.ortholog,
                            references=[
                                [
                                    wdi_core.WDItemID(OMA_WD_ITEM_ID, prop.stated_in, is_reference=True),
                                    wdi_core.WDExternalID(oma_url2, prop.reference_url, is_reference=True),
                                ]
                            ],
                            qualifiers=[
                                wdi_core.WDItemID(taxon_to_wdid[taxon1], prop.found_in_taxon, is_qualifier=True)
                            ],
                        )
                    ],
                    append_value=[prop.ortholog],
                    global_ref_mode=GLOBAL_REF_MODE,
                )
                wdi_helpers.try_write(
                    item2,
                    gene2,
                    prop.ensembl_gene_id,
                    edit_summary=EDIT_SUMMARY,
                    login=wd_login,
                    write=write,
                )

                # print(json.dumps(item1.get_wd_json_representation(), indent=2))
                # if "omabrowser" in json.dumps(item1.get_wd_json_representation(), indent=2):
                #     print("omabrowser reference added in item1")
                # if "omabrowser" in json.dumps(item2.get_wd_json_representation(), indent=2):
                #     print("omabrowser reference added in item2")

                valid_orthos.append(
                    {
                        "gene1": gene1,
                        "gene2": gene2,
                        "gene1_wdid": geneid_to_wdid[gene1],
                        "gene2_wdid": geneid_to_wdid[gene2],
                        "prot1_wdid": prot1_wdid,
                        "prot2_wdid": prot2_wdid,
                        "prot1_uniprot": prot1_uniprot,
                        "prot2_uniprot": prot2_uniprot,
                        "taxon1_wdid": taxon_to_wdid[taxon1],
                        "taxon2_wdid": taxon_to_wdid[taxon2],
                        "oma_url1": oma_url1,
                        "oma_url2": oma_url2,
                    }
                )

    for error_type, error_list in errors.items():
        print(f"⚠️ {error_type}: {len(error_list)}")
    print(f"✅ {len(valid_orthos)} orthologs pairs found with proper OMA URL")

    valid_orthos_df = pd.DataFrame(valid_orthos)
    valid_orthos_df.to_csv("valid_ortholog_pairs.csv", index=False)

    print(f"⏱️ Total time taken: {(time.time() - start_time)/60:.2f} minutes")


if __name__ == "__main__":
    main()
