# 🧬 Wikidata Orthologs Bot 🤖

> This repository uses [`hatch`](https://hatch.pypa.io/latest/) to easily handle scripts and virtual environments. Checkout the `pyproject.toml` file for more details on the scripts available.

## Add orthologs from OMA

Extension of [SuLab OrthologBot.py](https://github.com/SuLab/scheduled-bots/blob/main/scheduled_bots/geneprotein/OrthologBot.py) to include OMA orthologs.

1. Download orthologs data:

    ```bash
    mkdir -p data
    cd data
    wget https://www.bgee.org/ftp/current/homologous_genes/OMA_orthologs.zip
    unzip OMA_orthologs.zip
    ```

2. Run mapping script in parallel, will generate a CSV file with all orthologs:

    ```bash
    hatch run oma
    ```

    > You can manually choose the number of workers with `hatch run oma 32`

## See also

We use [WikidataIntegrator](https://github.com/SuLab/WikidataIntegrator) to interact with WikiData.

https://github.com/BgeeDB/Wikidata_BgeeDB-bot
