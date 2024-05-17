# 🧬 Wikidata Orthologs Bot 🤖

> This repository uses [`hatch`](https://hatch.pypa.io/latest/) to easily handle scripts and virtual environments. Checkout the `pyproject.toml` file for more details on the scripts available. 
>
> You can also just install dependencies with `pip install .` and run the python script in `src`

## Add orthologs from OMA

Extension of [SuLab OrthologBot.py](https://github.com/SuLab/scheduled-bots/blob/main/scheduled_bots/geneprotein/OrthologBot.py) to include [OMA orthologs](https://omabrowser.org/oma/home/) and references to the OMA browser (e.g. https://omabrowser.org/oma/vps/P04637/).

Define the Wikidata bot username and password in a `.env` file at the root of the repository:

```bash
WDUSER=BOT_USERNAME
WDPASS=BOT_PASSWORD
```

Run mapping script without writing to Wikidata, will generate a CSV file with all orthologs:

```bash
hatch run oma
```

> [!NOTE]
>
> Takes about 34h to run

Run the mapping script with writing to Wikidata enabled:

```bash
hatch run oma --write
```

> [!WARNING]
>
> It currently does not check if the OMA browser reference already exists in Wikidata, so it might create duplicates references if ran multiple times with `--write` enabled

## See also

We use [WikidataIntegrator](https://github.com/SuLab/WikidataIntegrator) to interact with WikiData.

https://github.com/BgeeDB/Wikidata_BgeeDB-bot

