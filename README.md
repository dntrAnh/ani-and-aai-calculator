# ANI and AAI Comparisons

This calculator compares genomes/proteomes and produces identity heatmaps for:

- **ANI (Average Nucleotide Identity)** from genome FASTA files
- **AAI (Average Amino Acid Identity)** from protein FASTA files

The project is designed to run **locally** on your machine (not as a hosted production app).

## What This Project Does

You can use this repo in two ways:

- **Command line scripts** for ANI or AAI analysis
- **Flask web app** for upload-and-visualize workflow in a browser

Outputs are PNG heatmaps saved under `outputs/`.

## How ANI Works Here

ANI is handled by `ani_calculation.py`, which wraps `pyANI-plus`.

Workflow:

1. Collect genome FASTA files from paths/directories.
2. Stage files in a temporary input folder.
3. Run `pyani-plus` (default method: `fastani`) into an SQLite database.
4. Read identity values from the database and generate a heatmap via `generate_heatmap.py`.

Supported ANI methods in this wrapper include:

- `fastani` (default)
- `anib`
- `anim`
- `dnadiff`
- `sourmash`

## How AAI Works Here

AAI is implemented in `aai_calculation.py` and is phage-compatible.

Workflow:

1. Read protein FASTA files into proteomes.
2. Use a **k-mer prefilter** (k=4) to find top candidate protein matches.
3. Run local protein alignments (BioPython `PairwiseAligner`) on candidates.
4. Compute **bidirectional best hits (BBH)** between proteomes.
5. Average BBH percent identities to get pairwise AAI.
6. Build all-vs-all AAI matrix and draw a heatmap.

This avoids external alignment tool requirements for AAI and runs fully in Python.

## Local Setup

From the project root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install flask biopython numpy pandas matplotlib seaborn pyani-plus
```

### ANI tool dependency note

For ANI with `fastani`, `pyANI-plus` also relies on ANI backend binaries (for example, `fastANI`).

This code attempts to discover tools on your `PATH` and also checks:

- `~/miniconda3/envs/ani-tools/bin`

If ANI runs fail, ensure those backend binaries are installed and available on `PATH`.

## Run Locally (CLI)

### ANI from genome FASTAs

```bash
source .venv/bin/activate
python ani_calculation.py fasta_files/ -m fastani -d ani_results.db
```

ANI heatmap output:

- `outputs/ani_heatmap.png`

### AAI from protein FASTAs

```bash
source .venv/bin/activate
python aai_calculation.py protein_files/ --output outputs/aai_heatmap.png
```

AAI heatmap output:

- `outputs/aai_heatmap.png`

## Run Locally (Web App)

Start Flask app:

```bash
source .venv/bin/activate
python app.py
```

Then open:

- `http://localhost:8080`

From the web UI, upload at least 2 FASTA files for ANI or AAI.

## Project Structure

- `app.py`: Flask frontend and analysis endpoints
- `ani_calculation.py`: ANI wrapper around `pyANI-plus`
- `aai_calculation.py`: AAI implementation (k-mer + BBH)
- `generate_heatmap.py`: Shared heatmap generation utilities
- `fasta_files/`: Example genome nucleotide FASTAs
- `protein_files/`: Example protein FASTAs
- `outputs/`: Generated heatmaps

## Notes for Public Repo

When you publish this repository:

- Keep it as a local-analysis toolkit for class use
- Mention expected file formats (`.fa`, `.fas`, `.fasta`, `.fna`, `.faa`, `.pep`)
- Optionally add a `requirements.txt` for one-command dependency install

## Sources and References
- Jain, C., Rodriguez-R, L.M., Phillippy, A.M., Konstantinidis, K.T. and Aluru, S. (2018) High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries, Nature Communications, 9(1). doi: 10.1038/s41467-018-07641-9.
- Ruiz-Perez, C. (2024) FastAAI. Available at: https://github.com/cruizperez/FastAAI (Accessed: 14 April 2026).
- Jain, C. (2023) FastANI. Available at: https://github.com/ParBLiSS/FastANI (Accessed: 14 April 2026).