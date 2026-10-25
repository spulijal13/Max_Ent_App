# MaxEnt Splice Site Scoring Tool

A Streamlit app for calculating MaxEntScan scores for 5′ donor (`5ss`) and 3′
acceptor (`3ss`) splice-site sequences. Supply sequences directly, or supply
human hg38 genomic intervals and let the app extract and score the sequences.

Hosted app: [max-ent-hertel.streamlit.app](https://max-ent-hertel.streamlit.app/).
The hosted deployment may differ from the current local checkout.

## What the app supports

- Individual sequences and batches of sequences.
- Uploaded CSV, TSV, or TXT files, or pasted text.
- Inputs with headers or fixed-order headerless rows.
- 5ss or 3ss sequence scoring; 5ss, 3ss, or Both for coordinates.
- Result tables, CSV and tab-separated TXT downloads, and copy buttons.

MaxEnt scores are model scores, not probabilities. Negative values are valid.
Interpret predictions alongside biological or experimental evidence; this app
uses no classification threshold.

## Local setup

Run commands from the repository root. You need Python 3, pandas, Streamlit,
and Perl. Coordinate mode additionally needs R and the genome packages below.
Use a recent Streamlit release supporting `st.container(horizontal=True)`.

### Python environment

On macOS or Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install streamlit pandas
```

On Windows, activate the environment with `.venv\Scripts\activate` instead.
The repository's `requirements.txt` is currently empty, so installing from that
file alone will not install the app's dependencies. Dependency versions are not
currently pinned.

### Perl

Install Perl if it is not already available, and ensure it is on your `PATH`:

```bash
perl -v
```

Keep the entire `scripts/` directory, including `splicemodels/`, alongside the
Python code. The Perl scorers read these bundled model tables.

### R and hg38 (coordinate mode only)

Install R and confirm that `Rscript` is available in the same terminal used to
launch Streamlit:

```bash
Rscript --version
```

Run the following in an R console to install the required Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(c(
    "Biostrings",
    "GenomicRanges",
    "BSgenome.Hsapiens.UCSC.hg38"
))
```

The genome package is a substantial download. It must be installed in an R
library visible to `Rscript`. Direct sequence scoring does not execute R.

### Launch and stop

```bash
python -m streamlit run app.py
```

Open [localhost:8501](http://localhost:8501). Stop the server with `Ctrl+C`.
If that port is occupied:

```bash
python -m streamlit run app.py --server.port 8502
```

## Using the app

1. Choose **From sequences** or **From raw coordinates**.
2. Select the site type using the horizontal radio buttons.
3. Select **Upload file** or **Paste text**. Only the selected input method is
   displayed and used.
4. Supply the input and click the tab's **Run** button.
5. Inspect the table, download the full results, or copy a sequence-score block.

## Sequence inputs

| Choice | Required length |
| --- | --- |
| 5ss | 9 bases |
| 3ss | 23 bases |

Supply the complete scoring window in the orientation to be scored. The direct
sequence route does not reverse-complement sequences or extract a window.
Leading/trailing whitespace is removed and lowercase bases become uppercase.
Use `A`, `C`, `G`, and `T`; other characters yield a missing score if the length
is otherwise correct. FASTA headers and multiline FASTA records are not supported.

### Paste examples

Select **5ss**:

```text
CAGGTAAGT
AAGGTAAGT
caggtaagt
```

The bundled scorers return `10.86`, `11.00`, and `10.86`, respectively.

Select **3ss**:

```text
TTCCAAACGAACTTTTGTAGGGA
TTCCAAACGAACTTTTGCAGGGA
ttccaaacgaacttttgtaggga
```

The bundled scorers return `2.89`, `3.84`, and `2.89`, respectively. These are
regression examples from the bundled model, not independent biological validation.

### Files and explicit IDs

A CSV with headers:

```csv
row_id,sequence
101,CAGGTAAGT
102,AAGGTAAGT
```

A headerless two-column input:

```text
101 CAGGTAAGT
102 AAGGTAAGT
```

For the most predictable header-based parsing, use the exact column names
`row_id` and `sequence`. Uploaded header-based files use delimiter detection;
headerless inputs and pasted sequence rows accept commas, tabs, or whitespace.
Keep row widths consistent and avoid extra columns in pasted/headerless input.

`row_id` is optional. If omitted, the app creates IDs starting at `0`. If supplied,
use unique integers. The current code does not enforce uniqueness and converts
numeric IDs to integers; fractional IDs can be truncated and duplicate IDs can
multiply rows during merges. Text labels are not supported as IDs.

## Coordinate inputs

Coordinate mode uses **human hg38** through `BSgenome.Hsapiens.UCSC.hg38`.
There is no genome-build selector or coordinate conversion.

- Coordinates are **1-based and inclusive**.
- Supply the complete interval: 9 bases for 5ss or 23 bases for 3ss.
- Keep `start <= end` on both strands.
- Use chromosome identifiers such as `chr1`.
- Use **`1` for the positive strand and `-1` for the negative strand**.

**Current strand limitation:** the R scripts convert `1` to `+` and every other
value to `-`. A literal `+` is therefore handled incorrectly, despite appearing
in the UI placeholder. Use numeric strands as shown here.

The R scripts extract the supplied interval directly, using the requested strand.
They do not locate splice junctions or construct windows around a single position.

### Required headers

| Choice | Required columns besides optional `row_id` |
| --- | --- |
| 5ss | `Chromosome.scaffold.name`, `Strand`, `X5ss_start`, `X5ss_end` |
| 3ss | `Chromosome.scaffold.name`, `Strand`, `X3ss_start`, `X3ss_end` |
| Both | `Chromosome.scaffold.name`, `Strand`, `X5ss_start`, `X5ss_end`, `X3ss_start`, `X3ss_end` |

Use the exact spelling and capitalization. Both mode scores two intervals in
**each input row**; it does not alternate between site types across rows.

### Headerless paste examples

These are arbitrary hg38 intervals for functional testing, not curated splice
sites. Low or negative scores are expected possibilities.

**5ss** — chromosome, strand, 5ss start, 5ss end:

```text
chr1 1 100001 100009
chr1 -1 100001 100009
```

**3ss** — chromosome, strand, 3ss start, 3ss end:

```text
chr1 1 100101 100123
chr1 -1 100101 100123
```

**Both** — chromosome, strand, 5ss start/end, 3ss start/end:

```text
chr1 1 100001 100009 100101 100123
chr1 -1 100001 100009 100101 100123
```

To supply IDs in a headerless file, prepend an integer to every row. Use exactly
4 or 5 columns for single-site inputs, and 6 or 7 columns for Both.

The same Both example as CSV:

```csv
row_id,Chromosome.scaffold.name,Strand,X5ss_start,X5ss_end,X3ss_start,X3ss_end
201,chr1,1,100001,100009,100101,100123
202,chr1,-1,100001,100009,100101,100123
```

## Results, downloads, and copying

- **Direct sequences:** `row_id`, normalized `sequence`, and `score`.
- **Single-site coordinates:** original input columns plus `row_id`, `sequence`,
  and `score`.
- **Both coordinates:** original input columns plus `row_id`, `five_ss_sequence`,
  `five_ss_score`, `three_ss_sequence`, and `three_ss_score`.

The two adjacent buttons immediately below the table download all result columns:

| Button | File contents |
| --- | --- |
| Download results as TXT | Tab-separated text with a header row (`.txt`) |
| Download results as CSV | Comma-separated text with a header row (`.csv`) |

Files are UTF-8 and omit the pandas index. Filenames include input mode and site
type, for example `coordinates_Both_results.csv`. The table also has Streamlit's
built-in CSV export control.

Below the buttons, each labeled text block has a copy icon. These blocks contain
sequence-score columns and headers, separated by tabs; they exclude row IDs and
coordinate metadata. Both mode provides separate 5ss and 3ss blocks and a combined
block.

Perl formats scores to two decimal places, but pandas or the UI may display fewer
trailing zeros. Perl's `NA` becomes a missing numeric value and may appear blank
in downloaded files.

## How the pipeline works

```text
Sequence input -> parsing/normalization -> length check -> Perl scorer -> merge -> results
Coordinate input -> parsing/validation -> R hg38 extraction -> Perl scorer -> merge -> results
```

Python assigns or normalizes row IDs and uses them to join output back to input.
For Both, it runs extraction and scoring once per site type and merges both
outputs. Each run uses a temporary directory that is cleaned up when processing
finishes, including when an exception occurs.

R and Perl run as subprocesses with `scripts/` as their working directory so
relative model-table paths resolve correctly. Nonzero subprocess exits become
errors containing the command, standard output, and standard error. Parsing and
processing errors are displayed in the app.

## Project structure

```text
app.py                  Streamlit entry point, tabs, selectors, and input widgets
maxent/
    __init__.py         Python package marker
    inputs.py           Parsing, sequence normalization, row IDs, coordinate validation
    scoring.py          Script paths, subprocess execution, temporary files, result joins
    results.py          Result table, download buttons, and copyable output blocks
scripts/
    5prime.r            Extract 5ss intervals from hg38
    3prime.r            Extract 3ss intervals from hg38
    score5.pl           Score 9-base sequences
    score3.pl           Score 23-base sequences
    splicemodels/       Bundled model tables and supporting files
requirements.txt        Currently empty; install dependencies explicitly
README.md               This guide
test_data/              Optional local manual fixtures (currently Git-ignored)
```

UI changes belong in `app.py` or `maxent/results.py`; parsing changes belong in
`maxent/inputs.py`; workflow changes belong in `maxent/scoring.py`. The Python
module split preserved the R and Perl scripts.

## Manual testing

If your checkout includes `test_data/`, start with these five files:

| File | Tab and choice |
| --- | --- |
| `sequences_5ss.csv` | From sequences / 5ss |
| `sequences_3ss.csv` | From sequences / 3ss |
| `coordinates_5ss.csv` | From raw coordinates / 5ss |
| `coordinates_3ss.csv` | From raw coordinates / 3ss |
| `coordinates_Both.csv` | From raw coordinates / Both |

Other fixture files cover additional input paths:

- `.tsv` and `_headerless.txt`: alternate delimiters, pasted rows, and generated IDs.
- `expected_sequences_*`: expected outputs from the bundled Perl scorers.
- `invalid_*`: wrong lengths, invalid IDs, missing columns, and non-numeric coordinates.
- `na_*`: correctly sized sequences containing an ambiguous base.

See `test_data/README.md` for the detailed checklist when that directory is present.
It is ignored by Git, so a fresh clone may not include it; the inline examples
above can also be used for testing. These files are input data, not executable
scripts, and the app does not need them to run.

For each successful case, check row counts, IDs, sequence lengths, copy output,
and downloads. Opposite-strand rows for the same interval should yield reverse
complements. Both-mode results should agree with separate single-site runs.
Paste extracted sequences into the sequence tab to compare the two scoring paths.

## Troubleshooting and current limitations

| Symptom | What to check |
| --- | --- |
| `No module named streamlit` or `pandas` | Activate the virtual environment and install the Python dependencies. |
| Unexpected `horizontal` argument | Upgrade Streamlit in the environment used to launch the app. |
| `perl` or `Rscript` not found | Install the runtime and ensure it is on that process's `PATH`. |
| R package missing | Install the Bioconductor packages in the library used by `Rscript`. |
| Missing scripts or model files | Restore the complete `scripts/` directory. Startup checks paths, not every model table. |
| Expected sequence length error | Use 9 bases for 5ss or 23 for 3ss; do not include IDs in the sequence itself. |
| Missing score | Check for non-ACGT characters or a wrong extracted interval length. |
| Missing coordinate columns | Match the exact header names and selected site type. |
| Coordinate extraction fails | Check hg38 chromosome names, integer bounds, interval order, and installed genome package. |
| Unexpected row count | Check for duplicate or fractional row IDs. |
| Results disappear after interaction | Results currently render only on the Run button's execution; run again. |

Additional limitations:

- Coordinate validation checks required columns and numeric values, but not
  integer coordinates, bounds, interval length, or valid strands before R runs.
- Header detection uses keyword matching and delimiter inference. Simple,
  consistent formats with documented headers are the safest inputs.
- Scoring is synchronous and subprocesses have no timeout or cancellation control.
- Results are not persisted in Streamlit session state. Changes to widgets and
  interactions that trigger reruns can remove the displayed output.
- This repository does not currently provide a pinned deployment environment
  or an automated test suite. The manual fixtures are not a substitute for
  independent validation of the model or genome extraction.
