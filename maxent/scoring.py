"""Run the existing R/Perl pipeline and assemble scored results."""

import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from .inputs import ensure_row_id, normalize_sequence_column, validate_coordinate_columns

BASE_DIR = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = BASE_DIR / "scripts"

R_3PRIME = SCRIPTS_DIR / "3prime.r"
R_5PRIME = SCRIPTS_DIR / "5prime.r"
PERL_3PRIME = SCRIPTS_DIR / "score3.pl"
PERL_5PRIME = SCRIPTS_DIR / "score5.pl"

def run_command(cmd, cwd=None):
    """Run a subprocess and raise RuntimeError if it exits unsuccessfully.

    Args:
        cmd (list[str | Path]): Executable and command-line arguments.
        cwd (str | Path | None): Working directory, or the current directory.

    Returns:
        subprocess.CompletedProcess[str]: Exit status and captured text output.
    """
    result = subprocess.run(
        cmd,
        cwd=cwd,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed:\n{' '.join(map(str, cmd))}\n\n"
            f"STDOUT:\n{result.stdout}\n\n"
            f"STDERR:\n{result.stderr}"
        )
    return result

def validate_required_files():
    """Check that the R/Perl scripts and model directory exist.

    Args:
        None.

    Returns:
        None: Raises FileNotFoundError when a required path is missing.
    """
    needed = [
        R_3PRIME,
        R_5PRIME,
        PERL_3PRIME,
        PERL_5PRIME,
        SCRIPTS_DIR / "splicemodels",
    ]
    missing = [str(p) for p in needed if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required files or folders:\n" + "\n".join(missing))

def read_sequence_file(path: Path, seq_col_name: str) -> pd.DataFrame:
    """Read the headerless tab-separated output of an extraction script.

    Args:
        path (Path): Sequence output file.
        seq_col_name (str): Name to assign to the sequence column.

    Returns:
        pd.DataFrame: Integer row_id values and the named sequence column.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["row_id", seq_col_name],
    )
    df["row_id"] = pd.to_numeric(df["row_id"], errors="raise").astype(int)
    return df

def read_score_file(path: Path, seq_col_name: str, score_col_name: str) -> pd.DataFrame:
    """Read Perl output and convert scores to numbers, coercing NA to missing.

    Args:
        path (Path): Headerless tab-separated score file.
        seq_col_name (str): Name to assign to the sequence column.
        score_col_name (str): Name to assign to the score column.

    Returns:
        pd.DataFrame: Integer row_id, sequence, and numeric score columns.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["row_id", seq_col_name, score_col_name],
    )
    df["row_id"] = pd.to_numeric(df["row_id"], errors="raise").astype(int)
    df[score_col_name] = pd.to_numeric(df[score_col_name], errors="coerce")
    return df

def run_r_script(r_script: Path, input_csv: Path, output_txt: Path):
    """Run an R extraction script from the scripts directory.

    Args:
        r_script (Path): R script to execute.
        input_csv (Path): Coordinate CSV to read.
        output_txt (Path): Destination for extracted row IDs and sequences.

    Returns:
        None: The subprocess writes the output file.
    """
    run_command(
        ["Rscript", str(r_script), str(input_csv), str(output_txt)],
        cwd=SCRIPTS_DIR,
    )

def run_perl_script(perl_script: Path, input_txt: Path, output_txt: Path):
    """Run a Perl scoring script from the scripts directory.

    Args:
        perl_script (Path): Perl script to execute.
        input_txt (Path): Tab-separated row IDs and sequences.
        output_txt (Path): Destination for row IDs, sequences, and scores.

    Returns:
        None: The subprocess writes the output file.
    """
    run_command(
        ["perl", str(perl_script), str(input_txt), str(output_txt)],
        cwd=SCRIPTS_DIR,
    )

def write_sequence_input(df: pd.DataFrame, seq_col: str, output_txt: Path):
    """Write row IDs and sequences as headerless tab-separated text.

    Args:
        df (pd.DataFrame): Rows containing row_id and the requested sequence column.
        seq_col (str): Sequence column to write.
        output_txt (Path): Destination file.

    Returns:
        None: Writes the selected columns to disk.
    """
    out = df[["row_id", seq_col]].copy()
    out.to_csv(output_txt, sep="\t", index=False, header=False)

def score_sequences_direct(df: pd.DataFrame, site_type: str, tmpdir: Path) -> pd.DataFrame:
    """Normalize and validate sequences, run Perl, and join scores by row ID.

    Args:
        df (pd.DataFrame): Sequence input with optional row_id values.
        site_type (str): "5ss" or "3ss".
        tmpdir (Path): Existing directory for intermediate input and output files.

    Returns:
        pd.DataFrame: row_id, normalized sequence, and score columns.
    """
    df = ensure_row_id(df)
    df = normalize_sequence_column(df)

    expected_len = 9 if site_type == "5ss" else 23
    bad_len = df[df["sequence"].str.len() != expected_len]
    if not bad_len.empty:
        raise ValueError(
            f"Expected sequence length {expected_len} for {site_type}. "
            f"Problem row_ids: {bad_len['row_id'].head(10).tolist()}"
        )

    input_txt = tmpdir / f"{site_type}_input_sequences.txt"
    score_txt = tmpdir / f"{site_type}_scores.txt"

    write_sequence_input(df, "sequence", input_txt)

    perl_script = PERL_5PRIME if site_type == "5ss" else PERL_3PRIME
    run_perl_script(perl_script, input_txt, score_txt)

    score_df = read_score_file(score_txt, "sequence_scored", "score")
    result = df.merge(score_df[["row_id", "score"]], on="row_id", how="left")

    return result[["row_id", "sequence", "score"]]

def extract_and_score_from_coordinates(df: pd.DataFrame, site_type: str, tmpdir: Path) -> pd.DataFrame:
    """Extract one site type with R, score with Perl, and join the outputs.

    Args:
        df (pd.DataFrame): Validated coordinates with optional row_id values.
        site_type (str): "5ss" or "3ss".
        tmpdir (Path): Existing directory for intermediate input and output files.

    Returns:
        pd.DataFrame: Input columns plus row_id, sequence, and score.
    """
    df = ensure_row_id(df)

    input_csv = tmpdir / "coordinate_input.csv"
    seq_txt = tmpdir / f"{site_type}_sequences.txt"
    score_txt = tmpdir / f"{site_type}_scores.txt"

    df.to_csv(input_csv, index=False)

    if site_type == "5ss":
        run_r_script(R_5PRIME, input_csv, seq_txt)
        run_perl_script(PERL_5PRIME, seq_txt, score_txt)
    else:
        run_r_script(R_3PRIME, input_csv, seq_txt)
        run_perl_script(PERL_3PRIME, seq_txt, score_txt)

    seq_df = read_sequence_file(seq_txt, "sequence")
    score_df = read_score_file(score_txt, "sequence_scored", "score")

    result = df.merge(seq_df, on="row_id", how="left")
    result = result.merge(score_df[["row_id", "score"]], on="row_id", how="left")

    return result

def run_coordinate_mode(df: pd.DataFrame, site_choice: str) -> pd.DataFrame:
    """Validate coordinates and run one or both site workflows in temporary files.

    Args:
        df (pd.DataFrame): Coordinate input with optional row_id values.
        site_choice (str): "5ss", "3ss", or "Both".

    Returns:
        pd.DataFrame: Input rows with sequences and scores; Both uses site-specific columns.
    """
    df = ensure_row_id(df)
    df = validate_coordinate_columns(df, site_choice)

    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)

        if site_choice == "5ss":
            return extract_and_score_from_coordinates(df, "5ss", tmpdir)

        if site_choice == "3ss":
            return extract_and_score_from_coordinates(df, "3ss", tmpdir)

        result5 = extract_and_score_from_coordinates(df, "5ss", tmpdir)[["row_id", "sequence", "score"]]
        result3 = extract_and_score_from_coordinates(df, "3ss", tmpdir)[["row_id", "sequence", "score"]]

        result5 = result5.rename(
            columns={"sequence": "five_ss_sequence", "score": "five_ss_score"}
        )
        result3 = result3.rename(
            columns={"sequence": "three_ss_sequence", "score": "three_ss_score"}
        )

        merged = df.merge(result5, on="row_id", how="left")
        merged = merged.merge(result3, on="row_id", how="left")
        return merged


def run_sequence_mode(df: pd.DataFrame, site_type: str) -> pd.DataFrame:
    """Score sequences and clean up the temporary working directory afterward.

    Args:
        df (pd.DataFrame): Sequence input with optional row_id values.
        site_type (str): "5ss" or "3ss".

    Returns:
        pd.DataFrame: row_id, normalized sequence, and score columns.
    """
    with tempfile.TemporaryDirectory() as td:
        return score_sequences_direct(df, site_type, Path(td))
