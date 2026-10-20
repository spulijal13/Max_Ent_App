import io
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

BASE_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = BASE_DIR / "scripts"

R_3PRIME = SCRIPTS_DIR / "3prime.r"
R_5PRIME = SCRIPTS_DIR / "5prime.r"
PERL_3PRIME = SCRIPTS_DIR / "score3.pl"
PERL_5PRIME = SCRIPTS_DIR / "score5.pl"

st.set_page_config(page_title="MaxEnt Splice Site Tool", layout="wide")
st.title("MaxEnt Splice Site Tool")


def run_command(cmd, cwd=None):
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


def ensure_row_id(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "row_id" not in df.columns:
        df = df.reset_index(drop=True)
        df["row_id"] = df.index
    else:
        df["row_id"] = pd.to_numeric(df["row_id"], errors="coerce")
        if df["row_id"].isna().any():
            raise ValueError("row_id column exists but contains non-numeric values.")
        df["row_id"] = df["row_id"].astype(int)
    return df


def looks_like_header(first_line: str, expected_names=None):
    if not first_line:
        return False
    line = first_line.strip().lower()
    if expected_names:
        return any(name.lower() in line for name in expected_names)
    return False


def read_sequence_file(path: Path, seq_col_name: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["row_id", seq_col_name],
    )
    df["row_id"] = pd.to_numeric(df["row_id"], errors="raise").astype(int)
    return df


def read_score_file(path: Path, seq_col_name: str, score_col_name: str) -> pd.DataFrame:
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
    run_command(
        ["Rscript", str(r_script), str(input_csv), str(output_txt)],
        cwd=SCRIPTS_DIR,
    )


def run_perl_script(perl_script: Path, input_txt: Path, output_txt: Path):
    run_command(
        ["perl", str(perl_script), str(input_txt), str(output_txt)],
        cwd=SCRIPTS_DIR,
    )


def write_sequence_input(df: pd.DataFrame, seq_col: str, output_txt: Path):
    out = df[["row_id", seq_col]].copy()
    out.to_csv(output_txt, sep="\t", index=False, header=False)


def make_copy_text(df: pd.DataFrame, columns):
    cols = [c for c in columns if c in df.columns]
    if not cols:
        return ""
    return df[cols].to_csv(sep="\t", index=False)


def show_copy_boxes(result_df: pd.DataFrame, prefix: str):
    st.markdown("### Copy-friendly output")

    if {"sequence", "score"}.issubset(result_df.columns):
        st.text_area(
            "Sequence and score",
            value=make_copy_text(result_df, ["sequence", "score"]),
            height=220,
            key=f"{prefix}_sequence_score",
        )

    if {"five_ss_sequence", "five_ss_score"}.issubset(result_df.columns):
        st.text_area(
            "5ss sequence and score",
            value=make_copy_text(result_df, ["five_ss_sequence", "five_ss_score"]),
            height=220,
            key=f"{prefix}_5ss_sequence_score",
        )

    if {"three_ss_sequence", "three_ss_score"}.issubset(result_df.columns):
        st.text_area(
            "3ss sequence and score",
            value=make_copy_text(result_df, ["three_ss_sequence", "three_ss_score"]),
            height=220,
            key=f"{prefix}_3ss_sequence_score",
        )

    if {
        "five_ss_sequence",
        "five_ss_score",
        "three_ss_sequence",
        "three_ss_score",
    }.issubset(result_df.columns):
        st.text_area(
            "Both 5ss and 3ss sequence-score pairs",
            value=make_copy_text(
                result_df,
                ["five_ss_sequence", "five_ss_score", "three_ss_sequence", "three_ss_score"],
            ),
            height=260,
            key=f"{prefix}_both_sequence_score",
        )


def normalize_sequence_column(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    if "sequence" in df.columns:
        pass
    else:
        non_row_cols = [c for c in df.columns if c != "row_id"]
        if len(non_row_cols) == 1:
            df = df.rename(columns={non_row_cols[0]: "sequence"})
        else:
            raise ValueError(
                "Sequence input must have either a 'sequence' column or exactly one non-row_id column."
            )

    df["sequence"] = df["sequence"].astype(str).str.strip().str.upper()
    return df


def parse_pasted_sequences(text: str) -> pd.DataFrame:
    lines = [line.strip() for line in text.strip().splitlines() if line.strip()]
    if not lines:
        raise ValueError("No pasted sequence content found.")

    first = lines[0]
    has_header = looks_like_header(first, ["row_id", "sequence"])

    rows = lines[1:] if has_header else lines
    parsed = []

    for line in rows:
        parts = line.replace(",", " ").split()
        if len(parts) == 1:
            parsed.append([parts[0]])
        elif len(parts) >= 2:
            parsed.append([parts[0], parts[1]])

    if not parsed:
        raise ValueError("Could not parse pasted sequences.")

    if all(len(r) == 1 for r in parsed):
        return pd.DataFrame(parsed, columns=["sequence"])
    return pd.DataFrame(parsed, columns=["row_id", "sequence"])


def parse_uploaded_sequences(file_bytes: bytes) -> pd.DataFrame:
    text = file_bytes.decode("utf-8-sig", errors="replace").strip()
    if not text:
        raise ValueError("Uploaded sequence file is empty.")

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    first = lines[0]
    has_header = looks_like_header(first, ["row_id", "sequence"])

    if has_header:
        try:
            return pd.read_csv(io.StringIO(text), sep=None, engine="python")
        except Exception:
            for sep in [",", "\t"]:
                try:
                    return pd.read_csv(io.StringIO(text), sep=sep)
                except Exception:
                    pass
            raise ValueError("Could not parse uploaded sequence file with header.")

    parsed = []
    for line in lines:
        parts = line.replace(",", " ").split()
        if len(parts) == 1:
            parsed.append([parts[0]])
        elif len(parts) >= 2:
            parsed.append([parts[0], parts[1]])

    if not parsed:
        raise ValueError("Could not parse uploaded sequence file.")

    if all(len(r) == 1 for r in parsed):
        return pd.DataFrame(parsed, columns=["sequence"])
    return pd.DataFrame(parsed, columns=["row_id", "sequence"])


def parse_coordinate_input(uploaded_file, pasted_text, site_choice):
    if uploaded_file is not None:
        text = uploaded_file.getvalue().decode("utf-8-sig", errors="replace").strip()
    elif pasted_text and pasted_text.strip():
        text = pasted_text.strip()
    else:
        raise ValueError("Provide either an uploaded file or pasted text.")

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("No coordinate content found.")

    possible_headers = [
        "row_id",
        "chromosome.scaffold.name",
        "strand",
        "x5ss_start",
        "x5ss_end",
        "x3ss_start",
        "x3ss_end",
    ]
    has_header = looks_like_header(lines[0], possible_headers)

    if has_header:
        try:
            df = pd.read_csv(io.StringIO(text), sep=None, engine="python")
        except Exception:
            df = None
            for sep in [",", "\t"]:
                try:
                    df = pd.read_csv(io.StringIO(text), sep=sep)
                    break
                except Exception:
                    pass
            if df is None:
                raise ValueError("Could not parse coordinate input with headers.")
        return df

    parsed = []
    for line in lines:
        parts = line.replace(",", " ").split()
        parsed.append(parts)

    if site_choice == "5ss":
        if all(len(r) == 4 for r in parsed):
            df = pd.DataFrame(
                parsed,
                columns=["Chromosome.scaffold.name", "Strand", "X5ss_start", "X5ss_end"],
            )
        elif all(len(r) >= 5 for r in parsed):
            df = pd.DataFrame(
                parsed,
                columns=["row_id", "Chromosome.scaffold.name", "Strand", "X5ss_start", "X5ss_end"],
            )
        else:
            raise ValueError(
                "Headerless 5ss coordinate input must have either 4 columns "
                "(chr strand X5ss_start X5ss_end) or 5 columns including row_id."
            )

    elif site_choice == "3ss":
        if all(len(r) == 4 for r in parsed):
            df = pd.DataFrame(
                parsed,
                columns=["Chromosome.scaffold.name", "Strand", "X3ss_start", "X3ss_end"],
            )
        elif all(len(r) >= 5 for r in parsed):
            df = pd.DataFrame(
                parsed,
                columns=["row_id", "Chromosome.scaffold.name", "Strand", "X3ss_start", "X3ss_end"],
            )
        else:
            raise ValueError(
                "Headerless 3ss coordinate input must have either 4 columns "
                "(chr strand X3ss_start X3ss_end) or 5 columns including row_id."
            )

    else:
        if all(len(r) == 6 for r in parsed):
            df = pd.DataFrame(
                parsed,
                columns=[
                    "Chromosome.scaffold.name",
                    "Strand",
                    "X5ss_start",
                    "X5ss_end",
                    "X3ss_start",
                    "X3ss_end",
                ],
            )
        elif all(len(r) >= 7 for r in parsed):
            df = pd.DataFrame(
                parsed,
                columns=[
                    "row_id",
                    "Chromosome.scaffold.name",
                    "Strand",
                    "X5ss_start",
                    "X5ss_end",
                    "X3ss_start",
                    "X3ss_end",
                ],
            )
        else:
            raise ValueError(
                "Headerless Both coordinate input must have either 6 columns "
                "(chr strand X5ss_start X5ss_end X3ss_start X3ss_end) "
                "or 7 columns including row_id."
            )

    return df


def validate_coordinate_columns(df: pd.DataFrame, site_type: str):
    base_required = ["Chromosome.scaffold.name", "Strand"]

    if site_type == "5ss":
        site_required = ["X5ss_start", "X5ss_end"]
    elif site_type == "3ss":
        site_required = ["X3ss_start", "X3ss_end"]
    elif site_type == "Both":
        site_required = ["X5ss_start", "X5ss_end", "X3ss_start", "X3ss_end"]
    else:
        raise ValueError("Invalid site type.")

    missing = [c for c in base_required + site_required if c not in df.columns]
    if missing:
        raise ValueError("Coordinate input is missing required columns:\n" + "\n".join(missing))

    for col in site_required:
        df[col] = pd.to_numeric(df[col], errors="coerce")
        if df[col].isna().any():
            raise ValueError(f"Column {col} contains missing or non-numeric values.")

    return df


def score_sequences_direct(df: pd.DataFrame, site_type: str, tmpdir: Path) -> pd.DataFrame:
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


validate_required_files()

tab_seq, tab_coord = st.tabs(["From sequences", "From raw coordinates"])

with tab_seq:
    st.subheader("Run from sequences")
    st.write("Paste one sequence per line, or two columns as row_id and sequence. Headers are optional.")

    seq_site_choice = st.selectbox("Choose sequence type", ["5ss", "3ss"], key="seq_site_choice")
    seq_upload = st.file_uploader("Upload sequence file", type=["csv", "txt", "tsv"], key="seq_upload")
    seq_paste = st.text_area(
        "Or paste sequences here",
        height=240,
        key="seq_paste",
        placeholder="Examples:\nCAGGTAAGT\nAAGGTAAGT\n\nor\n1 CAGGTAAGT\n2 AAGGTAAGT",
    )

    if st.button("Run from sequences", key="run_seq"):
        try:
            if seq_upload is not None:
                seq_df = parse_uploaded_sequences(seq_upload.getvalue())
            elif seq_paste and seq_paste.strip():
                seq_df = parse_pasted_sequences(seq_paste)
            else:
                raise ValueError("Provide either an uploaded file or pasted sequences.")

            with tempfile.TemporaryDirectory() as td:
                tmpdir = Path(td)
                result_df = score_sequences_direct(seq_df, seq_site_choice, tmpdir)

            st.success("Done.")
            st.dataframe(result_df, use_container_width=True)
            show_copy_boxes(result_df, "seq")

            st.download_button(
                "Download results as TXT",
                data=result_df.to_csv(sep="\t", index=False).encode("utf-8"),
                file_name=f"sequences_{seq_site_choice}_results.txt",
                mime="text/plain",
            )
            st.download_button(
                "Download results as CSV",
                data=result_df.to_csv(index=False).encode("utf-8"),
                file_name=f"sequences_{seq_site_choice}_results.csv",
                mime="text/csv",
            )

        except Exception as e:
            st.error(str(e))

with tab_coord:
    st.subheader("Run from raw coordinates")
    st.write("Headers are optional. If there are no headers, use the fixed column order shown below.")

    coord_site_choice = st.selectbox(
        "Choose what to generate and score",
        ["5ss", "3ss", "Both"],
        key="coord_site_choice",
    )

    if coord_site_choice == "5ss":
        st.caption("Headerless format: chr strand X5ss_start X5ss_end  or  row_id chr strand X5ss_start X5ss_end")
    elif coord_site_choice == "3ss":
        st.caption("Headerless format: chr strand X3ss_start X3ss_end  or  row_id chr strand X3ss_start X3ss_end")
    else:
        st.caption(
            "Headerless format: chr strand X5ss_start X5ss_end X3ss_start X3ss_end  "
            "or  row_id chr strand X5ss_start X5ss_end X3ss_start X3ss_end"
        )

    coord_upload = st.file_uploader("Upload coordinate file", type=["csv", "txt", "tsv"], key="coord_upload")
    coord_paste = st.text_area(
        "Or paste coordinate rows here",
        height=240,
        key="coord_paste",
        placeholder="Examples:\nchr1 + 100 108\nchr2 - 500 508",
    )

    if st.button("Run from coordinates", key="run_coord"):
        try:
            coord_df = parse_coordinate_input(coord_upload, coord_paste, coord_site_choice)
            result_df = run_coordinate_mode(coord_df, coord_site_choice)

            st.success("Done.")
            st.dataframe(result_df, use_container_width=True)
            show_copy_boxes(result_df, "coord")

            st.download_button(
                "Download results as TXT",
                data=result_df.to_csv(sep="\t", index=False).encode("utf-8"),
                file_name=f"coordinates_{coord_site_choice}_results.txt",
                mime="text/plain",
            )
            st.download_button(
                "Download results as CSV",
                data=result_df.to_csv(index=False).encode("utf-8"),
                file_name=f"coordinates_{coord_site_choice}_results.csv",
                mime="text/csv",
            )

        except Exception as e:
            st.error(str(e))