"""Parse and validate sequence and coordinate inputs."""

import io

import pandas as pd


def ensure_row_id(df: pd.DataFrame) -> pd.DataFrame:
    """Copy the input and create or coerce integer row IDs.

    Args:
        df (pd.DataFrame): Input rows with an optional row_id column.

    Returns:
        pd.DataFrame: A copy with integer row_id values.
    """
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
    """Check whether a line contains any expected header name.

    Args:
        first_line (str): First nonempty input line.
        expected_names (list[str] | None): Header keywords to search for.

    Returns:
        bool: Whether a keyword occurs, ignoring case.
    """
    if not first_line:
        return False
    line = first_line.strip().lower()
    if expected_names:
        return any(name.lower() in line for name in expected_names)
    return False


def normalize_sequence_column(df: pd.DataFrame) -> pd.DataFrame:
    """Identify the sequence column and trim and uppercase its values.

    Args:
        df (pd.DataFrame): Rows containing sequence or one non-row_id column.

    Returns:
        pd.DataFrame: A copy with a normalized sequence column.
    """
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
    """Parse pasted sequence rows with optional headers and IDs.

    Args:
        text (str): Comma-, tab-, or whitespace-separated sequence rows.

    Returns:
        pd.DataFrame: Parsed sequence values and optional row_id values.
    """
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
    """Decode and parse an uploaded sequence file.

    Args:
        file_bytes (bytes): UTF-8 file contents with optional headers and IDs.

    Returns:
        pd.DataFrame: Parsed file columns, or sequence and optional row_id columns.
    """
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
    """Parse uploaded or pasted coordinates, preferring an upload if provided.

    Args:
        uploaded_file (UploadedFile | None): Streamlit upload exposing getvalue().
        pasted_text (str | None): Coordinate text used when no upload is provided.
        site_choice (str): "5ss", "3ss", or "Both"; determines headerless columns.

    Returns:
        pd.DataFrame: Parsed coordinates with an optional row_id column.
    """
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
    """Check required columns and convert coordinates to numbers in place.

    Args:
        df (pd.DataFrame): Coordinate rows to validate and update.
        site_type (str): "5ss", "3ss", or "Both".

    Returns:
        pd.DataFrame: The same dataframe with numeric coordinate columns.
    """
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
