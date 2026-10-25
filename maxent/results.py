"""Streamlit result tables, copy-friendly text, and downloads."""

import pandas as pd
import streamlit as st


def make_copy_text(df: pd.DataFrame, columns):
    """Serialize the requested available columns as tab-separated text.

    Args:
        df (pd.DataFrame): Results to serialize.
        columns (list[str]): Column names in the desired output order.

    Returns:
        str: Text including headers, or an empty string if no columns match.
    """
    cols = [c for c in columns if c in df.columns]
    if not cols:
        return ""
    return df[cols].to_csv(sep="\t", index=False)


def show_copy_boxes(result_df: pd.DataFrame, prefix: str):
    """Render available sequence-score blocks with built-in copy buttons.

    Args:
        result_df (pd.DataFrame): Results containing sequence and score columns.
        prefix (str): Legacy output identifier; currently unused.

    Returns:
        None: Adds labeled text blocks to the Streamlit page.
    """
    outputs = [
        ("Sequence and score", ["sequence", "score"]),
        ("5ss sequence and score", ["five_ss_sequence", "five_ss_score"]),
        ("3ss sequence and score", ["three_ss_sequence", "three_ss_score"]),
        (
            "Both 5ss and 3ss sequence-score pairs",
            ["five_ss_sequence", "five_ss_score", "three_ss_sequence", "three_ss_score"],
        ),
    ]
    for label, columns in outputs:
        if set(columns).issubset(result_df.columns):
            st.caption(label)
            st.code(make_copy_text(result_df, columns), language=None)


def show_results(result_df: pd.DataFrame, prefix: str, file_stem: str):
    """Render the result table, adjacent download buttons, and copy blocks.

    Args:
        result_df (pd.DataFrame): Full results to display and export.
        prefix (str): Output identifier forwarded to show_copy_boxes.
        file_stem (str): Download filename without the extension.

    Returns:
        None: Adds output controls to the Streamlit page.
    """
    st.success("Done.")
    st.dataframe(result_df, use_container_width=True)
    with st.container(horizontal=True, gap="small"):
        st.download_button(
            "Download results as TXT",
            data=result_df.to_csv(sep="\t", index=False).encode("utf-8"),
            file_name=f"{file_stem}.txt",
            mime="text/plain",
        )
        st.download_button(
            "Download results as CSV",
            data=result_df.to_csv(index=False).encode("utf-8"),
            file_name=f"{file_stem}.csv",
            mime="text/csv",
        )
    show_copy_boxes(result_df, prefix)
