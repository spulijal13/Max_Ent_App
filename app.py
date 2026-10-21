"""Streamlit entry point for the MaxEnt Splice Site Tool."""

import streamlit as st

from maxent.inputs import (
    parse_coordinate_input,
    parse_pasted_sequences,
    parse_uploaded_sequences,
)
from maxent.results import show_results
from maxent.scoring import run_coordinate_mode, run_sequence_mode, validate_required_files


def render_sequence_tab():
    st.subheader("Run from sequences")
    st.write("Paste one sequence per line, or two columns as row_id and sequence. Headers are optional.")

    seq_site_choice = st.radio("Choose sequence type", ["5ss", "3ss"], horizontal=True, key="seq_site_choice")
    seq_input_method = st.radio(
        "Input method", ["Upload file", "Paste text"],
        horizontal=True, key="seq_input_method",
    )
    seq_upload = None
    seq_paste = ""
    if seq_input_method == "Upload file":
        seq_upload = st.file_uploader("Upload sequence file", type=["csv", "txt", "tsv"], key="seq_upload")
    else:
        seq_paste = st.text_area(
            "Paste sequences here",
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

            result_df = run_sequence_mode(seq_df, seq_site_choice)

            show_results(result_df, "seq", f"sequences_{seq_site_choice}_results")

        except Exception as e:
            st.error(str(e))


def render_coordinate_tab():
    st.subheader("Run from raw coordinates")
    st.write("Headers are optional. If there are no headers, use the fixed column order shown below.")

    coord_site_choice = st.radio(
        "Choose what to generate and score",
        ["5ss", "3ss", "Both"],
        horizontal=True,
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

    coord_input_method = st.radio(
        "Input method", ["Upload file", "Paste text"],
        horizontal=True, key="coord_input_method",
    )
    coord_upload = None
    coord_paste = ""
    if coord_input_method == "Upload file":
        coord_upload = st.file_uploader("Upload coordinate file", type=["csv", "txt", "tsv"], key="coord_upload")
    else:
        coord_paste = st.text_area(
            "Paste coordinate rows here",
            height=240,
            key="coord_paste",
            placeholder="Examples:\nchr1 + 100 108\nchr2 - 500 508",
        )

    if st.button("Run from coordinates", key="run_coord"):
        try:
            coord_df = parse_coordinate_input(coord_upload, coord_paste, coord_site_choice)
            result_df = run_coordinate_mode(coord_df, coord_site_choice)

            show_results(result_df, "coord", f"coordinates_{coord_site_choice}_results")

        except Exception as e:
            st.error(str(e))


def main():
    st.set_page_config(page_title="MaxEnt Splice Site Tool", layout="wide")
    st.title("MaxEnt Splice Site Tool")
    validate_required_files()
    tab_seq, tab_coord = st.tabs(["From sequences", "From raw coordinates"])
    with tab_seq:
        render_sequence_tab()
    with tab_coord:
        render_coordinate_tab()


if __name__ == "__main__":
    main()
