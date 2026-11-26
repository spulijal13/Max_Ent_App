"""Streamlit entry point for the MaxEnt Splice Site Tool."""

import hashlib

import streamlit as st

from maxent.inputs import (
    parse_coordinate_input,
    parse_pasted_sequences,
    parse_uploaded_sequences,
)
from maxent.results import show_results
from maxent.scoring import run_coordinate_mode, run_sequence_mode, validate_required_files


def site_label(site):
    return site


def input_signature(site, method, upload, text, strand=None):
    content = upload.getvalue() if upload is not None else text.encode("utf-8")
    return (site, method, hashlib.sha256(content).hexdigest(), strand)


def show_saved_results(key, signature):
    saved = st.session_state.get(key)
    if saved is None:
        return
    if saved["signature"] != signature:
        st.info("Inputs have changed. Showing the last successful run; submit the current inputs to update.")
    st.subheader("Results", anchor=False)
    st.write(f"Last successful run: {saved['label']}")
    with st.expander("Input used for these results"):
        st.dataframe(saved["input"], use_container_width=True)
    show_results(saved["results"], key, saved["filename"])


def render_sequence_tab():
    """Render sequence inputs and handle scoring and result display.

    Args:
        None.

    Returns:
        None: Updates the Streamlit page.
    """
    st.header("Sequences", anchor=False)

    seq_site_choice = st.radio("Site type", ["5ss", "3ss"], format_func=site_label, horizontal=True, key="seq_site_choice")
    seq_input_method = st.radio(
        "Input method", ["Upload file", "Paste text"],
        horizontal=True, key="seq_input_method",
    )
    example = "CAGGTAAGT" if seq_site_choice == "5ss" else "TTCCAAACGAACTTTTGTAGGGA"
    seq_upload = None
    seq_paste = ""
    if seq_input_method == "Upload file":
        seq_upload = st.file_uploader("Upload sequence file", type=["csv", "txt", "tsv"], key="seq_upload")
    else:
        seq_paste = st.text_area(
            "Paste sequences here",
            height=240,
            key="seq_paste",
            placeholder=example,
        )

    with st.expander("Formats & Details"):
        st.write("5ss = 5′ donor site; 3ss = 3′ acceptor site.")
        st.subheader("Sequence format", anchor=False)
        st.write("Use 9 bases for 5ss or 23 bases for 3ss. One sequence per line; A, C, G, T only.")
        st.code(example, language=None)
        st.subheader("Optional row IDs", anchor=False)
        st.code(f"row_id,sequence\n1,{example}", language=None)
        st.write("Use two columns: row_id, then sequence. Headers are optional. "
                 "IDs must be unique integers. FASTA records are not supported.")

    signature = input_signature(seq_site_choice, seq_input_method, seq_upload, seq_paste)
    if st.button("Score sequences", key="run_seq", type="primary"):
        try:
            if seq_upload is not None:
                seq_df = parse_uploaded_sequences(seq_upload.getvalue())
            elif seq_paste and seq_paste.strip():
                seq_df = parse_pasted_sequences(seq_paste)
            else:
                raise ValueError("Provide either an uploaded file or pasted sequences.")

            with st.spinner("Scoring sequences…"):
                result_df = run_sequence_mode(seq_df, seq_site_choice)

            st.session_state["seq_results"] = {
                "signature": signature, "input": seq_df.copy(), "results": result_df,
                "label": f"Sequences / {seq_site_choice}",
                "filename": f"sequences_{seq_site_choice}_results",
            }

        except ValueError as e:
            st.error(f"Check your input: {e}")
        except Exception as e:
            st.error("Processing failed. See technical details below.")
            with st.expander("Technical details"):
                st.code(str(e), language=None)

    show_saved_results("seq_results", signature)


def render_coordinate_tab():
    """Render coordinate inputs and handle extraction and scoring.

    Args:
        None.

    Returns:
        None: Updates the Streamlit page.
    """
    st.header("Coordinates", anchor=False)
    st.caption("HG38 (GRCh38) · 1-based · Inclusive")

    coord_site_choice = st.radio(
        "Site type",
        ["5ss", "3ss", "Both"],
        horizontal=True,
        key="coord_site_choice",
        format_func=site_label,
    )

    examples = {
        "5ss": ("Chromosome   Strand   Start      End", "chr12        -        25215431   25215439"),
        "3ss": ("Chromosome   Strand   Start      End", "chr12        -        25215431   25215453"),
        "Both": (
            "Chromosome   Strand   5ss start   5ss end    3ss start   3ss end",
            "chr12        -        25215431    25215439   25215431    25215453",
        ),
    }
    columns, example = examples[coord_site_choice]
    ucsc_example = "chr12:25,215,431-" + ("25,215,439" if coord_site_choice == "5ss" else "25,215,453")
    ucsc_strand = None
    if coord_site_choice != "Both":
        ucsc_strand = st.radio(
            "UCSC Genome Browser strand (Optional)", [None, "+", "-"],
            format_func=lambda value: {None: "None", "+": "+ (positive)", "-": "− (negative)"}[value],
            horizontal=True, key="ucsc_strand_choice",
            help="Choose positive or negative for UCSC positions. Column-format rows supply their own strand.",
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
            placeholder=(" ".join(example.split()) + (f"\n\nOr UCSC format:\n{ucsc_example}" if coord_site_choice != "Both" else "")),
        )

    with st.expander("Formats & Details"):
        st.write("5ss = 5′ donor site; 3ss = 3′ acceptor site.")
        st.subheader("Column format", anchor=False)
        st.text(columns)
        st.code(example, language=None)
        st.write("One row per site. Separate columns with spaces, tabs, or commas. "
                 "An optional unique integer row ID goes first.")
        st.write("Strand: + or 1 = positive; - or -1 = negative.")
        if coord_site_choice != "Both":
            st.subheader("UCSC browser position", anchor=False)
            st.code(ucsc_example, language=None)
            st.write("Choose + or - under UCSC Genome Browser strand for this format. Choose None for column format; "
                     "each row supplies its own strand. Use one format per run.")
        else:
            st.write("Both uses separate donor and acceptor intervals in each row.")
        sites = ["5ss", "3ss"] if coord_site_choice == "Both" else [coord_site_choice]
        headers = ["Chromosome.scaffold.name", "Strand"]
        for site in sites:
            headers.extend([f"X{site}_start", f"X{site}_end"])
        st.subheader("Optional file headers", anchor=False)
        st.code(",".join(headers), language=None)
        st.caption("BED files use 0-based starts, add 1 to the start before entering coordinates.")

    signature = input_signature(coord_site_choice, coord_input_method, coord_upload, coord_paste, ucsc_strand)
    if st.button("Extract & score", key="run_coord", type="primary"):
        try:
            coord_df = parse_coordinate_input(coord_upload, coord_paste, coord_site_choice, ucsc_strand)
            with st.spinner("Extracting HG38 sequences and scoring…"):
                result_df = run_coordinate_mode(coord_df, coord_site_choice)

            st.session_state["coord_results"] = {
                "signature": signature, "input": coord_df.copy(), "results": result_df,
                "label": f"Coordinates / HG38 / {coord_site_choice}",
                "filename": f"coordinates_{coord_site_choice}_results",
            }

        except ValueError as e:
            st.error(f"Check your input: {e}")
        except Exception as e:
            st.error("Processing failed. See technical details below.")
            with st.expander("Technical details"):
                st.code(str(e), language=None)

    show_saved_results("coord_results", signature)


def main():
    """Configure the page, check required files, and render both tabs.

    Args:
        None.

    Returns:
        None: Updates the Streamlit page.
    """
    st.set_page_config(page_title="MaxEnt Splice Site Tool", layout="wide")
    st.title("MaxEnt Splice Site Tool", anchor=False)
    st.caption(
        "Built on MaxEntScan, originally developed by Gene Yeo and Christopher B. Burge "
        "in the [Burge Lab](https://www.genes.mit.edu/)."
    )
    try:
        validate_required_files()
    except FileNotFoundError as error:
        st.error("The scoring files are missing. Restore the scripts folder before running the app.")
        with st.expander("Technical details"):
            st.code(str(error), language=None)
        st.stop()
    tab_seq, tab_coord = st.tabs(["Sequences", "Coordinates"])
    with tab_seq:
        render_sequence_tab()
    with tab_coord:
        render_coordinate_tab()

    st.divider()
    with st.expander("References"):
        st.header("References", anchor=False)
        st.subheader("MaxEntScan", anchor=False)
        st.markdown(
            "[Burge Lab](https://www.genes.mit.edu/)"
        )
        st.subheader("Research paper", anchor=False)
        st.markdown(
            "**Yeo G, Burge CB (2004).** "
            "[Maximum entropy modeling of short sequence motifs with applications to RNA splicing signals.]"
            "(https://pubmed.ncbi.nlm.nih.gov/15285897/) "
            "*Journal of Computational Biology*, 11(2–3), 377–394. "
            "[doi:10.1089/1066527041410418](https://doi.org/10.1089/1066527041410418)"
        )


if __name__ == "__main__":
    main()
