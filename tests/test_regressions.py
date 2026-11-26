import io
import shutil
import unittest

import pandas as pd

from maxent.inputs import ensure_row_id, parse_coordinate_input, validate_coordinate_columns
from maxent.scoring import run_coordinate_mode, run_sequence_mode


class InputTests(unittest.TestCase):
    def test_strands(self):
        for value, expected in [("+", "+"), ("1", "+"), ("-", "-"), ("-1", "-")]:
            frame = parse_coordinate_input(None, f"chr12 {value} 25215431 25215439", "5ss")
            self.assertEqual(validate_coordinate_columns(frame, "5ss").Strand.iloc[0], expected)

    def test_invalid_intervals_and_strands(self):
        for row in ["chr12 + 0 8", "chr12 + 10 9", "chr12 + 1 9.5",
                    "chr12 + 1 10", "chr12 unknown 1 9", "chr12 + nan 9"]:
            with self.subTest(row=row), self.assertRaises(ValueError):
                validate_coordinate_columns(parse_coordinate_input(None, row, "5ss"), "5ss")

    def test_ids(self):
        for ids in [[1, 1], [1.5, 2], ["a", 2], [float("inf"), 2]]:
            with self.subTest(ids=ids), self.assertRaises(ValueError):
                ensure_row_id(pd.DataFrame({"row_id": ids}))
        self.assertEqual(ensure_row_id(pd.DataFrame({"sequence": ["A", "C"]})).row_id.tolist(), [0, 1])

    def test_known_scores(self):
        for site, sequences, expected in [
            ("5ss", ["caggtaagt", "AAGGTAAGT"], [10.86, 11.00]),
            ("3ss", ["TTCCAAACGAACTTTTGTAGGGA", "TTCCAAACGAACTTTTGCAGGGA"], [2.89, 3.84]),
        ]:
            self.assertEqual(run_sequence_mode(pd.DataFrame({"sequence": sequences}), site).score.tolist(), expected)


class UCSCInputTests(unittest.TestCase):
    def test_ucsc_matches_tabular_input(self):
        for site, end in [("5ss", 25215439), ("3ss", 25215453)]:
            for strand in ["+", "-", "1", "-1"]:
                expected = validate_coordinate_columns(parse_coordinate_input(
                    None, f"chr12 {strand} 25215431 {end}", site), site)
                for commas in [False, True]:
                    position = (f"chr12:{25215431:,}-{end:,}" if commas
                                else f"chr12:25215431-{end}")
                    for uploaded in [False, True]:
                        with self.subTest(site=site, strand=strand, commas=commas, uploaded=uploaded):
                            upload = io.BytesIO(position.encode()) if uploaded else None
                            actual = validate_coordinate_columns(parse_coordinate_input(
                                upload, "" if uploaded else position, site, strand), site)
                            pd.testing.assert_frame_equal(actual, expected)

    def test_ucsc_multiple_rows_and_bom_upload(self):
        text = "\ufeffchr12:25,215,431-25,215,439\r\n\r\nchr1:100-108\r\n"
        actual = ensure_row_id(validate_coordinate_columns(parse_coordinate_input(
            io.BytesIO(text.encode()), "", "5ss", "-"), "5ss"))
        self.assertEqual(actual.row_id.tolist(), [0, 1])
        self.assertEqual(actual.X5ss_start.tolist(), [25215431, 100])
        self.assertEqual(actual.Strand.tolist(), ["-", "-"])

    def test_ucsc_requires_valid_strand(self):
        for strand in [None, "None", "", "unknown", "0"]:
            with self.subTest(strand=strand), self.assertRaises(ValueError):
                parse_coordinate_input(None, "chr12:25215431-25215439", "5ss", strand)

    def test_ucsc_rejects_malformed_or_mixed_rows(self):
        for text in [
            "chr12:25,21,431-25,215,439", "chr12:25215431.5-25215439",
            "chr12:25215431", "chr12:25215431-25215439 +",
            "chr12:25215431-25215439\nchr12 - 25215431 25215439",
        ]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_coordinate_input(None, text, "5ss", "+")

    def test_ucsc_validates_interval_lengths_and_order(self):
        for site, position in [
            ("5ss", "chr12:0-8"), ("5ss", "chr12:108-100"),
            ("5ss", "chr12:100-109"), ("3ss", "chr12:100-108"),
        ]:
            with self.subTest(site=site, position=position), self.assertRaises(ValueError):
                validate_coordinate_columns(parse_coordinate_input(None, position, site, "+"), site)

    def test_ucsc_both_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "one site type"):
            parse_coordinate_input(None, "chr12:25215431-25215439", "Both", "+")

    def test_ucsc_selector_does_not_override_tabular_strand(self):
        for selection in [None, "+", "-"]:
            actual = validate_coordinate_columns(parse_coordinate_input(
                None, "chr12 + 25215431 25215439\nchr12 - 25215431 25215439",
                "5ss", selection), "5ss")
            self.assertEqual(actual.Strand.tolist(), ["+", "-"])


@unittest.skipUnless(shutil.which("Rscript"), "Rscript required")
class GenomeTests(unittest.TestCase):
    def test_kras_opposite_strands(self):
        frame = parse_coordinate_input(None, "chr12 + 25215431 25215439\nchr12 - 25215431 25215439", "5ss")
        result = run_coordinate_mode(frame, "5ss")
        positive, negative = result.sequence.tolist()
        self.assertEqual(positive.translate(str.maketrans("ACGT", "TGCA"))[::-1], negative)
        self.assertNotEqual(positive, negative)
        direct = run_sequence_mode(result[["row_id", "sequence"]], "5ss")
        self.assertEqual(result.score.tolist(), direct.score.tolist())

    def test_genome_bounds(self):
        for row in ["chrUnknown + 1 9", "chr12 + 999999991 999999999"]:
            with self.subTest(row=row), self.assertRaises(RuntimeError):
                run_coordinate_mode(parse_coordinate_input(None, row, "5ss"), "5ss")


if __name__ == "__main__":
    unittest.main()
