import csv
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


TREKKER_SCRIPT_DIR = Path(__file__).parents[1] / "scripts" / "trekker"
sys.path.insert(0, str(TREKKER_SCRIPT_DIR))

from generateSummaryFiles import generate_summary  # noqa: E402


class TrekkerSummaryTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.libraries = self.root / "libraries.csv"
        self.libraries.write_text(
            "Name,Flowcell,Sample,Type,BarcodeFile,Profile,Subsample,Cores\n"
            "sample1,FLOWCELL,gex,Gene Expression,,,,\n"
            "sample1,FLOWCELL,spatial,Trekker,/tile.txt,conda,no,8\n"
        )
        self.cellranger = self.root / "cellranger" / "sample1" / "outs"
        self.cellranger.mkdir(parents=True)
        with (self.cellranger / "metrics_summary.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["Estimated Number of Cells", "Valid Barcodes"])
            writer.writerow(["1,234", "98.5%"])
        (self.cellranger / "web_summary.html").write_text("cellranger")

        self.trekker = (
            self.root / "trekker_out" / "_sample1" / "trekker_sample1" / "output"
        )
        self.trekker.mkdir(parents=True)
        (self.trekker / "sample1_summary_metrics.csv").write_text(
            "Metrics,Value\n"
            "Sample_ID,sample1\n"
            "Tile_ID,U0102_001\n"
            "Nuclei_confidently_positioned,321\n"
            "Pct_useful_reads,45.67\n"
        )
        (self.trekker / "sample1_Trekker_Report.html").write_text("trekker")

    def tearDown(self):
        self.temporary_directory.cleanup()

    def test_generates_workbook_and_copies_both_reports(self):
        destination = self.root / "finalreport"
        workbook = generate_summary(
            self.libraries,
            self.root / "cellranger",
            self.root / "trekker_out",
            destination,
        )

        metrics = pd.read_excel(workbook, sheet_name="metrics_summary")
        self.assertEqual(metrics.loc[0, "Sample"], "sample1")
        self.assertEqual(metrics.loc[0, "Estimated Number of Cells"], 1234)
        self.assertEqual(metrics.loc[0, "Tile_ID"], "U0102_001")
        self.assertEqual(metrics.loc[0, "Nuclei_confidently_positioned"], 321)
        self.assertTrue((destination / "summaries" / "sample1_web_summary.html").is_file())
        self.assertTrue(
            (destination / "summaries" / "sample1_Trekker_Report.html").is_file()
        )


if __name__ == "__main__":
    unittest.main()
