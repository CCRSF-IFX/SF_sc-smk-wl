import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path


TREKKER_SCRIPT_DIR = Path(__file__).parents[1] / "scripts" / "trekker"
sys.path.insert(0, str(TREKKER_SCRIPT_DIR))

from trekker_utils import (  # noqa: E402
    TREKKER_COLUMNS,
    TREKKER_LIBRARY_COLUMNS,
    TrekkerConfigError,
    cellranger_fastqs,
    load_libraries,
    load_samplesheet,
    prepare_runtime,
    run_trekker,
    stage_fastqs,
    write_vendor_samplesheet,
)


class TrekkerInputTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.matrix_dir = self.root / "filtered_feature_bc_matrix"
        self.matrix_dir.mkdir()
        for name in ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"):
            (self.matrix_dir / name).touch()
        for name in ("tile_BeadBarcodes.txt", "sample_R1.fastq.gz", "sample_R2.fastq.gz"):
            (self.root / name).touch()
        self.fastq_root = self.root / "FLOWCELL01"
        self.gex_fastq_dir = self.fastq_root / "gex"
        self.trekker_fastq_dir = self.fastq_root / "trekker"
        self.gex_fastq_dir.mkdir(parents=True)
        self.trekker_fastq_dir.mkdir()
        for read in ("R1", "R2"):
            (self.gex_fastq_dir / f"gex_S1_L001_{read}_001.fastq.gz").touch()
            (self.trekker_fastq_dir / f"spatial_S2_L001_{read}_001.fastq.gz").touch()
        self.transcriptome = self.root / "refdata-gex-test"
        self.transcriptome.mkdir()
        (self.transcriptome / "reference.json").touch()

    def tearDown(self):
        self.temporary_directory.cleanup()

    def record(self, sample="sample_1"):
        return {
            "sample": sample,
            "sc_sample": f"{sample}_GEX",
            "experiment_date": "20260924",
            "barcode_file": str((self.root / "tile_BeadBarcodes.txt").resolve()),
            "fastq_1": str((self.root / "sample_R1.fastq.gz").resolve()),
            "fastq_2": str((self.root / "sample_R2.fastq.gz").resolve()),
            "sc_outdir": str(self.matrix_dir.resolve()),
            "sc_platform": "TrekkerU_CX",
            "profile": "Conda",
            "subsample": "NO",
            "cores": "8",
        }

    def write_sheet(self, records, fieldnames=TREKKER_COLUMNS):
        sheet = self.root / "samplesheet.csv"
        with sheet.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)
        return sheet

    def library_rows(self, sample="sample_1"):
        return [
            {
                "Name": sample,
                "Flowcell": "FLOWCELL01",
                "Sample": "gex",
                "Type": "Gene Expression",
                "BarcodeFile": "",
                "Profile": "",
                "Subsample": "",
                "Cores": "",
            },
            {
                "Name": sample,
                "Flowcell": "FLOWCELL01",
                "Sample": "spatial",
                "Type": "Trekker",
                "BarcodeFile": str((self.root / "tile_BeadBarcodes.txt").resolve()),
                "Profile": "Conda",
                "Subsample": "NO",
                "Cores": "8",
            },
        ]

    def write_libraries(self, rows, fieldnames=TREKKER_LIBRARY_COLUMNS):
        libraries = self.root / "libraries.csv"
        with libraries.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        return libraries

    def load_library_records(self, rows=None, unaligned=None):
        return load_libraries(
            self.write_libraries(rows or self.library_rows()),
            unaligned=unaligned or [self.fastq_root],
            transcriptome=self.transcriptome,
            cellranger_output_root=self.root / "cellranger",
            trekker_fastq_root=self.root / "trekker_fastqs",
        )

    def test_loads_and_normalizes_multiple_project_rows(self):
        sheet = self.write_sheet([self.record("sample_1"), self.record("sample_2")])
        records = load_samplesheet(sheet)
        self.assertEqual(list(records), ["sample_1", "sample_2"])
        self.assertEqual(records["sample_1"]["profile"], "conda")
        self.assertEqual(records["sample_1"]["subsample"], "no")

    def test_rejects_non_cx_platform(self):
        record = self.record()
        record["sc_platform"] = "TrekkerU_C"
        sheet = self.write_sheet([record])
        with self.assertRaisesRegex(TrekkerConfigError, "only sc_platform=TrekkerU_CX"):
            load_samplesheet(sheet)

    def test_rejects_missing_filtered_matrix_component(self):
        (self.matrix_dir / "features.tsv.gz").unlink()
        sheet = self.write_sheet([self.record()])
        with self.assertRaisesRegex(TrekkerConfigError, "features.tsv.gz"):
            load_samplesheet(sheet)

    def test_libraries_csv_builds_cellranger_and_trekker_inputs(self):
        record = self.load_library_records()["sample_1"]
        self.assertEqual(record["profile"], "conda")
        self.assertEqual(record["subsample"], "no")
        self.assertEqual(record["experiment_date"], "")
        self.assertEqual(
            record["sc_outdir"],
            str(self.root / "cellranger" / "sample_1" / "outs" / "filtered_feature_bc_matrix"),
        )
        self.assertEqual(len(cellranger_fastqs(record)), 2)
        self.assertEqual(record["cellranger_fastq_samples"], ["gex"])
        self.assertEqual(len(record["trekker_r1_fastqs"]), 1)
        self.assertEqual(
            record["fastq_1"],
            str(self.root / "trekker_fastqs" / "sample_1_R1.fastq.gz"),
        )

    def test_libraries_csv_combines_trekker_rows_from_multiple_runs(self):
        second_root = self.root / "FLOWCELL02"
        second_fastq_dir = second_root / "trekker"
        second_fastq_dir.mkdir(parents=True)
        for read in ("R1", "R2"):
            (second_fastq_dir / f"spatial2_S1_L001_{read}_001.fastq.gz").touch()
        rows = self.library_rows()
        second_trekker_row = dict(rows[1])
        second_trekker_row.update({"Flowcell": "FLOWCELL02", "Sample": "spatial2"})
        rows.append(second_trekker_row)

        record = self.load_library_records(
            rows,
            unaligned=[self.fastq_root, second_root],
        )["sample_1"]
        self.assertEqual(len(record["trekker_r1_fastqs"]), 2)
        self.assertEqual(len(record["trekker_r2_fastqs"]), 2)

    def test_libraries_csv_combines_both_library_types_across_flowcells(self):
        second_root = self.root / "FLOWCELL02"
        second_gex_dir = second_root / "gex"
        second_trekker_dir = second_root / "trekker"
        second_gex_dir.mkdir(parents=True)
        second_trekker_dir.mkdir()
        for read in ("I1", "R1", "R2"):
            (second_gex_dir / f"gex_S1_L001_{read}_001.fastq.gz").touch()
        for read in ("R1", "R2"):
            (second_trekker_dir / f"spatial_S2_L001_{read}_001.fastq.gz").touch()

        rows = self.library_rows()
        second_gex_row = dict(rows[0], Flowcell="FLOWCELL02")
        second_trekker_row = dict(rows[1], Flowcell="FLOWCELL02")
        rows.extend([second_gex_row, second_trekker_row])

        record = self.load_library_records(
            rows,
            unaligned=[self.fastq_root, second_root],
        )["sample_1"]
        self.assertEqual(len(record["cellranger_fastqs"]), 4)
        self.assertEqual(len(record["cellranger_fastq_dirs"]), 2)
        self.assertEqual(record["cellranger_fastq_samples"], ["gex"])
        self.assertEqual(len(record["trekker_r1_fastqs"]), 2)
        self.assertEqual(len(record["trekker_r2_fastqs"]), 2)
        self.assertEqual(len(record["source_fastqs_by_sample"]["gex"]), 5)
        self.assertEqual(len(record["source_fastqs_by_sample"]["spatial"]), 4)

    def test_libraries_csv_rejects_ambiguous_flowcell_substring(self):
        ambiguous_root = self.root / "FLOWCELL010"
        ambiguous_root.mkdir()
        with self.assertRaisesRegex(TrekkerConfigError, "matched 0"):
            self.load_library_records(unaligned=[ambiguous_root])

    def test_libraries_csv_rejects_unconfigured_absolute_flowcell(self):
        rows = self.library_rows()
        rows[0]["Flowcell"] = str(self.fastq_root.resolve())
        with self.assertRaisesRegex(TrekkerConfigError, "not listed in config.unaligned"):
            self.load_library_records(rows, unaligned=[self.root])

    def test_libraries_csv_rejects_prefix_reused_by_another_sample(self):
        rows = self.library_rows()
        rows.append(dict(rows[0], Name="sample_2"))
        with self.assertRaisesRegex(TrekkerConfigError, "already assigned"):
            self.load_library_records(rows)

    def test_libraries_csv_requires_both_library_types(self):
        with self.assertRaisesRegex(TrekkerConfigError, "requires at least one Trekker"):
            self.load_library_records(self.library_rows()[:1])

    def test_libraries_csv_rejects_unpaired_fastqs(self):
        (self.gex_fastq_dir / "gex_S1_L001_R2_001.fastq.gz").unlink()
        with self.assertRaisesRegex(TrekkerConfigError, "paired Gene Expression"):
            self.load_library_records()

    def test_libraries_csv_rejects_unknown_library_type(self):
        rows = self.library_rows()
        rows[1]["Type"] = "Spatial"
        with self.assertRaisesRegex(TrekkerConfigError, "Gene Expression.*Trekker"):
            self.load_library_records(rows)

    def test_libraries_csv_rejects_same_fastqs_for_both_library_types(self):
        rows = self.library_rows()
        rows[1]["Sample"] = "gex"
        with self.assertRaisesRegex(
            TrekkerConfigError,
            "already assigned.*library type",
        ):
            self.load_library_records(rows)

    def test_libraries_csv_rejects_vendor_unsafe_path(self):
        unsafe_barcode = self.root / "tile,unsafe_BeadBarcodes.txt"
        unsafe_barcode.touch()
        rows = self.library_rows()
        rows[1]["BarcodeFile"] = str(unsafe_barcode.resolve())
        with self.assertRaisesRegex(TrekkerConfigError, "cannot contain commas"):
            self.load_library_records(rows)

    def test_libraries_csv_requires_exact_header(self):
        columns = list(TREKKER_LIBRARY_COLUMNS)
        columns[-2], columns[-1] = columns[-1], columns[-2]
        libraries = self.write_libraries(self.library_rows(), columns)
        with self.assertRaisesRegex(TrekkerConfigError, "exact order"):
            load_libraries(
                libraries,
                unaligned=[self.fastq_root],
                transcriptome=self.transcriptome,
                cellranger_output_root=self.root / "cellranger",
                trekker_fastq_root=self.root / "trekker_fastqs",
            )

    def test_writes_one_vendor_row(self):
        record = self.record()
        destination = self.root / "generated" / "sample_1.csv"
        write_vendor_samplesheet(record, destination)
        with destination.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["sample"], "sample_1")
        self.assertEqual(tuple(rows[0]), TREKKER_COLUMNS)

    def test_library_record_writes_vendor_schema_only(self):
        record = self.load_library_records()["sample_1"]
        destination = self.root / "generated" / "sample_1.csv"
        write_vendor_samplesheet(record, destination)
        with destination.open(newline="") as handle:
            reader = csv.DictReader(handle)
            row = next(reader)
            self.assertEqual(tuple(reader.fieldnames), TREKKER_COLUMNS)
        self.assertEqual(row["sc_outdir"], record["sc_outdir"])
        self.assertEqual(row["experiment_date"], "")

    def test_stages_single_and_multiple_flowcell_fastqs_without_overwriting_source(self):
        first = self.root / "first_R1.fastq.gz"
        second = self.root / "second_R1.fastq.gz"
        destination = self.root / "staged" / "sample_R1.fastq.gz"
        with gzip.open(str(first), "wb") as handle:
            handle.write(b"first\n")
        with gzip.open(str(second), "wb") as handle:
            handle.write(b"second\n")
        first_bytes = first.read_bytes()
        second_bytes = second.read_bytes()

        stage_fastqs([first], destination)
        self.assertTrue(destination.is_symlink())
        with gzip.open(str(destination), "rb") as handle:
            self.assertEqual(handle.read(), b"first\n")

        stage_fastqs([first, second], destination)
        self.assertFalse(destination.is_symlink())
        with gzip.open(str(destination), "rb") as handle:
            self.assertEqual(handle.read(), b"first\nsecond\n")
        self.assertEqual(first.read_bytes(), first_bytes)
        self.assertEqual(second.read_bytes(), second_bytes)

        stage_fastqs([second], destination)
        self.assertTrue(destination.is_symlink())
        with gzip.open(str(destination), "rb") as handle:
            self.assertEqual(handle.read(), b"second\n")


class TrekkerRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.installation = self.root / "trekker-v1.4.11"
        self.installation.mkdir()
        (self.installation / "common").mkdir()
        (self.installation / "cellbarcode_whitelists").mkdir()
        (self.installation / "nuclei_locater_toplevel.sh").write_text(
            "#!/usr/bin/env bash\n"
            "set -e\n"
            "OUT_DIR=/home/trekker_out/\n"
            "SCRIPT_DIR=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"\n"
            "source \"$SCRIPT_DIR/nuclei_locater_conda.sh\"\n"
        )
        (self.installation / "nuclei_locater_conda.sh").write_text(
            "#!/usr/bin/env bash\n"
            "PROFILE_CONDA_PATH=/home/tools/miniconda3/envs/trekker/\n"
            "source /home/tools/miniconda3/etc/profile.d/conda.sh\n"
            "conda activate ${PROFILE_CONDA_PATH}\n"
        )
        (self.installation / "nuclei_locater_singularity.sh").touch()
        (self.installation / "trekker-v1.4.11.sif").touch()
        self.environment = self.root / "env"
        (self.environment / "bin").mkdir(parents=True)
        (self.environment / "bin" / "python").touch()
        (self.environment / "bin" / "Rscript").touch()

    def tearDown(self):
        self.temporary_directory.cleanup()

    def test_runtime_overrides_paths_without_editing_installation(self):
        output_root = self.root / "output"
        runtime = self.root / "runtime"
        runtime.mkdir()
        # A previous profile may have left launcher symlinks in the overlay.
        # Reconfiguring for conda must replace them, never follow them.
        (runtime / "nuclei_locater_toplevel.sh").symlink_to(
            self.installation / "nuclei_locater_toplevel.sh"
        )
        (runtime / "nuclei_locater_conda.sh").symlink_to(
            self.installation / "nuclei_locater_conda.sh"
        )
        launcher = prepare_runtime(
            self.installation,
            output_root,
            runtime,
            "conda",
            self.environment,
        )
        launcher_text = launcher.read_text()
        conda_text = (runtime / "nuclei_locater_conda.sh").read_text()
        self.assertIn(f"OUT_DIR={output_root}", launcher_text)
        self.assertIn(f"SCRIPT_DIR={self.installation}", launcher_text)
        self.assertIn(f"source {runtime / 'nuclei_locater_conda.sh'}", launcher_text)
        self.assertIn(f"PROFILE_CONDA_PATH={self.environment}", conda_text)
        self.assertIn("export R_LIBS_SITE=", conda_text)
        self.assertIn("OUT_DIR=/home/trekker_out/", (self.installation / launcher.name).read_text())
        self.assertIn(
            "PROFILE_CONDA_PATH=/home/tools/miniconda3/envs/trekker/",
            (self.installation / "nuclei_locater_conda.sh").read_text(),
        )
        self.assertFalse(launcher.is_symlink())
        self.assertFalse((runtime / "nuclei_locater_conda.sh").is_symlink())
        self.assertFalse((runtime / "common").exists())

    def test_singularity_runtime_loads_module(self):
        launcher = prepare_runtime(
            self.installation,
            self.root / "output",
            self.root / "runtime-singularity",
            "singularity",
        )
        launcher_text = launcher.read_text()
        self.assertIn("set -e\n\nmodule load singularity", launcher_text)
        self.assertIn("command -v singularity", launcher_text)
        self.assertNotIn("module load singularity", (
            self.installation / "nuclei_locater_toplevel.sh"
        ).read_text())

    def test_runner_executes_overlay_and_requires_expected_report(self):
        barcode = self.root / "tile_BeadBarcodes.txt"
        fastq_1 = self.root / "sample_R1.fastq.gz"
        fastq_2 = self.root / "sample_R2.fastq.gz"
        matrix_dir = self.root / "filtered_feature_bc_matrix"
        matrix_dir.mkdir()
        barcode.touch()
        fastq_1.touch()
        fastq_2.touch()
        for name in ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"):
            (matrix_dir / name).touch()

        # Extend the fake vendor launcher so its final action creates the report
        # in the same layout as Trekker v1.4.11.
        with (self.installation / "nuclei_locater_toplevel.sh").open("a") as handle:
            handle.write(
                "SAMPLE_DATA=$(tail -n +2 \"$1\" | head -n 1)\n"
                "SAMPLE_ID=$(echo \"$SAMPLE_DATA\" | awk -F ',' '{print $1}')\n"
                "ANALYSIS_DATE=$(echo \"$SAMPLE_DATA\" | awk -F ',' '{print $3}')\n"
                "REPORT_DIR=\"${OUT_DIR}/${ANALYSIS_DATE}_${SAMPLE_ID}/"
                "trekker_${SAMPLE_ID}/output\"\n"
                "mkdir -p \"$REPORT_DIR\"\n"
                "touch \"$REPORT_DIR/${SAMPLE_ID}_Trekker_Report.html\"\n"
            )

        record = {
            "sample": "sample_1",
            "sc_sample": "sample_1_GEX",
            "experiment_date": "",
            "barcode_file": str(barcode),
            "fastq_1": str(fastq_1),
            "fastq_2": str(fastq_2),
            "sc_outdir": str(matrix_dir),
            "sc_platform": "TrekkerU_CX",
            "profile": "conda",
            "subsample": "no",
            "cores": "8",
        }
        samplesheet = self.root / "samplesheet.csv"
        write_vendor_samplesheet(record, samplesheet)
        report = run_trekker(
            samplesheet,
            self.installation,
            self.root / "output",
            self.root / "runtime",
            self.environment,
        )
        self.assertTrue(report.is_file())


if __name__ == "__main__":
    unittest.main()
