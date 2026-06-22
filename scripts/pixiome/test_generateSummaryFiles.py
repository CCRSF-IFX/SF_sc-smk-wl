import importlib.util
import json
import os
import zipfile
from pathlib import Path


def load_module():
    path = Path(__file__).with_name("generateSummaryFiles.py")
    spec = importlib.util.spec_from_file_location("pixiome_summary", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_generate_summary_files(tmp_path):
    module = load_module()
    os.chdir(tmp_path)
    Path("pixelator/analysis").mkdir(parents=True)
    Path("pixelator/experiment-summary.html").write_text("<html>summary</html>")
    Path("samplesheet.pixiome.csv").write_text(
        "sample,sample_alias,condition,design,panel,fastq_1,fastq_2\n"
        "1_UNT,S1_UNT,UNT,proxiome-v1,proxiome-v1-immuno-155-v1.1,a_R1.fastq.gz,a_R2.fastq.gz\n"
    )
    Path("params.pixiome.yaml").write_text('pixelator_container: "quay.io/pixelgen-technologies/pixelator:0.27.2"\n')
    report = {
        "sample_id": "1_UNT",
        "product_id": "single-cell-pna",
        "report_type": "analysis",
        "k_cores": {"median_average_k_core": 3.42},
        "svd": {"median_variance_explained_3d": 0.11},
    }
    Path("pixelator/analysis/1_UNT.report.json").write_text(json.dumps(report))

    module.main()

    workbook = Path("finalreport/metric_summary.xlsx")
    assert workbook.exists()
    assert Path("finalreport/summaries/experiment-summary.html").exists()
    with zipfile.ZipFile(workbook) as zipped:
        shared_strings = zipped.read("xl/sharedStrings.xml").decode()
    assert "1_UNT" in shared_strings
    assert "Pixelator" in shared_strings
