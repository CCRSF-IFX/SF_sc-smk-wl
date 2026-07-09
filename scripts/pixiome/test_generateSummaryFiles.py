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
    Path("pixelator").mkdir()
    quality_payload = {
        "x": {
            "data": [
                ["S1_UNT"],
                ["42"],
                ["0.44"],
                ["0.26"],
                ["3.35"],
                ["0.16"],
                ["92.69"],
                ["4.35"],
                ["92.04"],
                ["85.71"],
                ["2.68"],
            ],
            "container": (
                "<table><thead><tr>"
                "<th>Sample ID</th>"
                "<th>Number of cells in sample <i class=\"bi bi-info-circle-fill\"></i></th>"
                "<th>Median isotype % counts</th>"
                "<th>Median proteins per cell [k]</th>"
                "<th>Median reads per cell [k]</th>"
                "<th>Total reads [M]</th>"
                "<th>% Graph reads</th>"
                "<th>% Denoised UMIs</th>"
                "<th>Graph node saturation [%]</th>"
                "<th>Graph edge saturation [%]</th>"
                "<th>Median mean coreness</th>"
                "</tr></thead></table>"
            ),
        }
    }
    Path("pixelator/experiment-summary.html").write_text(
        '<html><script type="application/json" data-for="quality">'
        + json.dumps(quality_payload)
        + "</script></html>"
    )
    Path("params.pixiome.yaml").write_text('pixelator_container: "quay.io/pixelgen-technologies/pixelator:0.27.2"\n')

    module.main()

    workbook = Path("finalreport/metric_summary.xlsx")
    assert workbook.exists()
    assert Path("finalreport/summaries/experiment-summary.html").exists()
    with zipfile.ZipFile(workbook) as zipped:
        workbook_xml = zipped.read("xl/workbook.xml").decode()
        shared_strings = zipped.read("xl/sharedStrings.xml").decode()
    assert "metrics_summary" in workbook_xml
    assert "run_metadata" in workbook_xml
    assert "quality_metrics" not in workbook_xml
    assert "json_metrics" not in workbook_xml
    assert "S1_UNT" in shared_strings
    assert "0.44" in shared_strings
    assert "92.69" in shared_strings
    assert "Median isotype % counts" in shared_strings
    assert "Pixelator" in shared_strings
