import importlib.util
from pathlib import Path


def load_prepare_module():
    path = Path(__file__).with_name("prepare_samplesheet.py")
    spec = importlib.util.spec_from_file_location("prepare_samplesheet", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_validate_accepts_panel_file_without_panel(tmp_path):
    samplesheet = tmp_path / "samplesheet.csv"
    samplesheet.write_text(
        "sample,sample_alias,condition,design,panel_file,fastq_1,fastq_2\n"
        "sample1,S1,treated,proxiome-v1,/tmp/panel.csv,R1.fastq.gz,R2.fastq.gz\n"
    )

    load_prepare_module().validate(samplesheet)
