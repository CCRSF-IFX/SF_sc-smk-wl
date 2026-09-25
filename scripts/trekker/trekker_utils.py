#!/usr/bin/env python3
"""Validation and runtime helpers for the Takara Trekker pipeline."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import shutil
import subprocess
import tempfile
from collections import OrderedDict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence


TREKKER_COLUMNS = (
    "sample",
    "sc_sample",
    "experiment_date",
    "barcode_file",
    "fastq_1",
    "fastq_2",
    "sc_outdir",
    "sc_platform",
    "profile",
    "subsample",
    "cores",
)
TREKKER_LIBRARY_COLUMNS = (
    "Name",
    "Flowcell",
    "Sample",
    "Type",
    "BarcodeFile",
    "Profile",
    "Subsample",
    "Cores",
)
SUPPORTED_PLATFORM = "TrekkerU_CX"
SUPPORTED_PROFILES = {"conda", "singularity", "docker"}
MATRIX_FILES = ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
SAFE_SAMPLE_RE = re.compile(r"^[A-Za-z0-9_-]+$")
SUPPORTED_LIBRARY_TYPES = {"Gene Expression", "Trekker"}


class TrekkerConfigError(ValueError):
    """Raised when a Trekker samplesheet or installation is invalid."""


def _require_absolute_existing_path(
    value: str, field: str, sample: str, *, directory: bool = False
) -> str:
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise TrekkerConfigError(
            f"Sample '{sample}': '{field}' must be an absolute path: {value!r}"
        )
    if directory and not path.is_dir():
        raise TrekkerConfigError(
            f"Sample '{sample}': '{field}' directory does not exist: {path}"
        )
    if not directory and not path.is_file():
        raise TrekkerConfigError(
            f"Sample '{sample}': '{field}' file does not exist: {path}"
        )
    return str(path.resolve())


def _discover_fastq_pair(
    fastq_root: Path,
    fastq_sample: str,
    final_sample: str,
    library_type: str,
) -> Dict[str, List[str]]:
    """Find and validate one demultiplexed R1/R2 FASTQ set."""

    reads: Dict[str, List[Path]] = {
        "R1": [],
        "R2": [],
        "R3": [],
        "I1": [],
        "I2": [],
    }
    pattern = re.compile(
        rf"^{re.escape(fastq_sample)}(?:_S[0-9]+)?_L[0-9]{{3}}_"
        rf"(R[123]|I[12])_[0-9]{{3}}\.fastq\.gz$"
    )
    for path in fastq_root.rglob("*.fastq.gz"):
        match = pattern.fullmatch(path.name)
        if match:
            reads[match.group(1)].append(path.resolve())

    if not reads["R1"] or not reads["R2"]:
        raise TrekkerConfigError(
            f"Sample '{final_sample}': no paired {library_type} R1/R2 FASTQs were "
            f"found for prefix '{fastq_sample}' below {fastq_root}"
        )

    def pair_key(path: Path) -> str:
        name = re.sub(r"_R[12]_([0-9]{3}\.fastq\.gz)$", r"_R?_\1", path.name)
        return str(path.parent / name)

    r1_keys = {pair_key(path) for path in reads["R1"]}
    r2_keys = {pair_key(path) for path in reads["R2"]}
    if r1_keys != r2_keys:
        missing_r1 = sorted(r2_keys - r1_keys)
        missing_r2 = sorted(r1_keys - r2_keys)
        details = []
        if missing_r1:
            details.append("missing R1 for " + ", ".join(missing_r1))
        if missing_r2:
            details.append("missing R2 for " + ", ".join(missing_r2))
        raise TrekkerConfigError(
            f"Sample '{final_sample}': unpaired {library_type} FASTQs: "
            + "; ".join(details)
        )
    return {
        read: [str(path) for path in sorted(paths)]
        for read, paths in reads.items()
        if paths
    }


def _resolve_fastq_root(
    flowcell: str,
    unaligned: Sequence[os.PathLike],
    final_sample: str,
) -> Path:
    """Resolve a multi-style Flowcell value against configured FASTQ roots."""

    configured_roots = [Path(path).expanduser().resolve() for path in unaligned]
    candidate = Path(flowcell).expanduser()
    if candidate.is_absolute():
        if not candidate.is_dir():
            raise TrekkerConfigError(
                f"Sample '{final_sample}': Flowcell path does not exist: {candidate}"
            )
        candidate = candidate.resolve()
        if candidate not in configured_roots:
            raise TrekkerConfigError(
                f"Sample '{final_sample}': absolute Flowcell path is not listed in "
                f"config.unaligned: {candidate}"
            )
        return candidate

    # Match a complete path component or a token in an Illumina-style run name.
    # Plain substring matching can silently resolve FC1 to an FC10 directory.
    token = re.compile(rf"(?:^|[_-]){re.escape(flowcell)}(?:$|[_-])")
    matches = [
        root
        for root in configured_roots
        if any(component == flowcell or token.search(component) for component in root.parts)
    ]
    if len(matches) != 1:
        raise TrekkerConfigError(
            f"Sample '{final_sample}': Flowcell '{flowcell}' must match exactly one "
            f"configured FASTQ root; matched {len(matches)}"
        )
    if not matches[0].is_dir():
        raise TrekkerConfigError(
            f"Sample '{final_sample}': FASTQ root does not exist: {matches[0]}"
        )
    return matches[0]


def _unique(values: Iterable[str]) -> List[str]:
    return list(OrderedDict.fromkeys(values))


def load_libraries(
    libraries: os.PathLike,
    *,
    unaligned: Sequence[os.PathLike],
    transcriptome: os.PathLike,
    cellranger_output_root: os.PathLike,
    trekker_fastq_root: os.PathLike,
) -> "OrderedDict[str, Dict[str, Any]]":
    """Load a multi-style libraries.csv and build internal Trekker records."""

    path = Path(libraries).expanduser().resolve()
    if not path.is_file():
        raise TrekkerConfigError(f"Trekker libraries CSV does not exist: {path}")

    transcriptome_path = Path(transcriptome).expanduser().resolve()
    if not transcriptome_path.is_dir():
        raise TrekkerConfigError(
            f"Cell Ranger transcriptome directory does not exist: {transcriptome_path}"
        )
    if not (transcriptome_path / "reference.json").is_file():
        raise TrekkerConfigError(
            "Cell Ranger transcriptome is missing reference.json: "
            f"{transcriptome_path}"
        )

    configured_roots = [Path(item).expanduser().resolve() for item in unaligned]
    if not configured_roots:
        raise TrekkerConfigError("At least one configured FASTQ root is required")
    if len(configured_roots) != len(set(configured_roots)):
        raise TrekkerConfigError("Configured FASTQ roots must be unique")

    groups: "OrderedDict[str, Dict[str, List[Dict[str, Any]]]]" = OrderedDict()
    source_sample_owners: Dict[str, tuple] = {}
    used_roots = set()
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        actual_columns = tuple(reader.fieldnames or ())
        if actual_columns != TREKKER_LIBRARY_COLUMNS:
            raise TrekkerConfigError(
                "Trekker libraries.csv columns must appear in this exact order: "
                + ",".join(TREKKER_LIBRARY_COLUMNS)
                + f". Received: {','.join(actual_columns)}"
            )

        for row_number, raw_row in enumerate(reader, start=2):
            row = {key: (value or "").strip() for key, value in raw_row.items()}
            for column in ("Name", "Flowcell", "Sample", "Type"):
                if not row[column]:
                    raise TrekkerConfigError(
                        f"Row {row_number}: column '{column}' is required"
                    )
            final_sample = row["Name"]
            fastq_sample = row["Sample"]
            for field, value in (("Name", final_sample), ("Sample", fastq_sample)):
                if not SAFE_SAMPLE_RE.fullmatch(value):
                    raise TrekkerConfigError(
                        f"Row {row_number}: '{field}' must contain only letters, "
                        "numbers, underscores, and hyphens"
                    )
            if row["Type"] not in SUPPORTED_LIBRARY_TYPES:
                raise TrekkerConfigError(
                    f"Row {row_number}: Type must be 'Gene Expression' or 'Trekker'; "
                    f"received {row['Type']!r}"
                )

            owner = (final_sample, row["Type"])
            previous_owner = source_sample_owners.setdefault(fastq_sample, owner)
            if previous_owner != owner:
                raise TrekkerConfigError(
                    f"Row {row_number}: FASTQ Sample prefix '{fastq_sample}' is "
                    f"already assigned to Name '{previous_owner[0]}' as "
                    f"'{previous_owner[1]}'; source prefixes must identify one "
                    "final sample and library type across all flowcells"
                )

            fastq_root = _resolve_fastq_root(row["Flowcell"], unaligned, final_sample)
            used_roots.add(fastq_root)
            row["fastq_root"] = str(fastq_root)
            row["fastqs"] = _discover_fastq_pair(
                fastq_root,
                fastq_sample,
                final_sample,
                row["Type"],
            )
            group = groups.setdefault(
                final_sample,
                {"Gene Expression": [], "Trekker": []},
            )
            group[row["Type"]].append(row)

    if not groups:
        raise TrekkerConfigError("Trekker libraries.csv must contain at least one row")

    unused_roots = set(configured_roots) - used_roots
    if unused_roots:
        raise TrekkerConfigError(
            "Not all configured FASTQ roots are represented in libraries.csv: "
            + ", ".join(str(item) for item in sorted(unused_roots))
        )

    cellranger_root = Path(cellranger_output_root).expanduser().resolve()
    concat_root = Path(trekker_fastq_root).expanduser().resolve()
    records: "OrderedDict[str, Dict[str, Any]]" = OrderedDict()
    for final_sample, group in groups.items():
        gex_rows = group["Gene Expression"]
        trekker_rows = group["Trekker"]
        if not gex_rows or not trekker_rows:
            missing = "Gene Expression" if not gex_rows else "Trekker"
            raise TrekkerConfigError(
                f"Sample '{final_sample}' requires at least one {missing} row"
            )

        metadata_fields = (
            "BarcodeFile",
            "Profile",
            "Subsample",
            "Cores",
        )
        for row in trekker_rows:
            for field in metadata_fields:
                if not row[field]:
                    raise TrekkerConfigError(
                        f"Sample '{final_sample}': '{field}' is required on every "
                        "Trekker row"
                    )
        metadata = {field: trekker_rows[0][field] for field in metadata_fields}
        for row in trekker_rows[1:]:
            for field, expected in metadata.items():
                if row[field] != expected:
                    raise TrekkerConfigError(
                        f"Sample '{final_sample}': all Trekker rows must use the same "
                        f"'{field}' value"
                    )

        profile = metadata["Profile"].lower()
        if profile not in SUPPORTED_PROFILES:
            raise TrekkerConfigError(
                f"Sample '{final_sample}': Profile must be one of "
                f"{sorted(SUPPORTED_PROFILES)}"
            )
        subsample = metadata["Subsample"].lower()
        if subsample not in {"yes", "no"}:
            raise TrekkerConfigError(
                f"Sample '{final_sample}': Subsample must be 'yes' or 'no'"
            )
        try:
            cores = int(metadata["Cores"])
        except ValueError as exc:
            raise TrekkerConfigError(
                f"Sample '{final_sample}': Cores must be a positive integer"
            ) from exc
        if cores < 1:
            raise TrekkerConfigError(
                f"Sample '{final_sample}': Cores must be a positive integer"
            )
        barcode_file = _require_absolute_existing_path(
            metadata["BarcodeFile"], "BarcodeFile", final_sample
        )

        gex_fastqs = [
            fastq
            for row in gex_rows
            for read in ("R1", "R2")
            for fastq in row["fastqs"][read]
        ]
        trekker_r1 = [fastq for row in trekker_rows for fastq in row["fastqs"]["R1"]]
        trekker_r2 = [fastq for row in trekker_rows for fastq in row["fastqs"]["R2"]]
        for label, files in (
            ("Gene Expression", gex_fastqs),
            ("Trekker R1", trekker_r1),
            ("Trekker R2", trekker_r2),
        ):
            if len(files) != len(set(files)):
                raise TrekkerConfigError(
                    f"Sample '{final_sample}': duplicate {label} FASTQ input detected"
                )
        shared_fastqs = set(gex_fastqs) & set(trekker_r1 + trekker_r2)
        if shared_fastqs:
            raise TrekkerConfigError(
                f"Sample '{final_sample}': Gene Expression and Trekker rows cannot "
                "refer to the same FASTQ files: " + ", ".join(sorted(shared_fastqs))
            )

        source_fastqs_by_sample: "OrderedDict[str, List[str]]" = OrderedDict()
        for row in gex_rows + trekker_rows:
            source_fastqs = source_fastqs_by_sample.setdefault(row["Sample"], [])
            for read in sorted(row["fastqs"]):
                source_fastqs.extend(row["fastqs"][read])
        source_fastqs_by_sample = OrderedDict(
            (sample, _unique(fastqs))
            for sample, fastqs in source_fastqs_by_sample.items()
        )

        record = {
            "sample": final_sample,
            "sc_sample": gex_rows[0]["Sample"],
            "experiment_date": "",
            "barcode_file": barcode_file,
            "fastq_1": str(concat_root / f"{final_sample}_R1.fastq.gz"),
            "fastq_2": str(concat_root / f"{final_sample}_R2.fastq.gz"),
            "sc_outdir": str(
                cellranger_root
                / final_sample
                / "outs"
                / "filtered_feature_bc_matrix"
            ),
            "sc_platform": SUPPORTED_PLATFORM,
            "profile": profile,
            "subsample": subsample,
            "cores": str(cores),
            "transcriptome": str(transcriptome_path),
            "cellranger_fastq_dirs": _unique(
                str(Path(fastq).parent)
                for row in gex_rows
                for read in ("R1", "R2")
                for fastq in row["fastqs"][read]
            ),
            "cellranger_fastq_samples": _unique(row["Sample"] for row in gex_rows),
            "cellranger_fastqs": gex_fastqs,
            "trekker_r1_fastqs": trekker_r1,
            "trekker_r2_fastqs": trekker_r2,
            "source_fastqs_by_sample": source_fastqs_by_sample,
        }
        for column in TREKKER_COLUMNS:
            value = str(record[column])
            if "," in value or "\n" in value or "\r" in value:
                raise TrekkerConfigError(
                    f"Sample '{final_sample}': generated vendor field '{column}' "
                    "cannot contain commas or newlines"
                )
        records[final_sample] = record
    return records


def cellranger_fastqs(record: Mapping[str, Any]) -> List[str]:
    """Return the GEX FASTQs tracked for a normalized libraries.csv record."""

    return list(record["cellranger_fastqs"])


def stage_fastqs(fastqs: Sequence[os.PathLike], destination: os.PathLike) -> Path:
    """Atomically symlink one FASTQ or concatenate several gzip streams."""

    sources = [Path(path).expanduser().resolve() for path in fastqs]
    if not sources:
        raise TrekkerConfigError("At least one FASTQ is required for staging")
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise TrekkerConfigError("FASTQ input does not exist: " + ", ".join(missing))

    # Do not resolve the destination: an earlier single-flowcell run may have
    # left a symlink here, and resolving it would target the source FASTQ.
    destination = Path(os.path.abspath(Path(destination).expanduser()))
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=str(destination.parent),
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        if len(sources) == 1:
            temporary.unlink()
            temporary.symlink_to(sources[0])
        else:
            with temporary.open("wb") as output_handle:
                for source in sources:
                    with source.open("rb") as input_handle:
                        shutil.copyfileobj(input_handle, output_handle, length=16 * 1024 * 1024)
        os.replace(temporary, destination)
    except Exception:
        if temporary.exists() or temporary.is_symlink():
            temporary.unlink()
        raise
    return destination


def _validate_row(
    row: MutableMapping[str, str],
    row_number: int,
) -> Dict[str, str]:
    row = {key: (value or "").strip() for key, value in row.items()}
    sample = row.get("sample", "")
    prefix = f"Row {row_number}" if not sample else f"Sample '{sample}'"

    for column in TREKKER_COLUMNS:
        if column != "experiment_date" and not row.get(column):
            raise TrekkerConfigError(f"{prefix}: column '{column}' is required")
        if "," in row[column] or "\n" in row[column] or "\r" in row[column]:
            raise TrekkerConfigError(
                f"{prefix}: column '{column}' cannot contain commas or newlines because "
                "the vendor launcher does not use a CSV-aware parser"
            )

    for field in ("sample", "sc_sample"):
        if not SAFE_SAMPLE_RE.fullmatch(row[field]):
            raise TrekkerConfigError(
                f"{prefix}: '{field}' must contain only letters, numbers, underscores, "
                "and hyphens"
            )

    if row["experiment_date"]:
        try:
            datetime.strptime(row["experiment_date"], "%Y%m%d")
        except ValueError as exc:
            raise TrekkerConfigError(
                f"{prefix}: when supplied, 'experiment_date' must be a valid date "
                "in YYYYMMDD format"
            ) from exc

    if row["sc_platform"] != SUPPORTED_PLATFORM:
        raise TrekkerConfigError(
            f"{prefix}: only sc_platform={SUPPORTED_PLATFORM} is currently supported; "
            f"received {row['sc_platform']!r}"
        )

    profile = row["profile"].lower()
    if profile not in SUPPORTED_PROFILES:
        raise TrekkerConfigError(
            f"{prefix}: 'profile' must be one of {sorted(SUPPORTED_PROFILES)}"
        )
    row["profile"] = profile

    subsample = row["subsample"].lower()
    if subsample not in {"yes", "no"}:
        raise TrekkerConfigError(f"{prefix}: 'subsample' must be 'yes' or 'no'")
    row["subsample"] = subsample

    try:
        cores = int(row["cores"])
    except ValueError as exc:
        raise TrekkerConfigError(f"{prefix}: 'cores' must be a positive integer") from exc
    if cores < 1:
        raise TrekkerConfigError(f"{prefix}: 'cores' must be a positive integer")
    row["cores"] = str(cores)

    row["barcode_file"] = _require_absolute_existing_path(
        row["barcode_file"], "barcode_file", sample
    )
    row["fastq_1"] = _require_absolute_existing_path(row["fastq_1"], "fastq_1", sample)
    row["fastq_2"] = _require_absolute_existing_path(row["fastq_2"], "fastq_2", sample)
    row["sc_outdir"] = _require_absolute_existing_path(
        row["sc_outdir"], "sc_outdir", sample, directory=True
    )

    fastq_1_name = Path(row["fastq_1"]).name
    fastq_2_name = Path(row["fastq_2"]).name
    if "R1" not in fastq_1_name or not fastq_1_name.endswith(".fastq.gz"):
        raise TrekkerConfigError(
            f"{prefix}: fastq_1 must contain 'R1' and end with '.fastq.gz'"
        )
    if "R2" not in fastq_2_name or not fastq_2_name.endswith(".fastq.gz"):
        raise TrekkerConfigError(
            f"{prefix}: fastq_2 must contain 'R2' and end with '.fastq.gz'"
        )

    for matrix_file in MATRIX_FILES:
        matrix_path = Path(row["sc_outdir"]) / matrix_file
        if not matrix_path.is_file():
            raise TrekkerConfigError(
                f"{prefix}: required quality-filtered 10x matrix file is missing: "
                f"{matrix_path}"
            )

    return dict(row)


def load_samplesheet(samplesheet: os.PathLike) -> "OrderedDict[str, Dict[str, str]]":
    """Load and validate a vendor-format TrekkerU_CX samplesheet."""

    path = Path(samplesheet).expanduser().resolve()
    if not path.is_file():
        raise TrekkerConfigError(f"Trekker samplesheet does not exist: {path}")

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        actual_columns = tuple(reader.fieldnames or ())
        if actual_columns != TREKKER_COLUMNS:
            raise TrekkerConfigError(
                "Trekker vendor samplesheet columns must appear in this exact order: "
                + ",".join(TREKKER_COLUMNS)
                + f". Received: {','.join(actual_columns)}"
            )

        records: "OrderedDict[str, Dict[str, str]]" = OrderedDict()
        output_keys = set()
        for row_number, row in enumerate(reader, start=2):
            record = _validate_row(row, row_number)
            sample = record["sample"]
            if sample in records:
                raise TrekkerConfigError(f"Duplicate Trekker sample name: {sample}")
            output_key = (record["experiment_date"], sample)
            if output_key in output_keys:
                raise TrekkerConfigError(
                    f"Duplicate Trekker output directory for sample '{sample}'"
                )
            records[sample] = record
            output_keys.add(output_key)

    if not records:
        raise TrekkerConfigError("Trekker samplesheet must contain at least one data row")
    return records


def write_vendor_samplesheet(record: Mapping[str, str], destination: os.PathLike) -> None:
    """Write exactly one row in the column order required by Takara's launcher."""

    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TREKKER_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerow({column: record[column] for column in TREKKER_COLUMNS})


def validate_installation(
    installation: os.PathLike,
    profiles: Iterable[str],
    conda_environment: Optional[os.PathLike] = None,
) -> Path:
    """Validate the vendor installation and selected execution profiles."""

    installation = Path(installation).expanduser().resolve()
    required = [
        installation / "nuclei_locater_toplevel.sh",
        installation / "common",
        installation / "cellbarcode_whitelists",
    ]
    profiles = set(profiles)
    for profile in profiles:
        required.append(installation / f"nuclei_locater_{profile}.sh")
    if "singularity" in profiles:
        required.append(installation / "trekker-v1.4.11.sif")

    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise TrekkerConfigError(
            "Trekker installation is incomplete; missing: " + ", ".join(missing)
        )

    if "conda" in profiles:
        if conda_environment is None:
            raise TrekkerConfigError("A Trekker conda environment is required for profile=conda")
        conda_environment = Path(conda_environment).expanduser().resolve()
        conda_required = [conda_environment / "bin/python", conda_environment / "bin/Rscript"]
        conda_missing = [str(path) for path in conda_required if not path.is_file()]
        if conda_missing:
            raise TrekkerConfigError(
                "Trekker conda environment is incomplete; missing: "
                + ", ".join(conda_missing)
            )
    return installation


def _replace_one(text: str, pattern: str, replacement: str, description: str) -> str:
    updated, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise TrekkerConfigError(
            f"Could not configure {description}; the vendor launcher format changed"
        )
    return updated


def _prepare_generated_file(destination: Path) -> None:
    """Ensure writing a generated file can never follow a stale symlink."""

    if destination.is_symlink():
        destination.unlink()
    elif destination.is_dir():
        raise TrekkerConfigError(
            f"Refusing to replace unexpected directory in generated runtime: {destination}"
        )


def prepare_runtime(
    installation: os.PathLike,
    output_root: os.PathLike,
    runtime_directory: os.PathLike,
    profile: str,
    conda_environment: Optional[os.PathLike] = None,
) -> Path:
    """Create a small runtime overlay without modifying Takara's installation."""

    profile = profile.lower()
    installation = validate_installation(
        installation, {profile}, conda_environment=conda_environment
    )
    output_root = Path(output_root).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    runtime = Path(runtime_directory).expanduser().resolve()
    runtime.mkdir(parents=True, exist_ok=True)

    top_source = installation / "nuclei_locater_toplevel.sh"
    top_destination = runtime / top_source.name
    _prepare_generated_file(top_destination)
    top_text = top_source.read_text(encoding="utf-8")
    top_text = _replace_one(
        top_text,
        r"^OUT_DIR=.*$",
        "OUT_DIR=" + shlex.quote(str(output_root)),
        "Trekker output directory",
    )
    top_text = _replace_one(
        top_text,
        r"^SCRIPT_DIR=.*$",
        "SCRIPT_DIR=" + shlex.quote(str(installation)),
        "Trekker installation directory",
    )
    if profile == "singularity":
        singularity_setup = """module load singularity
if ! command -v singularity >/dev/null 2>&1; then
   echo "Trekker requires the singularity command for profile=singularity" >&2
   exit 1
fi"""
        if re.search(r"^set -e\s*$", top_text, flags=re.MULTILINE):
            top_text = _replace_one(
                top_text,
                r"^set -e\s*$",
                r"\g<0>" + "\n\n" + singularity_setup,
                "Trekker Singularity module setup",
            )
        else:
            top_text = _replace_one(
                top_text,
                r"^#!/usr/bin/env bash\s*$",
                r"\g<0>" + "\n\n" + singularity_setup,
                "Trekker Singularity module setup",
            )

    conda_source = installation / "nuclei_locater_conda.sh"
    conda_destination = runtime / conda_source.name
    if profile == "conda":
        top_text = _replace_one(
            top_text,
            r'^\s*source "\$SCRIPT_DIR/nuclei_locater_conda\.sh"\s*$',
            "      source " + shlex.quote(str(conda_destination)),
            "Trekker conda launcher",
        )
    top_destination.write_text(top_text, encoding="utf-8")
    top_destination.chmod(0o755)

    if profile == "conda":
        _prepare_generated_file(conda_destination)
        conda_path = Path(conda_environment).expanduser().resolve()
        conda_text = conda_source.read_text(encoding="utf-8")
        conda_text = _replace_one(
            conda_text,
            r"^PROFILE_CONDA_PATH=.*$",
            "PROFILE_CONDA_PATH=" + shlex.quote(str(conda_path)),
            "Trekker conda path",
        )
        conda_text = _replace_one(
            conda_text,
            r"^source .*conda\.sh\s*$",
            ": # SF_scMaestro activates the supplied environment below",
            "Trekker conda initialization",
        )
        activation = """export CONDA_PREFIX=\"${PROFILE_CONDA_PATH}\"
export CONDA_DEFAULT_ENV=\"${PROFILE_CONDA_PATH}\"
export CONDA_SHLVL=1
export PATH=\"${PROFILE_CONDA_PATH}/bin:${PATH}\"
export R_HOME=\"${PROFILE_CONDA_PATH}/lib/R\"
export R_LIBS_SITE=\"${PROFILE_CONDA_PATH}/lib/R/library\"
for activation_script in \"${PROFILE_CONDA_PATH}\"/etc/conda/activate.d/*.sh; do
   [ -f \"${activation_script}\" ] && source \"${activation_script}\"
done
unset activation_script"""
        conda_text = _replace_one(
            conda_text,
            r"^conda activate .*?$",
            activation,
            "Trekker conda activation",
        )
        conda_destination.write_text(conda_text, encoding="utf-8")
        conda_destination.chmod(0o755)

    return top_destination


def expected_report(output_root: os.PathLike, record: Mapping[str, str]) -> Path:
    return (
        Path(output_root).expanduser().resolve()
        / f"{record['experiment_date']}_{record['sample']}"
        / f"trekker_{record['sample']}"
        / "output"
        / f"{record['sample']}_Trekker_Report.html"
    )


def run_trekker(
    samplesheet: os.PathLike,
    installation: os.PathLike,
    output_root: os.PathLike,
    runtime_directory: os.PathLike,
    conda_environment: Optional[os.PathLike] = None,
    stdout=None,
) -> Path:
    """Validate and execute a one-row vendor samplesheet."""

    records = load_samplesheet(samplesheet)
    if len(records) != 1:
        raise TrekkerConfigError(
            "The generated vendor samplesheet must contain exactly one data row"
        )
    record = next(iter(records.values()))
    launcher = prepare_runtime(
        installation=installation,
        output_root=output_root,
        runtime_directory=runtime_directory,
        profile=record["profile"],
        conda_environment=conda_environment,
    )
    subprocess.run(
        ["bash", str(launcher), str(Path(samplesheet).resolve())],
        check=True,
        stdout=stdout,
        stderr=subprocess.STDOUT if stdout is not None else None,
    )
    report = expected_report(output_root, record)
    if not report.is_file():
        raise RuntimeError(
            f"Trekker exited successfully but the expected report was not created: {report}"
        )
    return report


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samplesheet")
    parser.add_argument("--installation", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--runtime-directory", required=True)
    parser.add_argument("--conda-environment")
    args = parser.parse_args(argv)
    run_trekker(
        samplesheet=args.samplesheet,
        installation=args.installation,
        output_root=args.output_root,
        runtime_directory=args.runtime_directory,
        conda_environment=args.conda_environment,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
