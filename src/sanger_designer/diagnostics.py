from __future__ import annotations

import csv
from dataclasses import asdict
from datetime import datetime
import io
import json
from pathlib import Path
import re
import sys
import zipfile
from typing import Any, Mapping, Sequence

from . import defaults
from .core import DesignResult, DesignSettings, Interval, Primer, PrimerInput, format_interval, format_primer_list


SENSITIVE_BUNDLE_NOTICE = "Bundleにはプラスミド配列やプライマー配列が入ります。共有前に内容を確認してください。"


def build_support_bundle(
    *,
    result: DesignResult,
    sequence_length: int,
    settings: DesignSettings,
    run_parameters: Mapping[str, Any],
    sequence_file_name: str,
    sequence_text: str,
    uploaded_primer_files: Sequence[Mapping[str, str]],
    default_primer_text: str | None,
    merged_primers: Sequence[PrimerInput],
    masks: Sequence[Interval],
    warnings: Sequence[str],
    runtime_seconds: float,
    delimiter: str = "\t",
    created_at: datetime | None = None,
) -> bytes:
    created = created_at or datetime.now().astimezone()
    manifest = build_manifest(result, sequence_length, runtime_seconds, created)
    settings_payload = build_settings_payload(settings, run_parameters, uploaded_primer_files)
    defaults_payload = exported_defaults()
    run_log = build_run_log(
        result=result,
        sequence_length=sequence_length,
        uploaded_primer_files=uploaded_primer_files,
        default_primer_text=default_primer_text,
        merged_primers=merged_primers,
        masks=masks,
        warnings=warnings,
        runtime_seconds=runtime_seconds,
    )

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        write_json(bundle, "settings/manifest.json", manifest)
        write_json(bundle, "settings/settings.json", settings_payload)
        write_json(bundle, "settings/defaults.json", defaults_payload)

        bundle.writestr(f"inputs/{safe_filename(sequence_file_name, 'plasmid_input.txt')}", sequence_text)
        if default_primer_text is not None:
            bundle.writestr(f"settings/{defaults.DEFAULT_PRIMER_LIST_NAME}", default_primer_text)
        for index, primer_file in enumerate(uploaded_primer_files, start=1):
            name = safe_filename(primer_file["name"], f"primer_list_{index}.txt")
            bundle.writestr(f"inputs/primer_lists/{index:02d}_{name}", primer_file["text"])

        bundle.writestr("intermediate/merged_primers.txt", format_primer_list_from_inputs(merged_primers, delimiter))
        bundle.writestr("intermediate/masks.txt", "\n".join(format_interval(mask) for mask in masks) + ("\n" if masks else ""))

        bundle.writestr("results/sanger_primers.txt", format_primer_list(result.primers, delimiter))
        bundle.writestr("results/depth_report.txt", depth_report_text(result, sequence_length, settings))
        bundle.writestr("results/primer_table.csv", primer_table_csv(result.primers))
        bundle.writestr("results/depth_by_position.csv", depth_by_position_csv(result.depth))

        bundle.writestr("logs/run_log.txt", run_log)
        bundle.writestr("logs/warnings.txt", "\n".join(warnings) + ("\n" if warnings else ""))

    return buffer.getvalue()


def default_bundle_filename(created_at: datetime | None = None) -> str:
    created = created_at or datetime.now()
    return f"sanger_designer_bundle_{created.strftime('%Y%m%d_%H%M%S')}.zip"


def build_manifest(
    result: DesignResult,
    sequence_length: int,
    runtime_seconds: float,
    created_at: datetime,
) -> dict[str, Any]:
    return {
        "tool": "sanger-designer",
        "created_at": created_at.isoformat(),
        "mode": "GUI",
        "python_version": sys.version,
        "sequence_length": sequence_length,
        "target_depth": result.target_depth,
        "achieved": result.achieved,
        "min_depth": result.min_depth,
        "total_primers": len(result.primers),
        "existing_primers": sum(1 for primer in result.primers if primer.source == "existing"),
        "designed_primers": sum(1 for primer in result.primers if primer.source == "designed"),
        "runtime_seconds": runtime_seconds,
    }


def build_settings_payload(
    settings: DesignSettings,
    run_parameters: Mapping[str, Any],
    uploaded_primer_files: Sequence[Mapping[str, str]],
) -> dict[str, Any]:
    return {
        "run_parameters": {
            "target_depth": run_parameters.get("depth"),
            "mask_text": run_parameters.get("mask_text"),
            "primer_name_prefix": run_parameters.get("prefix"),
            "delimiter": run_parameters.get("delimiter"),
            "max_iterations": run_parameters.get("max_iterations"),
            "use_default_primers": run_parameters.get("use_default_primers"),
            "uploaded_primer_files": [primer_file["name"] for primer_file in uploaded_primer_files],
        },
        "design_settings": asdict(settings),
    }


def exported_defaults() -> dict[str, Any]:
    payload: dict[str, Any] = {}
    for name in dir(defaults):
        if not name.isupper():
            continue
        value = getattr(defaults, name)
        if isinstance(value, Path):
            payload[name] = str(value)
        elif isinstance(value, tuple):
            payload[name] = list(value)
        elif isinstance(value, (str, int, float, bool)) or value is None:
            payload[name] = value
    return payload


def build_run_log(
    *,
    result: DesignResult,
    sequence_length: int,
    uploaded_primer_files: Sequence[Mapping[str, str]],
    default_primer_text: str | None,
    merged_primers: Sequence[PrimerInput],
    masks: Sequence[Interval],
    warnings: Sequence[str],
    runtime_seconds: float,
) -> str:
    total_input_primers = count_primer_rows(default_primer_text or "") + sum(
        count_primer_rows(primer_file["text"]) for primer_file in uploaded_primer_files
    )
    duplicate_count = max(0, total_input_primers - len(merged_primers))
    lines = [
        "Sanger Designer support bundle",
        f"Sequence length: {sequence_length}",
        f"Default primer list included: {'yes' if default_primer_text is not None else 'no'}",
        f"Uploaded primer lists: {len(uploaded_primer_files)}",
        f"Input primer rows: {total_input_primers}",
        f"Merged primer rows: {len(merged_primers)}",
        f"Duplicate primer sequences skipped: {duplicate_count}",
        f"Masks: {', '.join(format_interval(mask) for mask in masks) if masks else 'none'}",
        f"Achieved: {result.achieved}",
        f"Minimum depth: {result.min_depth} / {result.target_depth}",
        f"Total output primers: {len(result.primers)}",
        f"Existing output primers: {sum(1 for primer in result.primers if primer.source == 'existing')}",
        f"Designed output primers: {sum(1 for primer in result.primers if primer.source == 'designed')}",
        f"Missing regions: {', '.join(format_interval(region) for region in result.missing_regions) if result.missing_regions else 'none'}",
        f"Runtime seconds: {runtime_seconds:.3f}",
    ]
    if warnings:
        lines.append("Warnings:")
        lines.extend(f"  {warning}" for warning in warnings)
    return "\n".join(lines) + "\n"


def count_primer_rows(text: str) -> int:
    return sum(1 for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#"))


def format_primer_list_from_inputs(primers: Sequence[PrimerInput], delimiter: str) -> str:
    rows = [delimiter.join([primer.name, primer.sequence, primer.memo]) for primer in primers]
    return "\n".join(rows) + ("\n" if rows else "")


def primer_table_csv(primers: Sequence[Primer]) -> str:
    rows = [
        {
            "name": primer.name,
            "sequence": primer.sequence,
            "direction": primer.direction,
            "position": primer.position,
            "binding": format_interval(primer.binding),
            "cover": format_interval(primer.cover),
            "source": primer.source,
        }
        for primer in sorted(primers, key=lambda item: item.position)
    ]
    return rows_to_csv(rows, ["name", "sequence", "direction", "position", "binding", "cover", "source"])


def depth_by_position_csv(depth: Sequence[int]) -> str:
    rows = [{"position": position, "depth": value} for position, value in enumerate(depth, start=1)]
    return rows_to_csv(rows, ["position", "depth"])


def rows_to_csv(rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def depth_report_text(
    result: DesignResult,
    sequence_length: int,
    settings: DesignSettings,
) -> str:
    lines = [
        f"Status: {'achieved' if result.achieved else 'incomplete'}",
        f"Target depth: {result.target_depth}",
        f"Minimum depth: {result.min_depth}",
        f"Total primers: {len(result.primers)}",
        f"Existing primers: {sum(1 for primer in result.primers if primer.source == 'existing')}",
        f"Designed primers: {sum(1 for primer in result.primers if primer.source == 'designed')}",
        f"Sequence length: {sequence_length}",
        "Design settings:",
        f"  Read length: {settings.read_length}",
        f"  Noisy bases: {settings.noise_length}",
        f"  Primer length: {settings.min_primer_len}-{settings.max_primer_len}",
        f"  Tm: {settings.min_tm}-{settings.max_tm}",
        f"  GC: {settings.min_gc}-{settings.max_gc}",
        f"  Minimum binding gap: {settings.min_binding_gap}",
        f"  Preferred pair distance: {settings.preferred_pair_min}-{settings.preferred_pair_max}",
    ]
    if result.missing_regions:
        lines.append("Missing regions: " + ", ".join(format_interval(region) for region in result.missing_regions))
    return "\n".join(lines) + "\n"


def write_json(bundle: zipfile.ZipFile, path: str, payload: Mapping[str, Any]) -> None:
    bundle.writestr(path, json.dumps(payload, ensure_ascii=False, indent=2) + "\n")


def safe_filename(name: str, fallback: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    return cleaned or fallback

