from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Sequence

from .defaults import (
    DEFAULT_DELIMITER_NAME,
    DEFAULT_MASK_TEXT,
    DEFAULT_MAX_ITERATIONS,
    DEFAULT_OUTPUT_PRIMER_LIST_PATH,
    DEFAULT_PRIMER_LIST_NAME,
    DEFAULT_PRIMER_LIST_PATH,
    DEFAULT_PRIMER_NAME_PREFIX,
    DEFAULT_TARGET_COVERAGE,
    DELIMITER_OPTIONS,
    MASK_TEXT_PLACEHOLDER,
    TARGET_COVERAGE_OPTIONS,
)
from .core import (
    DesignResult,
    PrimerInput,
    delimiter_from_name,
    design_primers,
    format_interval,
    merge_primer_inputs,
    parse_masks,
    parse_primer_list,
    read_sequence,
    write_primer_list,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sanger-designer",
        description="Design Sanger sequencing primers for circular plasmid coverage.",
    )
    parser.add_argument("sequence", help="Input FASTA or GenBank file.")
    parser.add_argument(
        "-p",
        "--primers",
        action="append",
        default=[],
        metavar="PATH",
        help=(
            "Existing primer list txt file. Can be repeated. "
            f"Default: {DEFAULT_PRIMER_LIST_NAME} when no --primers is supplied."
        ),
    )
    parser.add_argument(
        "--include-default-primers",
        action="store_true",
        help=f"Include {DEFAULT_PRIMER_LIST_NAME} before files supplied with --primers.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=DEFAULT_OUTPUT_PRIMER_LIST_PATH,
        help=f"Output primer list txt file. Default: {DEFAULT_OUTPUT_PRIMER_LIST_PATH}",
    )
    parser.add_argument(
        "-c",
        "--coverage",
        type=int,
        choices=TARGET_COVERAGE_OPTIONS,
        default=DEFAULT_TARGET_COVERAGE,
        help=f"Target coverage. Default: {DEFAULT_TARGET_COVERAGE}",
    )
    parser.add_argument(
        "-m",
        "--mask",
        default=DEFAULT_MASK_TEXT,
        help=f"Comma-separated mask intervals, e.g. {MASK_TEXT_PLACEHOLDER}",
    )
    parser.add_argument(
        "--delimiter",
        choices=DELIMITER_OPTIONS,
        default=DEFAULT_DELIMITER_NAME,
        help=f"Output delimiter. Default: {DEFAULT_DELIMITER_NAME}",
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=DEFAULT_MAX_ITERATIONS,
        help=f"Maximum new-design iterations. Default: {DEFAULT_MAX_ITERATIONS}",
    )
    parser.add_argument(
        "--prefix",
        default=DEFAULT_PRIMER_NAME_PREFIX,
        help=f"Name prefix for newly designed primers. Default: {DEFAULT_PRIMER_NAME_PREFIX}",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        sequence = read_sequence(args.sequence)
        existing = read_existing_primers(
            args.primers,
            include_default=args.include_default_primers,
        )
        masks = parse_masks(args.mask, len(sequence))
        result = design_primers(
            sequence,
            existing,
            target_coverage=args.coverage,
            masks=masks,
            max_iterations=args.max_iterations,
            primer_name_prefix=args.prefix,
        )
        write_primer_list(args.output, result.primers, delimiter_from_name(args.delimiter))
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print_summary(result, args.output)
    return 0 if result.primers else 1


def read_existing_primers(
    paths: Sequence[str | Path] | str | Path,
    *,
    include_default: bool = False,
) -> list[PrimerInput]:
    supplied_paths = normalize_primer_list_paths(paths)
    primer_paths = [DEFAULT_PRIMER_LIST_PATH, *supplied_paths] if include_default or not supplied_paths else supplied_paths
    primer_groups = []
    for primers_path in primer_paths:
        if not primers_path.exists():
            print(
                f"warning: primer list not found: {primers_path}; continuing without this list.",
                file=sys.stderr,
            )
            continue
        primer_groups.append(parse_primer_list(primers_path))
    return merge_primer_inputs(primer_groups)


def normalize_primer_list_paths(paths: Sequence[str | Path] | str | Path) -> list[Path]:
    if isinstance(paths, (str, Path)):
        return [Path(paths)]
    return [Path(path) for path in paths]


def print_summary(result: DesignResult, output: str) -> None:
    existing_count = sum(1 for primer in result.primers if primer.source == "existing")
    designed_count = sum(1 for primer in result.primers if primer.source == "designed")
    print(f"Output: {output}")
    print(f"Primers: {len(result.primers)} total ({existing_count} existing, {designed_count} designed)")
    print(f"Minimum coverage: {result.min_coverage} / target {result.target_coverage}")
    if result.achieved:
        print("Coverage target achieved.")
    else:
        regions = ", ".join(format_interval(region) for region in result.missing_regions)
        print(f"Coverage target not achieved. Missing regions: {regions}")


if __name__ == "__main__":
    raise SystemExit(main())
