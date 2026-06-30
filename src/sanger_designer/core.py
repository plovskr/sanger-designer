from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import re
from statistics import median
from typing import Iterable, Sequence

from .defaults import (
    DEFAULT_MAX_ITERATIONS,
    DEFAULT_PRIMER_NAME_PREFIX,
    DEFAULT_TARGET_DEPTH,
    DESIGN_BEAM_WIDTH,
    DESIGN_DIVERSE_POSITION_GAP,
    DESIGN_OPTIONS_PER_REGION,
    DESIGN_OPTIONS_PER_STATE,
    DESIGN_REGIONS_PER_STATE,
    EFFECTIVE_READ_LENGTH,
    MAX_GC,
    MAX_PRIMER_LEN,
    MAX_TM,
    MIN_BINDING_GAP,
    MIN_GC,
    MIN_PRIMER_LEN,
    MIN_TM,
    NOISE_LENGTH,
    PREFERRED_PAIR_MAX,
    PREFERRED_PAIR_MIN,
    READ_LENGTH,
    TARGET_DEPTH_OPTIONS,
)

DNA_RE = re.compile(r"^[ACGT]+$", re.IGNORECASE)


@dataclass(frozen=True)
class DesignSettings:
    read_length: int = READ_LENGTH
    noise_length: int = NOISE_LENGTH
    min_primer_len: int = MIN_PRIMER_LEN
    max_primer_len: int = MAX_PRIMER_LEN
    min_tm: float = MIN_TM
    max_tm: float = MAX_TM
    min_gc: float = MIN_GC
    max_gc: float = MAX_GC
    min_binding_gap: int = MIN_BINDING_GAP
    preferred_pair_min: int = PREFERRED_PAIR_MIN
    preferred_pair_max: int = PREFERRED_PAIR_MAX

    @property
    def effective_read_length(self) -> int:
        return self.read_length - self.noise_length

    def validate(self) -> None:
        if self.read_length < 1:
            raise ValueError("read_length must be at least 1")
        if self.noise_length < 0:
            raise ValueError("noise_length must be non-negative")
        if self.noise_length >= self.read_length:
            raise ValueError("noise_length must be smaller than read_length")
        if self.min_primer_len < 1 or self.max_primer_len < self.min_primer_len:
            raise ValueError("primer length range is invalid")
        if self.max_tm < self.min_tm:
            raise ValueError("Tm range is invalid")
        if self.max_gc < self.min_gc:
            raise ValueError("GC range is invalid")
        if self.min_binding_gap < 0:
            raise ValueError("min_binding_gap must be non-negative")
        if self.preferred_pair_max < self.preferred_pair_min:
            raise ValueError("preferred pair range is invalid")


DEFAULT_DESIGN_SETTINGS = DesignSettings()


@dataclass(frozen=True)
class Interval:
    start: int
    end: int


@dataclass(frozen=True)
class PrimerInput:
    name: str
    sequence: str
    memo: str = ""
    order: int = 0


@dataclass(frozen=True)
class Primer:
    name: str
    sequence: str
    direction: str
    position: int
    binding: Interval
    cover: Interval
    source: str
    input_memo: str = ""
    order: int = 0

    @property
    def output_memo(self) -> str:
        parts = [
            f"Position: {self.position}",
            f"Direction: {self.direction}",
            f"Binding: {format_interval(self.binding)}",
            f"Cover: {format_interval(self.cover)}",
        ]
        if self.input_memo:
            parts.append(self.input_memo)
        return ", ".join(parts)


@dataclass(frozen=True)
class DesignResult:
    primers: tuple[Primer, ...]
    depth: tuple[int, ...]
    target_depth: int
    missing_regions: tuple[Interval, ...]
    achieved: bool

    @property
    def min_depth(self) -> int:
        return min(self.depth) if self.depth else 0


class Primer3UnavailableError(RuntimeError):
    pass


def design_primers(
    sequence: str,
    existing_primers: Sequence[PrimerInput],
    *,
    target_depth: int = DEFAULT_TARGET_DEPTH,
    masks: Sequence[Interval] = (),
    max_iterations: int = DEFAULT_MAX_ITERATIONS,
    primer_name_prefix: str = DEFAULT_PRIMER_NAME_PREFIX,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> DesignResult:
    if target_depth not in TARGET_DEPTH_OPTIONS:
        allowed = ", ".join(str(value) for value in TARGET_DEPTH_OPTIONS)
        raise ValueError(f"target_depth must be one of {allowed}")
    settings.validate()

    seq = normalize_sequence(sequence)
    n = len(seq)
    if n == 0:
        raise ValueError("input sequence is empty")

    normalized_masks = tuple(normalize_interval(mask, n) for mask in masks)
    name_prefix = normalize_primer_name_prefix(primer_name_prefix)
    selected = select_existing_primers(seq, existing_primers, target_depth, normalized_masks, settings)
    depth = compute_depth(n, selected)

    if min(depth) < target_depth:
        selected = design_missing_primers_beam(
            seq,
            selected,
            target_depth,
            normalized_masks,
            max_iterations=max_iterations,
            primer_name_prefix=name_prefix,
            settings=settings,
        )

    selected = optimize_primers(seq, selected, target_depth)
    depth = compute_depth(n, selected)
    missing = depth_regions_below(depth, target_depth)
    return DesignResult(
        primers=tuple(selected),
        depth=tuple(depth),
        target_depth=target_depth,
        missing_regions=tuple(missing),
        achieved=not missing,
    )


def normalize_primer_name_prefix(prefix: str) -> str:
    normalized = prefix.strip() or DEFAULT_PRIMER_NAME_PREFIX
    if any(separator in normalized for separator in ("\t", "\n", "\r", ",", ";")):
        raise ValueError("primer name prefix must not contain tab, newline, comma, or semicolon")
    return normalized


def design_missing_primers_beam(
    sequence: str,
    initial_primers: Sequence[Primer],
    target_depth: int,
    masks: Sequence[Interval],
    *,
    max_iterations: int,
    primer_name_prefix: str,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
    beam_width: int = DESIGN_BEAM_WIDTH,
) -> list[Primer]:
    sequence_length = len(sequence)
    initial = tuple(initial_primers)
    if min(compute_depth(sequence_length, initial)) >= target_depth:
        return list(initial)

    states: dict[tuple[tuple[str, str, int], ...], tuple[Primer, ...]] = {state_key(initial): initial}
    best_state = initial
    best_score = design_state_score(sequence_length, initial, target_depth)

    for _ in range(max_iterations):
        candidates: dict[tuple[tuple[str, str, int], ...], tuple[Primer, ...]] = {}
        for state in states.values():
            depth = compute_depth(sequence_length, state)
            for option in generate_design_options(
                sequence,
                state,
                depth,
                target_depth,
                masks,
                primer_name_prefix,
                settings,
            ):
                next_state = (*state, *option)
                key = state_key(next_state)
                if key in candidates:
                    continue
                candidates[key] = next_state

        if not candidates:
            break

        ranked = sorted(
            candidates.values(),
            key=lambda state: design_state_score(sequence_length, state, target_depth),
        )
        states = {state_key(state): state for state in ranked[:beam_width]}

        current_best = ranked[0]
        current_score = design_state_score(sequence_length, current_best, target_depth)
        if current_score < best_score:
            best_state = current_best
            best_score = current_score
        if min(compute_depth(sequence_length, current_best)) >= target_depth:
            return list(current_best)

    return list(best_state)


def state_key(primers: Sequence[Primer]) -> tuple[tuple[str, str, int], ...]:
    return tuple(sorted((primer.source, primer.direction, primer.position) for primer in primers))


def design_state_score(
    sequence_length: int,
    primers: Sequence[Primer],
    target_depth: int,
) -> tuple[int | float, ...]:
    depth = compute_depth(sequence_length, primers)
    missing_bases = sum(1 for value in depth if value < target_depth)
    depth_deficit = sum(max(0, target_depth - value) for value in depth)
    designed_count = sum(1 for primer in primers if primer.source == "designed")
    placement = placement_score(sequence_length, primers)
    overage_squared = sum(max(0, value - target_depth) ** 2 for value in depth)
    high_depth_bases = sum(1 for value in depth if value >= target_depth + 2)
    return (
        missing_bases,
        depth_deficit,
        designed_count,
        placement["same_direction_adjacent"],
        overage_squared,
        high_depth_bases,
        placement["max_adjacent_gap"],
        placement["fr_gap_variance"],
        len(primers),
    )


def generate_design_options(
    sequence: str,
    selected: Sequence[Primer],
    depth: Sequence[int],
    target_depth: int,
    masks: Sequence[Interval],
    primer_name_prefix: str,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> list[tuple[Primer, ...]]:
    sequence_length = len(sequence)
    anchors = positions_below_depth(depth, target_depth)
    if not anchors:
        return []

    options: list[tuple[Primer, ...]] = []
    seen: set[tuple[tuple[str, int], ...]] = set()
    min_depth = min(depth)
    regions = positions_to_intervals(anchors, sequence_length)[:DESIGN_REGIONS_PER_STATE]
    for region in regions:
        region_options = design_options_for_anchor(
            sequence,
            region.start,
            selected,
            masks,
            primer_name_prefix,
            settings,
            prefer_single=min_depth > 0,
        )
        for option in region_options:
            key = tuple((primer.direction, primer.position) for primer in option)
            if key in seen:
                continue
            seen.add(key)
            options.append(option)
            if len(options) >= DESIGN_OPTIONS_PER_STATE:
                return options
    return options


def design_options_for_anchor(
    sequence: str,
    anchor: int,
    selected: Sequence[Primer],
    masks: Sequence[Interval],
    primer_name_prefix: str,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
    *,
    prefer_single: bool,
) -> list[tuple[Primer, ...]]:
    options: list[tuple[Primer, ...]] = []

    if prefer_single:
        options.extend(
            (primer,)
            for primer in design_single_options_covering_anchor(
                sequence,
                anchor,
                selected,
                masks,
                name_prefix=primer_name_prefix,
                settings=settings,
            )
        )

    options.extend(
        design_pair_options_from_anchor(
            sequence,
            anchor,
            selected,
            masks,
            name_prefix=primer_name_prefix,
            settings=settings,
        )
    )

    if not prefer_single:
        options.extend(
            (primer,)
            for primer in design_single_options_covering_anchor(
                sequence,
                anchor,
                selected,
                masks,
                name_prefix=primer_name_prefix,
                settings=settings,
            )
        )

    return options[:DESIGN_OPTIONS_PER_REGION]


def normalize_sequence(sequence: str) -> str:
    seq = re.sub(r"[\s\d]", "", sequence.upper())
    if not seq:
        return ""
    if not DNA_RE.match(seq):
        raise ValueError("sequence contains unsupported bases")
    return seq


def read_sequence(path: str | Path) -> str:
    text = Path(path).read_text(encoding="utf-8")
    return read_sequence_text(text)


def read_sequence_text(text: str) -> str:
    stripped = text.lstrip()
    if stripped.startswith(">"):
        return parse_fasta(text)
    if re.search(r"^ORIGIN\b", text, flags=re.MULTILINE):
        return parse_genbank(text)
    return normalize_sequence(text)


def parse_fasta(text: str) -> str:
    header_count = sum(1 for line in text.splitlines() if line.startswith(">"))
    if header_count > 1:
        raise ValueError("FASTA input contains multiple records; please provide a single plasmid sequence")
    lines = [line.strip() for line in text.splitlines() if line.strip() and not line.startswith(">")]
    return normalize_sequence("".join(lines))


def parse_genbank(text: str) -> str:
    origin_count = len(re.findall(r"^ORIGIN\b", text, flags=re.MULTILINE))
    if origin_count > 1:
        raise ValueError("GenBank input contains multiple ORIGIN records; please provide a single plasmid sequence")
    if origin_count == 0:
        raise ValueError("GenBank input contains no ORIGIN sequence")

    in_origin = False
    chunks: list[str] = []
    for line in text.splitlines():
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if in_origin and line.startswith("//"):
            break
        if in_origin:
            chunks.append("".join(re.findall(r"[A-Za-z]+", line)))
    sequence = normalize_sequence("".join(chunks))
    if not sequence:
        raise ValueError("GenBank input contains no sequence in ORIGIN")
    return sequence


def parse_primer_list(path: str | Path) -> list[PrimerInput]:
    return parse_primer_list_text(Path(path).read_text(encoding="utf-8"))


def parse_primer_list_text(text: str) -> list[PrimerInput]:
    primers: list[PrimerInput] = []
    for index, raw in enumerate(text.splitlines()):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        delimiter = "\t" if "\t" in line else ";" if ";" in line else ","
        fields = [field.strip() for field in line.split(delimiter)]
        if len(fields) < 2:
            raise ValueError(f"invalid primer row: {raw}")
        name, sequence = fields[0], normalize_sequence(fields[1])
        memo = fields[2] if len(fields) > 2 else ""
        primers.append(PrimerInput(name=name, sequence=sequence, memo=memo, order=index))
    return primers


def merge_primer_inputs(primer_groups: Sequence[Sequence[PrimerInput]]) -> list[PrimerInput]:
    merged: list[PrimerInput] = []
    seen_sequences: set[str] = set()
    for primer_group in primer_groups:
        for primer in primer_group:
            if primer.sequence in seen_sequences:
                continue
            seen_sequences.add(primer.sequence)
            merged.append(
                PrimerInput(
                    name=primer.name,
                    sequence=primer.sequence,
                    memo=primer.memo,
                    order=len(merged),
                )
            )
    return merged


def parse_masks(mask_text: str, sequence_length: int) -> tuple[Interval, ...]:
    if not mask_text.strip():
        return ()
    masks = []
    for item in mask_text.split(","):
        match = re.fullmatch(r"\s*(\d+)\.\.(\d+)\s*", item)
        if not match:
            raise ValueError(f"invalid mask interval: {item}")
        masks.append(normalize_interval(Interval(int(match.group(1)), int(match.group(2))), sequence_length))
    return tuple(masks)


def find_unique_binding(
    sequence: str,
    primer: PrimerInput,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> Primer | None:
    hits = find_bindings(sequence, primer.sequence)
    if len(hits) != 1:
        return None
    return hit_to_primer(sequence, primer, hits[0], "existing", settings)


def find_bindings(sequence: str, primer_sequence: str) -> list[tuple[str, Interval, int]]:
    seq = normalize_sequence(sequence)
    primer = normalize_sequence(primer_sequence)
    rev = reverse_complement(primer)
    n = len(seq)
    hits: list[tuple[str, Interval, int]] = []
    seen: set[tuple[str, int, int]] = set()
    for start, end in circular_pattern_hits(seq, primer):
        key = ("Forward", start, end)
        if key not in seen:
            seen.add(key)
            hits.append(("Forward", Interval(start, end), start))
    for start, end in circular_pattern_hits(seq, rev):
        key = ("Reverse", start, end)
        if key not in seen:
            seen.add(key)
            hits.append(("Reverse", Interval(start, end), end))
    return hits


def circular_pattern_hits(sequence: str, pattern: str) -> list[tuple[int, int]]:
    n = len(sequence)
    if not pattern or len(pattern) > n:
        return []
    haystack = sequence + sequence[: len(pattern) - 1]
    hits: list[tuple[int, int]] = []
    index = haystack.find(pattern)
    while index != -1 and index < n:
        start = index + 1
        end = wrap(start + len(pattern) - 1, n)
        hits.append((start, end))
        index = haystack.find(pattern, index + 1)
    return hits


def hit_to_primer(
    sequence: str,
    primer: PrimerInput,
    hit: tuple[str, Interval, int],
    source: str,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> Primer:
    direction, binding, position = hit
    cover = depth_interval(position, direction, len(sequence), settings)
    return Primer(
        name=primer.name,
        sequence=primer.sequence,
        direction=direction,
        position=position,
        binding=binding,
        cover=cover,
        source=source,
        input_memo=primer.memo,
        order=primer.order,
    )


def select_existing_primers(
    sequence: str,
    primer_inputs: Sequence[PrimerInput],
    target_depth: int,
    masks: Sequence[Interval],
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> list[Primer]:
    del target_depth
    selected: list[Primer] = []
    seen_bindings: set[tuple[str, int, int]] = set()
    for primer_input in primer_inputs:
        primer = find_unique_binding(sequence, primer_input, settings)
        if primer and not interval_overlaps_any(primer.binding, masks, len(sequence)):
            key = (primer.direction, primer.binding.start, primer.binding.end)
            if key in seen_bindings:
                continue
            seen_bindings.add(key)
            selected.append(primer)
    return selected


def design_next_primers(
    sequence: str,
    selected: Sequence[Primer],
    depth: Sequence[int],
    target_depth: int,
    masks: Sequence[Interval],
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> tuple[Primer, ...]:
    n = len(sequence)
    anchors = positions_below_depth(depth, target_depth)
    if not anchors:
        return ()
    regions = positions_to_intervals(anchors, n)
    if min(depth) > 0:
        for region in regions:
            single = design_single_covering_anchor(sequence, region.start, selected, masks, settings)
            if single:
                return (single,)

    search_span = min(n, settings.preferred_pair_max + 120)
    for region in regions:
        anchor = region.start
        for offset in range(search_span):
            start = wrap(anchor + offset, n)
            forward = design_candidate(sequence, start, "Forward", selected, masks, settings=settings)
            if not forward:
                continue
            reverse = design_reverse_partner(sequence, forward, [*selected, forward], masks, settings)
            if reverse:
                return forward, reverse
    for region in regions:
        single = design_single_covering_anchor(sequence, region.start, selected, masks, settings)
        if single:
            return (single,)
    return ()


def design_single_covering_anchor(
    sequence: str,
    anchor: int,
    selected: Sequence[Primer],
    masks: Sequence[Interval],
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> Primer | None:
    options = design_single_options_covering_anchor(sequence, anchor, selected, masks, limit=1, settings=settings)
    return options[0] if options else None


def design_single_options_covering_anchor(
    sequence: str,
    anchor: int,
    selected: Sequence[Primer],
    masks: Sequence[Interval],
    *,
    limit: int = 4,
    name_prefix: str = DEFAULT_PRIMER_NAME_PREFIX,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> list[Primer]:
    n = len(sequence)
    options: list[Primer] = []
    seen: set[tuple[str, int]] = set()
    for delta in range(settings.noise_length, settings.read_length):
        forward_position = wrap(anchor - delta, n)
        forward = design_candidate(
            sequence,
            forward_position,
            "Forward",
            selected,
            masks,
            name_prefix=name_prefix,
            settings=settings,
        )
        if forward and anchor in interval_positions(forward.cover, n):
            key = (forward.direction, forward.position)
            if key not in seen:
                seen.add(key)
                options.append(forward)
                if len(options) >= limit:
                    return options

        reverse_position = wrap(anchor + delta, n)
        reverse = design_candidate(
            sequence,
            reverse_position,
            "Reverse",
            selected,
            masks,
            name_prefix=name_prefix,
            settings=settings,
        )
        if reverse and anchor in interval_positions(reverse.cover, n):
            key = (reverse.direction, reverse.position)
            if key not in seen:
                seen.add(key)
                options.append(reverse)
                if len(options) >= limit:
                    return options
    return options


def design_pair_options_from_anchor(
    sequence: str,
    anchor: int,
    selected: Sequence[Primer],
    masks: Sequence[Interval],
    *,
    limit: int = 6,
    name_prefix: str = DEFAULT_PRIMER_NAME_PREFIX,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> list[tuple[Primer, Primer]]:
    n = len(sequence)
    options: list[tuple[Primer, Primer]] = []
    seen: set[tuple[str, int, str, int]] = set()
    used_primary_positions: list[tuple[str, int]] = []
    search_span = min(n, settings.preferred_pair_max + 120)
    for offset in range(search_span):
        forward_position = wrap(anchor + offset, n)
        forward = design_candidate(
            sequence,
            forward_position,
            "Forward",
            selected,
            masks,
            name_prefix=name_prefix,
            settings=settings,
        )
        if forward and position_is_diverse(forward, used_primary_positions, n):
            for reverse in design_reverse_partner_options(
                sequence,
                forward,
                [*selected, forward],
                masks,
                limit=1,
                name_prefix=name_prefix,
                settings=settings,
            ):
                key = (forward.direction, forward.position, reverse.direction, reverse.position)
                if key not in seen:
                    seen.add(key)
                    used_primary_positions.append((forward.direction, forward.position))
                    options.append((forward, reverse))
                    if len(options) >= limit:
                        return options

        reverse_position = wrap(anchor + offset, n)
        reverse = design_candidate(
            sequence,
            reverse_position,
            "Reverse",
            selected,
            masks,
            name_prefix=name_prefix,
            settings=settings,
        )
        if reverse and position_is_diverse(reverse, used_primary_positions, n):
            for forward_partner in design_forward_partner_options(
                sequence,
                reverse,
                [*selected, reverse],
                masks,
                limit=1,
                name_prefix=name_prefix,
                settings=settings,
            ):
                key = (reverse.direction, reverse.position, forward_partner.direction, forward_partner.position)
                if key not in seen:
                    seen.add(key)
                    used_primary_positions.append((reverse.direction, reverse.position))
                    options.append((reverse, forward_partner))
                    if len(options) >= limit:
                        return options
    return options


def position_is_diverse(
    primer: Primer,
    used_positions: Sequence[tuple[str, int]],
    sequence_length: int,
) -> bool:
    return all(
        direction != primer.direction
        or circular_base_distance(position, primer.position, sequence_length) >= DESIGN_DIVERSE_POSITION_GAP
        for direction, position in used_positions
    )


def design_reverse_partner(
    sequence: str,
    forward: Primer,
    selected: Sequence[Primer],
    masks: Sequence[Interval],
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> Primer | None:
    options = design_reverse_partner_options(sequence, forward, selected, masks, limit=1, settings=settings)
    return options[0] if options else None


def design_reverse_partner_options(
    sequence: str,
    forward: Primer,
    selected: Sequence[Primer],
    masks: Sequence[Interval],
    *,
    limit: int = 3,
    name_prefix: str = DEFAULT_PRIMER_NAME_PREFIX,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> list[Primer]:
    n = len(sequence)
    preferred_positions = [
        wrap(forward.position + delta, n)
        for delta in range(settings.preferred_pair_min, settings.preferred_pair_max + 1)
    ]
    fallback_positions = [
        wrap(forward.position + delta, n)
        for delta in range(settings.min_binding_gap + 1, settings.read_length + 1)
    ]
    seen: set[int] = set()
    options: list[Primer] = []
    for position in [*preferred_positions, *fallback_positions]:
        if position in seen:
            continue
        seen.add(position)
        candidate = design_candidate(
            sequence,
            position,
            "Reverse",
            selected,
            masks,
            name_prefix=name_prefix,
            settings=settings,
        )
        if candidate:
            options.append(candidate)
            if len(options) >= limit:
                return options
    return options


def design_forward_partner_options(
    sequence: str,
    reverse: Primer,
    selected: Sequence[Primer],
    masks: Sequence[Interval],
    *,
    limit: int = 3,
    name_prefix: str = DEFAULT_PRIMER_NAME_PREFIX,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> list[Primer]:
    n = len(sequence)
    preferred_positions = [
        wrap(reverse.position + delta, n)
        for delta in range(settings.preferred_pair_min, settings.preferred_pair_max + 1)
    ]
    fallback_positions = [
        wrap(reverse.position + delta, n)
        for delta in range(settings.min_binding_gap + 1, settings.read_length + 1)
    ]
    seen: set[int] = set()
    options: list[Primer] = []
    for position in [*preferred_positions, *fallback_positions]:
        if position in seen:
            continue
        seen.add(position)
        candidate = design_candidate(
            sequence,
            position,
            "Forward",
            selected,
            masks,
            name_prefix=name_prefix,
            settings=settings,
        )
        if candidate:
            options.append(candidate)
            if len(options) >= limit:
                return options
    return options


def design_candidate(
    sequence: str,
    position: int,
    direction: str,
    selected: Sequence[Primer],
    masks: Sequence[Interval],
    *,
    name_prefix: str = DEFAULT_PRIMER_NAME_PREFIX,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> Primer | None:
    n = len(sequence)
    for length in range(settings.min_primer_len, settings.max_primer_len + 1):
        if direction == "Forward":
            binding = Interval(position, wrap(position + length - 1, n))
            primer_sequence = circular_subsequence(sequence, position, length)
        else:
            start = wrap(position - length + 1, n)
            binding = Interval(start, position)
            primer_sequence = reverse_complement(circular_subsequence(sequence, start, length))

        if interval_overlaps_any(binding, masks, n):
            continue
        if not binding_respects_spacing(binding, selected, n, settings):
            continue
        if not primer_quality_ok(primer_sequence, settings):
            continue
        hits = find_bindings(sequence, primer_sequence)
        expected = (direction, binding.start, binding.end)
        if len(hits) != 1 or (hits[0][0], hits[0][1].start, hits[0][1].end) != expected:
            continue
        cover = depth_interval(position, direction, n, settings)
        primer = Primer(
            name=f"{normalize_primer_name_prefix(name_prefix)}_{'F' if direction == 'Forward' else 'R'}_{position}",
            sequence=primer_sequence.lower(),
            direction=direction,
            position=position,
            binding=binding,
            cover=cover,
            source="designed",
            order=10_000 + position,
        )
        return primer
    return None


def primer_quality_ok(
    sequence: str,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> bool:
    gc = gc_percent(sequence)
    return settings.min_gc <= gc <= settings.max_gc and settings.min_tm <= calc_tm(sequence) <= settings.max_tm


def calc_tm(sequence: str) -> float:
    return _calc_tm_cached(sequence.upper())


@lru_cache(maxsize=100_000)
def _calc_tm_cached(sequence: str) -> float:
    try:
        import primer3
    except ImportError as exc:
        raise Primer3UnavailableError(
            "primer3-py is required for Tm calculation. Install dependencies with uv sync."
        ) from exc
    if hasattr(primer3, "calc_tm"):
        return float(primer3.calc_tm(sequence))
    return float(primer3.bindings.calc_tm(sequence))


def optimize_primers(sequence: str, primers: Sequence[Primer], target_depth: int) -> list[Primer]:
    selected = list(primers)
    n = len(sequence)
    selected = remove_redundant_primers_by_source(n, selected, target_depth, "designed")

    existing = [primer for primer in selected if primer.source == "existing"]
    fixed = [primer for primer in selected if primer.source != "existing"]
    if existing:
        selected = optimize_existing_primer_layout(n, existing, fixed, target_depth)
    return selected


def remove_redundant_primers_by_source(
    sequence_length: int,
    primers: Sequence[Primer],
    target_depth: int,
    source: str,
) -> list[Primer]:
    selected = list(primers)
    changed = True
    while changed:
        changed = False
        for primer in list(selected):
            if primer.source != source:
                continue
            trial = [p for p in selected if p is not primer]
            if min(compute_depth(sequence_length, trial)) >= target_depth:
                selected = trial
                changed = True
                break
    return selected


def optimize_existing_primer_layout(
    sequence_length: int,
    existing: Sequence[Primer],
    fixed: Sequence[Primer],
    target_depth: int,
    beam_width: int = 512,
) -> list[Primer]:
    if min(compute_depth(sequence_length, [*fixed, *existing])) < target_depth:
        return [*fixed, *existing]

    scorer = ExistingLayoutScorer(sequence_length, existing, fixed, target_depth)
    initial_state = tuple(range(len(existing)))
    best_state = initial_state
    best_score = scorer.score(initial_state)
    states = {initial_state}

    while states:
        candidates: dict[tuple[int, ...], tuple[int | float, ...]] = {}
        for state in states:
            for primer_index in state:
                next_state = tuple(index for index in state if index != primer_index)
                if next_state in candidates:
                    continue
                score = scorer.score(next_state)
                if score is not None:
                    candidates[next_state] = score
        if not candidates:
            break

        ranked_states = sorted(candidates, key=candidates.__getitem__)[:beam_width]
        states = set(ranked_states)
        if candidates[ranked_states[0]] < best_score:
            best_state = ranked_states[0]
            best_score = candidates[ranked_states[0]]

    return [*fixed, *(existing[index] for index in best_state)]


class ExistingLayoutScorer:
    def __init__(
        self,
        sequence_length: int,
        existing: Sequence[Primer],
        fixed: Sequence[Primer],
        target_depth: int,
    ) -> None:
        self.sequence_length = sequence_length
        self.existing = tuple(existing)
        self.fixed = tuple(fixed)
        self.target_depth = target_depth
        self.segments = depth_segments(sequence_length, [*self.fixed, *self.existing])
        self.segment_lengths = tuple(end - start + 1 for start, end in self.segments)
        self.fixed_depth = self._depth_vector(self.fixed)
        self.existing_depth = tuple(self._depth_vector((primer,)) for primer in self.existing)
        self._cache: dict[tuple[int, ...], tuple[int | float, ...] | None] = {}

    def score(self, state: tuple[int, ...]) -> tuple[int | float, ...] | None:
        if state in self._cache:
            return self._cache[state]

        depth = list(self.fixed_depth)
        for index in state:
            primer_depth = self.existing_depth[index]
            for segment_index, value in enumerate(primer_depth):
                depth[segment_index] += value

        if any(value < self.target_depth for value in depth):
            self._cache[state] = None
            return None

        primers = [*self.fixed, *(self.existing[index] for index in state)]
        placement = placement_score(self.sequence_length, primers)
        overage_squared = sum(
            self.segment_lengths[index] * (value - self.target_depth) ** 2
            for index, value in enumerate(depth)
        )
        high_depth_bases = sum(
            self.segment_lengths[index]
            for index, value in enumerate(depth)
            if value >= self.target_depth + 2
        )
        order_sum = sum(primer.order for primer in primers)
        score = (
            placement["same_direction_adjacent"],
            overage_squared,
            high_depth_bases,
            placement["max_adjacent_gap"],
            len(primers),
            placement["fr_gap_variance"],
            order_sum,
        )
        self._cache[state] = score
        return score

    def _depth_vector(self, primers: Sequence[Primer]) -> list[int]:
        depth = [0] * len(self.segments)
        for primer in primers:
            for index, (start, end) in enumerate(self.segments):
                if intervals_overlap(primer.cover, Interval(start, end), self.sequence_length):
                    depth[index] += 1
        return depth


def depth_segments(sequence_length: int, primers: Sequence[Primer]) -> tuple[tuple[int, int], ...]:
    breakpoints = {1, sequence_length + 1}
    for primer in primers:
        for start, end in interval_segments(primer.cover, sequence_length):
            breakpoints.add(start)
            breakpoints.add(end + 1)
    ordered = sorted(point for point in breakpoints if 1 <= point <= sequence_length + 1)
    return tuple((start, end - 1) for start, end in zip(ordered, ordered[1:]) if start <= end - 1)


def placement_score(sequence_length: int, primers: Sequence[Primer]) -> dict[str, int | float]:
    if len(primers) < 2:
        return {
            "same_direction_adjacent": 0,
            "max_adjacent_gap": sequence_length,
            "fr_gap_variance": 0.0,
        }

    ordered = sorted(primers, key=lambda primer: primer.position)
    same_direction_adjacent = 0
    adjacent_gaps: list[int] = []
    forward_reverse_gaps: list[int] = []
    for index, primer in enumerate(ordered):
        next_primer = ordered[(index + 1) % len(ordered)]
        if primer.direction == next_primer.direction:
            same_direction_adjacent += 1
        gap = (next_primer.position - primer.position) % sequence_length or sequence_length
        adjacent_gaps.append(gap)
        if primer.direction == "Forward" and next_primer.direction == "Reverse":
            forward_reverse_gaps.append(gap)

    if forward_reverse_gaps:
        center = median(forward_reverse_gaps)
        fr_gap_variance = sum((gap - center) ** 2 for gap in forward_reverse_gaps)
    else:
        fr_gap_variance = float(sequence_length**2)

    return {
        "same_direction_adjacent": same_direction_adjacent,
        "max_adjacent_gap": max(adjacent_gaps),
        "fr_gap_variance": fr_gap_variance,
    }


def compute_depth(sequence_length: int, primers: Sequence[Primer]) -> list[int]:
    depth = [0] * sequence_length
    for primer in primers:
        add_depth(depth, primer.cover, 1)
    return depth


def add_depth(depth: list[int], interval: Interval, value: int) -> None:
    for pos in interval_positions(interval, len(depth)):
        depth[pos - 1] += value


def depth_gain(depth: Sequence[int], primer: Primer, target_depth: int, sequence_length: int) -> int:
    return sum(1 for pos in interval_positions(primer.cover, sequence_length) if depth[pos - 1] < target_depth)


def low_depth_gain(depth: Sequence[int], primer: Primer, sequence_length: int) -> int:
    min_depth = min(depth) if depth else 0
    return sum(1 for pos in interval_positions(primer.cover, sequence_length) if depth[pos - 1] == min_depth)


def depth_regions_below(depth: Sequence[int], target_depth: int) -> list[Interval]:
    positions = [index + 1 for index, value in enumerate(depth) if value < target_depth]
    return positions_to_intervals(positions, len(depth))


def positions_below_depth(depth: Sequence[int], target_depth: int) -> list[int]:
    min_depth = min(depth)
    return [index + 1 for index, value in enumerate(depth) if value == min_depth and value < target_depth]


def positions_to_intervals(positions: Sequence[int], sequence_length: int) -> list[Interval]:
    if not positions:
        return []
    position_set = set(positions)
    if len(position_set) == sequence_length:
        return [Interval(1, sequence_length)]

    starts = [pos for pos in sorted(position_set) if wrap(pos - 1, sequence_length) not in position_set]
    intervals = []
    for start in starts:
        end = start
        while wrap(end + 1, sequence_length) in position_set:
            end = wrap(end + 1, sequence_length)
        intervals.append(Interval(start, end))
    return intervals


def depth_interval(
    position: int,
    direction: str,
    sequence_length: int,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> Interval:
    if direction == "Forward":
        return Interval(
            wrap(position + settings.noise_length, sequence_length),
            wrap(position + settings.read_length - 1, sequence_length),
        )
    return Interval(
        wrap(position - settings.read_length + 1, sequence_length),
        wrap(position - settings.noise_length, sequence_length),
    )


def respects_spacing(
    primer: Primer,
    selected: Sequence[Primer],
    sequence_length: int,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> bool:
    return binding_respects_spacing(primer.binding, selected, sequence_length, settings)


def binding_respects_spacing(
    binding: Interval,
    selected: Sequence[Primer],
    sequence_length: int,
    settings: DesignSettings = DEFAULT_DESIGN_SETTINGS,
) -> bool:
    return all(binding_gap(binding, other.binding, sequence_length) >= settings.min_binding_gap for other in selected)


def binding_gap(a: Interval, b: Interval, sequence_length: int) -> int:
    best = sequence_length
    for a_start, a_end in interval_segments(a, sequence_length):
        for b_start, b_end in interval_segments(b, sequence_length):
            best = min(best, segment_gap(a_start, a_end, b_start, b_end, sequence_length))
            if best == 0:
                return 0
    return best


def segment_gap(a_start: int, a_end: int, b_start: int, b_end: int, sequence_length: int) -> int:
    if a_start <= b_end and b_start <= a_end:
        return 0
    if a_end < b_start:
        linear_gap = b_start - a_end - 1
        circular_gap = a_start + sequence_length - b_end - 1
        return min(linear_gap, circular_gap)
    linear_gap = a_start - b_end - 1
    circular_gap = b_start + sequence_length - a_end - 1
    return min(linear_gap, circular_gap)


def circular_base_distance(a: int, b: int, sequence_length: int) -> int:
    diff = abs(a - b)
    return min(diff, sequence_length - diff)


def interval_overlaps_any(interval: Interval, others: Sequence[Interval], sequence_length: int) -> bool:
    return any(intervals_overlap(interval, other, sequence_length) for other in others)


def intervals_overlap(a: Interval, b: Interval, sequence_length: int) -> bool:
    for a_start, a_end in interval_segments(a, sequence_length):
        for b_start, b_end in interval_segments(b, sequence_length):
            if a_start <= b_end and b_start <= a_end:
                return True
    return False


def interval_segments(interval: Interval, sequence_length: int) -> list[tuple[int, int]]:
    start = normalize_position(interval.start, sequence_length)
    end = normalize_position(interval.end, sequence_length)
    if start <= end:
        return [(start, end)]
    return [(start, sequence_length), (1, end)]


def interval_positions(interval: Interval, sequence_length: int) -> list[int]:
    start = normalize_position(interval.start, sequence_length)
    end = normalize_position(interval.end, sequence_length)
    if start <= end:
        return list(range(start, end + 1))
    return [*range(start, sequence_length + 1), *range(1, end + 1)]


def circular_subsequence(sequence: str, start: int, length: int) -> str:
    n = len(sequence)
    return "".join(sequence[wrap(start + offset, n) - 1] for offset in range(length))


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def gc_percent(sequence: str) -> float:
    seq = sequence.upper()
    if not seq:
        return 0.0
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def normalize_interval(interval: Interval, sequence_length: int) -> Interval:
    return Interval(normalize_position(interval.start, sequence_length), normalize_position(interval.end, sequence_length))


def normalize_position(position: int, sequence_length: int) -> int:
    if not 1 <= position <= sequence_length:
        raise ValueError(f"position {position} is outside sequence length {sequence_length}")
    return position


def wrap(position: int, sequence_length: int) -> int:
    return ((position - 1) % sequence_length) + 1


def format_interval(interval: Interval) -> str:
    return f"{interval.start}..{interval.end}"


def write_primer_list(path: str | Path, primers: Iterable[Primer], delimiter: str = "\t") -> None:
    Path(path).write_text(format_primer_list(primers, delimiter), encoding="utf-8")


def format_primer_list(primers: Iterable[Primer], delimiter: str = "\t") -> str:
    lines = [delimiter.join([primer.name, primer.sequence, primer.output_memo]) for primer in primers]
    return "\n".join(lines) + ("\n" if lines else "")


def delimiter_from_name(name: str) -> str:
    if name == "tab":
        return "\t"
    if name == "semicolon":
        return ";"
    if name == "comma":
        return ","
    raise ValueError(f"unsupported delimiter: {name}")

