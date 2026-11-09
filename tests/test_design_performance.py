from dataclasses import replace
from itertools import product
from random import Random

import pytest

from sanger_designer import core


def reference_candidate(sequence, position, direction, selected, masks, settings):
    n = len(sequence)
    for length in range(settings.min_primer_len, settings.max_primer_len + 1):
        if direction == "Forward":
            binding = core.Interval(position, core.wrap(position + length - 1, n))
            primer_sequence = core.circular_subsequence(sequence, position, length)
        else:
            start = core.wrap(position - length + 1, n)
            binding = core.Interval(start, position)
            primer_sequence = core.reverse_complement(core.circular_subsequence(sequence, start, length))
        if core.interval_overlaps_any(binding, masks, n):
            continue
        if not core.binding_respects_spacing(binding, selected, n, settings):
            continue
        if not core.primer_quality_ok(primer_sequence, settings):
            continue
        hits = core.find_bindings(sequence, primer_sequence)
        if len(hits) != 1 or (hits[0][0], hits[0][1].start, hits[0][1].end) != (direction, binding.start, binding.end):
            continue
        return core.Primer(
            name=f"{core.DEFAULT_PRIMER_NAME_PREFIX}_{'F' if direction == 'Forward' else 'R'}_{position}",
            sequence=primer_sequence.lower(),
            direction=direction,
            position=position,
            binding=binding,
            cover=core.depth_interval(position, direction, n, settings),
            source="designed",
            order=10_000 + position,
        )
    return None


@pytest.mark.parametrize("n", range(1, 9))
def test_binding_gap_matches_base_distance_on_all_circular_intervals(n):
    intervals = [core.Interval(start, end) for start, end in product(range(1, n + 1), repeat=2)]
    positions = {interval: set(core.interval_positions(interval, n)) for interval in intervals}
    for a, b in product(intervals, repeat=2):
        if positions[a] & positions[b]:
            expected = 0
        else:
            expected = min(min(abs(x - y), n - abs(x - y)) for x in positions[a] for y in positions[b]) - 1
        assert core.binding_gap(a, b, n) == expected, (n, a, b)


@pytest.mark.parametrize("endpoint", range(4))
@pytest.mark.parametrize("value", [0, 11])
def test_binding_gap_rejects_out_of_range_endpoints(endpoint, value):
    points = [1, 2, 5, 7]
    points[endpoint] = value
    with pytest.raises(ValueError, match="outside sequence length"):
        core.binding_gap(core.Interval(*points[:2]), core.Interval(*points[2:]), 10)


@pytest.mark.parametrize("n", [7, 23, 80, 131])
@pytest.mark.parametrize("gap", [0, 1, 10])
def test_candidate_pruning_preserves_results(monkeypatch, n, gap):
    random = Random(n)
    sequence = "".join(random.choices("ACGT", k=n))
    settings = replace(core.DEFAULT_DESIGN_SETTINGS, min_primer_len=8, max_primer_len=25, min_binding_gap=gap)
    monkeypatch.setattr(core, "primer_quality_ok", lambda seq, settings: len(seq) >= 12)
    existing = core.Primer(
        name="synthetic", sequence="acgt", direction="Forward", position=n - 1,
        binding=core.Interval(n - 1, 2), cover=core.Interval(1, n), source="existing",
    )
    for selected, masks in product([(), (existing,)], [(), (core.Interval(3, 5),), (core.Interval(n - 2, 2),)]):
        for direction, position in product(["Forward", "Reverse"], range(1, n + 1)):
            expected = reference_candidate(sequence, position, direction, selected, masks, settings)
            assert core.design_candidate(sequence, position, direction, selected, masks, settings=settings) == expected


@pytest.mark.parametrize("direction,position", [("Forward", 98), ("Reverse", 3)])
@pytest.mark.parametrize("constraint", ["mask", "spacing"])
def test_blocked_candidate_stops_before_sequence_work(monkeypatch, direction, position, constraint):
    binding = core.Interval(100, 1)
    masks = (binding,) if constraint == "mask" else ()
    selected = () if constraint == "mask" else (
        core.Primer("synthetic", "ac", "Forward", 100, binding, core.Interval(1, 100), "existing"),
    )
    calls = []
    original = core.interval_overlaps_any

    def overlaps(*args):
        calls.append(args[0])
        return original(*args)

    def unexpected_sequence_work(*args):
        pytest.fail("A mask or spacing conflict should stop before extracting a sequence")

    monkeypatch.setattr(core, "interval_overlaps_any", overlaps)
    monkeypatch.setattr(core, "circular_subsequence", unexpected_sequence_work)
    assert core.design_candidate("ACGT" * 25, position, direction, selected, masks) is None
    assert len(calls) == 1


@pytest.mark.parametrize("direction", ["Forward", "Reverse"])
def test_candidate_satisfies_default_quality_constraints(direction):
    sequence = "".join(Random(42).choices("ACGT", k=400))
    candidates = [core.design_candidate(sequence, position, direction, (), ()) for position in range(1, 401)]
    assert any(candidates)
    for candidate in filter(None, candidates):
        assert core.primer_quality_ok(candidate.sequence)
        assert core.find_bindings(sequence, candidate.sequence) == [(direction, candidate.binding, candidate.position)]
