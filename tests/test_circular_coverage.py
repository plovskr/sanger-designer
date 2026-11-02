import random

import pytest

from sanger_designer.core import (
    DesignResult,
    DesignSettings,
    Interval,
    Primer,
    PrimerInput,
    compute_depth,
    depth_interval,
    design_primers,
    interval_positions,
    reverse_complement,
)
from sanger_designer.gui import build_primer_placement_figure


def make_primer(direction: str, position: int, length: int, settings: DesignSettings) -> Primer:
    return Primer(
        name="existing",
        sequence="ACGTACGTACGTACGTACGT",
        direction=direction,
        position=position,
        binding=Interval(position, position),
        cover=depth_interval(position, direction, length, settings),
        source="existing",
    )


@pytest.mark.parametrize("length", [200, 600, 699, 700, 701, 1000])
@pytest.mark.parametrize("direction", ["Forward", "Reverse"])
@pytest.mark.parametrize("near_end", [False, True])
def test_coverage_matches_unique_read_positions(length: int, direction: str, near_end: bool) -> None:
    settings = DesignSettings()
    position = length - 20 if near_end else 1
    primer = make_primer(direction, position, length, settings)
    sign = 1 if direction == "Forward" else -1
    expected = {
        (position - 1 + sign * offset) % length + 1
        for offset in range(settings.noise_length, settings.read_length)
    }
    assert set(interval_positions(primer.cover, length)) == expected
    assert len(expected) == min(length, settings.effective_read_length)
    assert compute_depth(length, [primer]) == [int(position in expected) for position in range(1, length + 1)]
    if length <= settings.effective_read_length:
        assert primer.cover == Interval(1, length)


@pytest.mark.parametrize("direction", ["Forward", "Reverse"])
def test_multiple_laps_contribute_only_one_depth_per_primer(direction: str) -> None:
    settings = DesignSettings(read_length=2000, noise_length=100)
    primer = make_primer(direction, 500, 600, settings)
    assert primer.cover == Interval(1, 600)
    assert compute_depth(600, [primer]) == [1] * 600
    assert compute_depth(600, [primer, primer]) == [2] * 600


@pytest.mark.parametrize("length", [600, 699, 700])
def test_short_plasmid_design_achieves_depth_with_two_existing_primers(length: int) -> None:
    generator = random.Random(42)
    sequence = "".join(generator.choices("ACGT", k=length))
    inputs = [
        PrimerInput("forward", sequence[49:69]),
        PrimerInput("reverse", reverse_complement(sequence[299:319])),
    ]
    result = design_primers(sequence, inputs, target_depth=2, max_iterations=1)
    assert result.achieved
    assert not result.missing_regions
    assert len(result.primers) == 2
    assert {primer.direction for primer in result.primers} == {"Forward", "Reverse"}
    assert all(primer.cover == Interval(1, length) for primer in result.primers)
    assert result.depth == (2,) * length


@pytest.mark.parametrize(
    ("direction", "position", "noise", "expected"),
    [
        ("Forward", 1, 100, [(1, 101)]),
        ("Forward", 580, 100, [(580, 600), (1, 80)]),
        ("Reverse", 1, 100, [(501, 600), (1, 1)]),
        ("Reverse", 580, 100, [(480, 580)]),
        ("Forward", 200, 50, [(200, 250)]),
        ("Reverse", 200, 50, [(150, 200)]),
        ("Forward", 200, 0, []),
        ("Reverse", 200, 0, []),
        ("Forward", 200, 600, [(1, 600)]),
        ("Reverse", 200, 600, [(1, 600)]),
    ],
)
def test_full_circle_placement_keeps_noise_based_on_primer_position(
    direction: str, position: int, noise: int, expected: list[tuple[int, int]]
) -> None:
    settings = DesignSettings(read_length=1500, noise_length=noise)
    primer = make_primer(direction, position, 600, settings)
    result = DesignResult((primer,), (1,) * 600, 2, (Interval(1, 600),), False)
    figure = build_primer_placement_figure(result, 600, (), True, settings)
    dotted = [(shape.x0, shape.x1) for shape in figure.layout.shapes if shape.line.dash == "dot"]
    assert dotted == expected
    figure = build_primer_placement_figure(result, 600, (), False, settings)
    assert not any(shape.line.dash == "dot" for shape in figure.layout.shapes)
