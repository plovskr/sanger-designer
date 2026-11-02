import csv
import io
from pathlib import Path
import zipfile

import pytest

from sanger_designer.core import (
    DesignResult,
    DesignSettings,
    Interval,
    Primer,
    PrimerInput,
    format_primer_list,
    parse_primer_list_text,
)
from sanger_designer.diagnostics import build_support_bundle, format_primer_list_from_inputs


@pytest.fixture
def primer() -> Primer:
    return Primer(
        name="existing",
        sequence="ACGTACGTACGTACGTACGT",
        direction="Forward",
        position=1,
        binding=Interval(1, 20),
        cover=Interval(1, 600),
        source="existing",
        input_memo='one, two; three\tfour "quoted"\nfive',
    )


@pytest.mark.parametrize("delimiter", ["\t", ",", ";"])
def test_output_memo_is_quoted_and_can_be_read_back(primer: Primer, delimiter: str) -> None:
    text = format_primer_list([primer], delimiter)
    assert text.startswith(primer.name + delimiter + primer.sequence + delimiter + '"')
    assert '""quoted""' in text
    assert list(csv.reader(io.StringIO(text), delimiter=delimiter)) == [
        [primer.name, primer.sequence, primer.output_memo]
    ]
    restored = parse_primer_list_text(text)
    assert len(restored) == 1
    assert restored[0].name == primer.name
    assert restored[0].sequence == primer.sequence
    assert restored[0].memo == primer.output_memo


@pytest.mark.parametrize("delimiter", ["\t", ",", ";"])
@pytest.mark.parametrize("memo", ["", "plain", "one,two", "one;two", "one\ttwo", 'one "two"'])
def test_intermediate_memo_is_always_quoted(delimiter: str, memo: str) -> None:
    inputs = [PrimerInput("existing", "ACGT", memo)]
    text = format_primer_list_from_inputs(inputs, delimiter)
    assert text == f'existing{delimiter}ACGT{delimiter}"{memo.replace(chr(34), chr(34) * 2)}"\n'
    assert parse_primer_list_text(text) == inputs


@pytest.mark.parametrize("delimiter", ["\t", ",", ";"])
def test_optional_memos_blank_lines_and_comments(delimiter: str) -> None:
    text = f"\n# comment, with; other\tseparators\nfirst{delimiter}ac 1gt\n\nsecond{delimiter}TGCA{delimiter}plain\n"
    assert parse_primer_list_text(text) == [
        PrimerInput("first", "ACGT", order=2),
        PrimerInput("second", "TGCA", "plain", order=4),
    ]


@pytest.mark.parametrize("delimiter", ["\t", ",", ";"])
def test_quoted_names_and_multiline_memos(delimiter: str) -> None:
    inputs = [
        PrimerInput('name,with;tab\tand"quote', "ACGT", 'first\nsecond,third;fourth\tfifth"sixth'),
        PrimerInput("other", "TGCA", "ordinary"),
    ]
    text = format_primer_list_from_inputs(inputs, delimiter)
    restored = parse_primer_list_text(text)
    assert [(row.name, row.sequence, row.memo) for row in restored] == [
        (row.name, row.sequence, row.memo) for row in inputs
    ]


@pytest.mark.parametrize("delimiter", ["\t", ",", ";"])
def test_support_bundle_uses_the_same_quoted_format(primer: Primer, delimiter: str) -> None:
    merged = [PrimerInput(primer.name, primer.sequence, primer.input_memo)]
    input_text = format_primer_list_from_inputs(merged, delimiter)
    result = DesignResult((primer,), (1,) * 600, 2, (Interval(1, 600),), False)
    data = build_support_bundle(
        result=result,
        sequence_length=600,
        settings=DesignSettings(),
        run_parameters={},
        sequence_file_name="plasmid.fasta",
        sequence_text=">plasmid\n" + "ACGT" * 150,
        uploaded_primer_files=[{"name": "primers.txt", "text": input_text}],
        default_primer_text=None,
        merged_primers=merged,
        masks=(),
        warnings=(),
        runtime_seconds=0.1,
        delimiter=delimiter,
    )
    with zipfile.ZipFile(io.BytesIO(data)) as bundle:
        intermediate = bundle.read("intermediate/merged_primers.txt").decode()
        output = bundle.read("results/sanger_primers.txt").decode()
        run_log = bundle.read("logs/run_log.txt").decode()
    assert intermediate == format_primer_list_from_inputs(merged, delimiter)
    assert output == format_primer_list([primer], delimiter)
    assert parse_primer_list_text(intermediate)[0].memo == primer.input_memo
    assert parse_primer_list_text(output)[0].memo == primer.output_memo
    assert "Input primer rows: 1\n" in run_log
    assert "Duplicate primer sequences skipped: 0\n" in run_log


@pytest.mark.parametrize("text", ["", "\n# comment\n\n", "  \t \n"])
def test_empty_primer_list(text: str) -> None:
    assert parse_primer_list_text(text) == []
    assert format_primer_list([]) == ""


@pytest.mark.parametrize(
    "text",
    ["invalid row", "primer,NNNN", "primer,ACGT,unquoted,memo", "ok,ACGT\nbad,TGCA,memo,extra"],
)
def test_invalid_primer_rows_are_rejected(text: str) -> None:
    with pytest.raises(ValueError):
        parse_primer_list_text(text)


@pytest.mark.parametrize("text", ["first,ACGT\nsecond;TGCA\n", "first;ACGT\nsecond\tTGCA\n"])
def test_mixed_delimiters_are_rejected(text: str) -> None:
    with pytest.raises(ValueError):
        parse_primer_list_text(text)


def test_local_default_primer_list_is_readable(tmp_path: Path) -> None:
    path = tmp_path / "default_primers.txt"
    path.write_text("default_forward\tacgt\ndefault_reverse\tTGCA\n", encoding="utf-8")
    inputs = parse_primer_list_text(path.read_text(encoding="utf-8"))
    assert inputs == [
        PrimerInput("default_forward", "ACGT"),
        PrimerInput("default_reverse", "TGCA", order=1),
    ]
