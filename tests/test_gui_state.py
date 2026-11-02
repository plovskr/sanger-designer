import random

import pytest
from streamlit.testing.v1 import AppTest

from sanger_designer import defaults


@pytest.fixture
def app() -> AppTest:
    generator = random.Random(42)
    sequence = "".join(generator.choices("ACGT", k=1000))
    primer_text = "\n".join(
        f"existing_{position}\t{sequence[position - 1:position + 19]}"
        for position in (100, 400, 700)
    )
    app = AppTest.from_file("src/sanger_designer/gui.py", default_timeout=30).run()
    app.file_uploader[0].set_value(("plasmid.fasta", (">plasmid\n" + sequence).encode(), "text/plain"))
    app.file_uploader[1].set_value([("primers.txt", primer_text.encode(), "text/plain")])
    app.checkbox[0].set_value(False)
    app.run()
    app.button[0].click().run()
    assert not app.exception
    assert app.session_state["result"].achieved
    assert len(app.get("download_button")) == 3
    return app


@pytest.mark.parametrize(
    ("widget_type", "label", "value"),
    [
        ("file_uploader", "Plasmid sequence", ("replacement.fasta", b">replacement\nACGT", "text/plain")),
        ("file_uploader", "Plasmid sequence", ("plasmid.fasta", b">plasmid\nACGT", "text/plain")),
        ("file_uploader", "Plasmid sequence", None),
        ("file_uploader", "Existing primer lists", [("replacement.txt", b"other\tACGT", "text/plain")]),
        ("file_uploader", "Existing primer lists", None),
        ("checkbox", "Use default_primers.txt", True),
        ("button_group", "Target depth", 3),
        ("button_group", "Target depth", 4),
        ("text_input", "Mask regions", "1..30"),
        ("text_input", "New primer prefix", "replacement"),
        ("selectbox", "Output delimiter", "comma"),
        ("number_input", "Max iterations", 210),
        ("number_input", "Read length", 850),
        ("number_input", "Noisy bases from primer", 110),
        ("number_input", "Minimum binding gap", 160),
        ("number_input", "Preferred pair min", 310),
        ("number_input", "Preferred pair max", 410),
        ("number_input", "Primer length min", 19),
        ("number_input", "Primer length max", 26),
        ("number_input", "Tm min", 54.5),
        ("number_input", "Tm max", 65.5),
        ("number_input", "GC min (%)", 39.0),
        ("number_input", "GC max (%)", 61.0),
        ("number_input", "Read length", 50),
        ("number_input", "Noisy bases from primer", 1000),
        ("number_input", "Preferred pair min", 410),
    ],
)
def test_input_changes_clear_previous_result(app: AppTest, widget_type: str, label: str, value) -> None:
    widget = next(item for item in app.get(widget_type) if item.label == label)
    widget.set_value(value).run()
    assert not app.exception
    assert not app.metric
    assert not app.get("download_button")
    assert app.info
    for key in (
        "result", "sequence_text", "sequence_file_name", "run_parameters",
        "settings", "uploaded_primer_files", "default_primer_text",
    ):
        assert key not in app.session_state
    assert app.session_state["run_requested"] is False


@pytest.mark.parametrize(
    ("widget_type", "label", "value"),
    [
        ("checkbox", "Show cover bars", False),
        ("button_group", "Source", "Existing"),
        ("button_group", "Source", "Designed"),
    ],
)
def test_display_changes_keep_result(app: AppTest, widget_type: str, label: str, value) -> None:
    result = app.session_state["result"]
    widget = next(item for item in app.get(widget_type) if item.label == label)
    widget.set_value(value).run()
    assert not app.exception
    assert app.session_state["result"] == result
    assert app.metric
    assert len(app.get("download_button")) == 3


def test_target_depth_requires_selection() -> None:
    app = AppTest.from_file("src/sanger_designer/gui.py").run()
    depth = app.get("button_group")[0]
    assert depth.required is True
    assert depth.value == 2
    depth.set_value(4).run()
    assert not app.exception
    assert app.get("button_group")[0].value == 4
    app.get("button_group")[0].set_value(4).run()
    assert not app.exception
    assert app.get("button_group")[0].value == 4


def test_reverting_input_does_not_restore_previous_result(app: AppTest) -> None:
    app.text_input[1].set_value("replacement").run()
    app.text_input[1].set_value("designed").run()
    assert not app.exception
    assert "result" not in app.session_state
    assert not app.metric
    assert not app.get("download_button")


def test_design_can_run_again_after_input_change(app: AppTest) -> None:
    app.text_input[1].set_value("replacement").run()
    assert not app.get("download_button")
    app.button[0].click().run()
    assert not app.exception
    assert app.session_state["result"].achieved
    assert app.session_state["run_parameters"]["prefix"] == "replacement"
    assert app.session_state["sequence_file_name"] == "plasmid.fasta"
    assert len(app.get("download_button")) == 3


def test_missing_local_default_list_does_not_block_design(app: AppTest, monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(defaults, "DEFAULT_PRIMER_LIST_PATH", tmp_path / "default_primers.txt")
    use_default = next(item for item in app.checkbox if item.label == "Use default_primers.txt")
    use_default.set_value(True).run()
    app.button[0].click().run()
    assert not app.exception
    assert app.session_state["result"].achieved
    assert any("Default primer list not found:" in warning.value for warning in app.warning)
    assert len(app.get("download_button")) == 3


def test_upload_does_not_execute_a_previous_request_without_input() -> None:
    app = AppTest.from_file("src/sanger_designer/gui.py").run()
    app.button[0].click().run()
    app.file_uploader[0].set_value(("plasmid.fasta", b">plasmid\nACGT", "text/plain")).run()
    assert not app.exception
    assert app.session_state["run_requested"] is False
    assert "result" not in app.session_state
    assert not app.get("download_button")
