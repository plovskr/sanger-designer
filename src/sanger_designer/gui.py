from __future__ import annotations

from collections import Counter
import time

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from sanger_designer.defaults import (
    DEFAULT_COVERAGE_REPORT_PATH,
    DEFAULT_DELIMITER_NAME,
    DEFAULT_MASK_TEXT,
    DEFAULT_MAX_ITERATIONS,
    DEFAULT_OUTPUT_PRIMER_LIST_PATH,
    DEFAULT_PRIMER_LIST_NAME,
    DEFAULT_PRIMER_LIST_PATH,
    DEFAULT_PRIMER_NAME_PREFIX,
    DEFAULT_TARGET_COVERAGE,
    DEFAULT_USE_DEFAULT_PRIMER_LIST,
    DELIMITER_OPTIONS,
    MASK_TEXT_PLACEHOLDER,
    MAX_ITERATIONS_INPUT_MAX,
    MAX_ITERATIONS_INPUT_MIN,
    MAX_ITERATIONS_INPUT_STEP,
    TARGET_COVERAGE_OPTIONS,
)
from sanger_designer.core import (
    DEFAULT_DESIGN_SETTINGS,
    DesignSettings,
    DesignResult,
    Interval,
    Primer,
    delimiter_from_name,
    design_primers,
    format_interval,
    format_primer_list,
    interval_segments,
    merge_primer_inputs,
    parse_masks,
    parse_primer_list_text,
    read_sequence_text,
)
from sanger_designer.diagnostics import (
    SENSITIVE_BUNDLE_NOTICE,
    build_support_bundle,
    default_bundle_filename,
)


def main() -> None:
    st.set_page_config(page_title="Sanger Primer Designer", layout="wide")
    st.title("Sanger Primer Designer")

    params = render_sidebar()
    if not params["ready"]:
        render_empty_state()
        return

    if st.session_state.get("run_requested"):
        run_design(params)

    result = st.session_state.get("result")
    if result is None:
        render_empty_state()
        return

    render_result(result, st.session_state["sequence_length"], st.session_state["masks"], params["delimiter"])


def render_sidebar() -> dict:
    with st.sidebar:
        st.header("Input")
        sequence_file = st.file_uploader("Plasmid sequence", type=["gb", "gbk", "fa", "fasta", "txt"])
        primer_files = st.file_uploader("Existing primer lists", type=["txt"], accept_multiple_files=True)
        use_default_primers = st.checkbox(
            f"Use {DEFAULT_PRIMER_LIST_NAME}",
            value=DEFAULT_USE_DEFAULT_PRIMER_LIST,
        )

        st.header("Options")
        coverage = st.segmented_control(
            "Target coverage",
            options=list(TARGET_COVERAGE_OPTIONS),
            default=DEFAULT_TARGET_COVERAGE,
        )
        mask_text = st.text_input(
            "Mask regions",
            value=DEFAULT_MASK_TEXT,
            placeholder=MASK_TEXT_PLACEHOLDER,
        )
        prefix = st.text_input("New primer prefix", value=DEFAULT_PRIMER_NAME_PREFIX)
        delimiter_name = st.selectbox(
            "Output delimiter",
            options=list(DELIMITER_OPTIONS),
            index=DELIMITER_OPTIONS.index(DEFAULT_DELIMITER_NAME),
        )
        max_iterations = st.number_input(
            "Max iterations",
            min_value=MAX_ITERATIONS_INPUT_MIN,
            max_value=MAX_ITERATIONS_INPUT_MAX,
            value=DEFAULT_MAX_ITERATIONS,
            step=MAX_ITERATIONS_INPUT_STEP,
        )
        settings = render_advanced_settings()

        run_clicked = st.button("Design primers", type="primary", use_container_width=True)

    validation_errors = validate_gui_settings(settings)
    if validation_errors:
        for error in validation_errors:
            st.sidebar.error(error)

    if run_clicked:
        st.session_state["run_requested"] = True

    return {
        "ready": sequence_file is not None and not validation_errors,
        "sequence_file": sequence_file,
        "primer_files": primer_files,
        "use_default_primers": use_default_primers,
        "coverage": int(coverage),
        "mask_text": mask_text,
        "prefix": prefix,
        "delimiter": delimiter_from_name(delimiter_name),
        "max_iterations": int(max_iterations),
        "settings": settings,
    }


def render_advanced_settings() -> DesignSettings:
    defaults = DEFAULT_DESIGN_SETTINGS
    with st.expander("Advanced settings"):
        read_length = st.number_input(
            "Read length",
            min_value=1,
            max_value=2000,
            value=defaults.read_length,
            step=50,
        )
        noise_length = st.number_input(
            "Noisy bases from primer",
            min_value=0,
            max_value=1000,
            value=defaults.noise_length,
            step=10,
        )
        min_gap = st.number_input(
            "Minimum binding gap",
            min_value=0,
            max_value=1000,
            value=defaults.min_binding_gap,
            step=10,
        )

        pair_cols = st.columns(2)
        preferred_pair_min = pair_cols[0].number_input(
            "Preferred pair min",
            min_value=0,
            max_value=2000,
            value=defaults.preferred_pair_min,
            step=10,
        )
        preferred_pair_max = pair_cols[1].number_input(
            "Preferred pair max",
            min_value=0,
            max_value=2000,
            value=defaults.preferred_pair_max,
            step=10,
        )

        length_cols = st.columns(2)
        min_primer_len = length_cols[0].number_input(
            "Primer length min",
            min_value=1,
            max_value=100,
            value=defaults.min_primer_len,
            step=1,
        )
        max_primer_len = length_cols[1].number_input(
            "Primer length max",
            min_value=1,
            max_value=100,
            value=defaults.max_primer_len,
            step=1,
        )

        tm_cols = st.columns(2)
        min_tm = tm_cols[0].number_input("Tm min", value=defaults.min_tm, step=0.5)
        max_tm = tm_cols[1].number_input("Tm max", value=defaults.max_tm, step=0.5)

        gc_cols = st.columns(2)
        min_gc = gc_cols[0].number_input(
            "GC min (%)",
            min_value=0.0,
            max_value=100.0,
            value=defaults.min_gc,
            step=1.0,
        )
        max_gc = gc_cols[1].number_input(
            "GC max (%)",
            min_value=0.0,
            max_value=100.0,
            value=defaults.max_gc,
            step=1.0,
        )

    return DesignSettings(
        read_length=int(read_length),
        noise_length=int(noise_length),
        min_primer_len=int(min_primer_len),
        max_primer_len=int(max_primer_len),
        min_tm=float(min_tm),
        max_tm=float(max_tm),
        min_gc=float(min_gc),
        max_gc=float(max_gc),
        min_binding_gap=int(min_gap),
        preferred_pair_min=int(preferred_pair_min),
        preferred_pair_max=int(preferred_pair_max),
    )


def validate_gui_settings(settings: DesignSettings) -> list[str]:
    try:
        settings.validate()
    except ValueError as exc:
        return [str(exc)]
    return []


def render_empty_state() -> None:
    st.info("Upload a plasmid sequence file and choose primer settings, then run the design.")


def run_design(params: dict) -> None:
    st.session_state["run_requested"] = False
    start = time.perf_counter()
    try:
        sequence_text = decode_upload(params["sequence_file"])
        sequence = read_sequence_text(sequence_text)
        uploaded_primer_files = []
        warnings = []
        default_primer_text = None
        primer_groups = []
        if params["use_default_primers"]:
            if not DEFAULT_PRIMER_LIST_PATH.exists():
                warning = f"Default primer list not found: {DEFAULT_PRIMER_LIST_PATH}. Continuing without this list."
                warnings.append(warning)
                st.warning(warning)
            else:
                default_primer_text = DEFAULT_PRIMER_LIST_PATH.read_text(encoding="utf-8")
                primer_groups.append(parse_primer_list_text(default_primer_text))
        for primer_file in params["primer_files"] or []:
            primer_text = decode_upload(primer_file)
            uploaded_primer_files.append({"name": primer_file.name, "text": primer_text})
            primer_groups.append(parse_primer_list_text(primer_text))
        existing_primers = merge_primer_inputs(primer_groups)
        masks = parse_masks(params["mask_text"], len(sequence))

        with st.spinner("Designing primers..."):
            result = design_primers(
                sequence,
                existing_primers,
                target_coverage=params["coverage"],
                masks=masks,
                max_iterations=params["max_iterations"],
                primer_name_prefix=params["prefix"],
                settings=params["settings"],
            )
    except Exception as exc:
        st.session_state["result"] = None
        st.error(f"Design failed: {exc}")
        return

    st.session_state["result"] = result
    st.session_state["sequence_length"] = len(sequence)
    st.session_state["masks"] = masks
    st.session_state["settings"] = params["settings"]
    st.session_state["runtime_seconds"] = time.perf_counter() - start
    st.session_state["run_parameters"] = {
        "coverage": params["coverage"],
        "mask_text": params["mask_text"],
        "prefix": params["prefix"],
        "delimiter": params["delimiter"],
        "max_iterations": params["max_iterations"],
        "use_default_primers": params["use_default_primers"],
    }
    st.session_state["sequence_file_name"] = params["sequence_file"].name
    st.session_state["sequence_text"] = sequence_text
    st.session_state["uploaded_primer_files"] = uploaded_primer_files
    st.session_state["default_primer_text"] = default_primer_text
    st.session_state["merged_primers"] = existing_primers
    st.session_state["warnings"] = warnings


def decode_upload(uploaded_file) -> str:
    return uploaded_file.getvalue().decode("utf-8-sig")


def render_result(result: DesignResult, sequence_length: int, masks, delimiter: str) -> None:
    settings = st.session_state.get("settings", DEFAULT_DESIGN_SETTINGS)
    render_summary(result)

    overview_tab, coverage_tab, placement_tab, primers_tab, report_tab = st.tabs(
        ["Overview", "Coverage", "Placement", "Primers", "Report"]
    )
    with overview_tab:
        render_missing_regions(result)
        render_downloads(result, delimiter, settings)
    with coverage_tab:
        st.plotly_chart(build_coverage_figure(result, sequence_length, masks), use_container_width=True)
        st.dataframe(coverage_distribution(result), hide_index=True, use_container_width=True)
    with placement_tab:
        show_read_ranges = st.checkbox("Show cover bars", value=True)
        st.plotly_chart(
            build_primer_placement_figure(result, sequence_length, masks, show_read_ranges),
            use_container_width=True,
        )
    with primers_tab:
        render_primer_table(result.primers)
    with report_tab:
        render_report(result, sequence_length, settings)


def render_summary(result: DesignResult) -> None:
    existing_count = sum(1 for primer in result.primers if primer.source == "existing")
    designed_count = sum(1 for primer in result.primers if primer.source == "designed")
    runtime = st.session_state.get("runtime_seconds", 0.0)

    cols = st.columns(5)
    cols[0].metric("Status", "Achieved" if result.achieved else "Incomplete")
    cols[1].metric("Total primers", len(result.primers))
    cols[2].metric("Existing", existing_count)
    cols[3].metric("Designed", designed_count)
    cols[4].metric("Min coverage", f"{result.min_coverage} / {result.target_coverage}")
    st.caption(f"Runtime: {runtime:.1f} seconds")


def render_missing_regions(result: DesignResult) -> None:
    if result.achieved:
        st.success("Coverage target achieved.")
        return

    st.warning("Coverage target was not achieved.")
    st.write(", ".join(format_interval(region) for region in result.missing_regions))


def render_downloads(result: DesignResult, delimiter: str, settings: DesignSettings) -> None:
    primer_text = format_primer_list(result.primers, delimiter)
    st.download_button(
        "Download primer list (.txt)",
        data=primer_text,
        file_name=DEFAULT_OUTPUT_PRIMER_LIST_PATH,
        mime="text/plain",
    )
    st.download_button(
        "Download coverage report (.txt)",
        data=coverage_report_text(result, settings=settings),
        file_name=DEFAULT_COVERAGE_REPORT_PATH,
        mime="text/plain",
    )
    st.warning(SENSITIVE_BUNDLE_NOTICE)
    st.download_button(
        "Download support bundle (.zip)",
        data=build_support_bundle(
            result=result,
            sequence_length=st.session_state["sequence_length"],
            settings=settings,
            run_parameters=st.session_state.get("run_parameters", {}),
            sequence_file_name=st.session_state.get("sequence_file_name", "plasmid_input.txt"),
            sequence_text=st.session_state.get("sequence_text", ""),
            uploaded_primer_files=st.session_state.get("uploaded_primer_files", []),
            default_primer_text=st.session_state.get("default_primer_text"),
            merged_primers=st.session_state.get("merged_primers", []),
            masks=st.session_state.get("masks", ()),
            warnings=st.session_state.get("warnings", []),
            runtime_seconds=st.session_state.get("runtime_seconds", 0.0),
            delimiter=delimiter,
        ),
        file_name=default_bundle_filename(),
        mime="application/zip",
    )


def build_coverage_figure(result: DesignResult, sequence_length: int, masks) -> go.Figure:
    positions = list(range(1, sequence_length + 1))
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=positions,
            y=list(result.coverage),
            mode="lines",
            name="Coverage",
            line=dict(color="#2563eb", width=2),
        )
    )
    fig.add_hline(
        y=result.target_coverage,
        line_dash="dash",
        line_color="#16a34a",
        annotation_text="Target",
    )
    for mask in masks:
        add_interval_shape(fig, mask, sequence_length, "#f59e0b", "Mask")
    for region in result.missing_regions:
        add_interval_shape(fig, region, sequence_length, "#ef4444", "Missing")
    fig.update_layout(
        height=430,
        margin=dict(l=20, r=20, t=20, b=20),
        xaxis_title="Position",
        yaxis_title="Coverage",
        hovermode="x unified",
    )
    return fig


def add_interval_shape(fig: go.Figure, interval, sequence_length: int, color: str, label: str) -> None:
    segments = [(interval.start, interval.end)] if interval.start <= interval.end else [
        (interval.start, sequence_length),
        (1, interval.end),
    ]
    for start, end in segments:
        fig.add_vrect(
            x0=start,
            x1=end,
            fillcolor=color,
            opacity=0.16,
            line_width=0,
            annotation_text=label,
            annotation_position="top left",
        )


def build_primer_placement_figure(
    result: DesignResult,
    sequence_length: int,
    masks,
    show_read_ranges: bool,
) -> go.Figure:
    fig = go.Figure()
    for mask in masks:
        add_interval_shape(fig, mask, sequence_length, "#f59e0b", "Mask")
    for region in result.missing_regions:
        add_interval_shape(fig, region, sequence_length, "#ef4444", "Missing")

    for primer in sorted(result.primers, key=lambda item: item.position):
        add_primer_shapes(fig, primer, sequence_length, show_read_ranges)

    for source in ("existing", "designed"):
        for direction, lane, symbol in (("Forward", 2, "triangle-right"), ("Reverse", 1, "triangle-left")):
            primers = [
                primer
                for primer in result.primers
                if primer.source == source and primer.direction == direction
            ]
            if not primers:
                continue
            fig.add_trace(
                go.Scatter(
                    x=[primer.position for primer in primers],
                    y=[lane for _ in primers],
                    mode="markers",
                    text=[primer.name for primer in primers],
                    name=f"{source} {direction}",
                    marker=dict(
                        symbol=symbol,
                        size=14,
                        color=primer_color(source),
                        line=dict(color="#111827", width=1),
                    ),
                    customdata=[
                        [
                            primer.sequence,
                            format_interval(primer.binding),
                            format_interval(primer.cover),
                        ]
                        for primer in primers
                    ],
                    hovertemplate=(
                        "<b>%{text}</b><br>"
                        "Position: %{x}<br>"
                        "Sequence: %{customdata[0]}<br>"
                        "Binding: %{customdata[1]}<br>"
                        "Cover: %{customdata[2]}<extra></extra>"
                    ),
                )
            )

    fig.update_layout(
        height=470,
        margin=dict(l=20, r=20, t=20, b=20),
        xaxis=dict(title="Position", range=[1, sequence_length]),
        yaxis=dict(
            title="",
            range=[0.35, 2.65],
            tickmode="array",
            tickvals=[1, 2],
            ticktext=["Reverse", "Forward"],
        ),
        hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
    )
    return fig


def add_primer_shapes(
    fig: go.Figure,
    primer: Primer,
    sequence_length: int,
    show_read_ranges: bool,
) -> None:
    lane = 2 if primer.direction == "Forward" else 1
    color = primer_color(primer.source)
    if show_read_ranges:
        add_noisy_connector(fig, primer, sequence_length, lane, color)
        for start, end in interval_segments(primer.cover, sequence_length):
            fig.add_shape(
                type="line",
                x0=start,
                x1=end,
                y0=lane,
                y1=lane,
                line=dict(color=color, width=7),
                opacity=0.28,
            )
    for start, end in interval_segments(primer.binding, sequence_length):
        fig.add_shape(
            type="line",
            x0=start,
            x1=end,
            y0=lane,
            y1=lane,
            line=dict(color=color, width=3),
            opacity=0.85,
        )


def add_noisy_connector(
    fig: go.Figure,
    primer: Primer,
    sequence_length: int,
    lane: int,
    color: str,
) -> None:
    if primer.direction == "Forward":
        connector = Interval(primer.position, primer.cover.start)
    else:
        connector = Interval(primer.cover.end, primer.position)

    for start, end in interval_segments(connector, sequence_length):
        fig.add_shape(
            type="line",
            x0=start,
            x1=end,
            y0=lane,
            y1=lane,
            line=dict(color=color, width=2, dash="dot"),
            opacity=0.45,
        )


def primer_color(source: str) -> str:
    return "#2563eb" if source == "existing" else "#dc2626"


def coverage_distribution(result: DesignResult) -> pd.DataFrame:
    counts = Counter(result.coverage)
    return pd.DataFrame(
        [{"Coverage": coverage, "Bases": bases} for coverage, bases in sorted(counts.items())]
    )


def render_primer_table(primers: tuple[Primer, ...]) -> None:
    source_filter = st.segmented_control("Source", options=["All", "Existing", "Designed"], default="All")
    rows = primer_rows(primers)
    if source_filter != "All":
        rows = [row for row in rows if row["Source"] == source_filter.lower()]
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)


def primer_rows(primers: tuple[Primer, ...]) -> list[dict]:
    return [
        {
            "Name": primer.name,
            "Sequence": primer.sequence,
            "Direction": primer.direction,
            "Position": primer.position,
            "Binding": format_interval(primer.binding),
            "Cover": format_interval(primer.cover),
            "Source": primer.source,
        }
        for primer in sorted(primers, key=lambda item: item.position)
    ]


def render_report(result: DesignResult, sequence_length: int, settings: DesignSettings) -> None:
    st.code(coverage_report_text(result, sequence_length, settings), language="text")


def coverage_report_text(
    result: DesignResult,
    sequence_length: int | None = None,
    settings: DesignSettings | None = None,
) -> str:
    lines = [
        f"Status: {'achieved' if result.achieved else 'incomplete'}",
        f"Target coverage: {result.target_coverage}",
        f"Minimum coverage: {result.min_coverage}",
        f"Total primers: {len(result.primers)}",
        f"Existing primers: {sum(1 for primer in result.primers if primer.source == 'existing')}",
        f"Designed primers: {sum(1 for primer in result.primers if primer.source == 'designed')}",
    ]
    if sequence_length is not None:
        lines.append(f"Sequence length: {sequence_length}")
    if settings is not None:
        lines.extend(
            [
                "Design settings:",
                f"  Read length: {settings.read_length}",
                f"  Noisy bases: {settings.noise_length}",
                f"  Primer length: {settings.min_primer_len}-{settings.max_primer_len}",
                f"  Tm: {settings.min_tm}-{settings.max_tm}",
                f"  GC: {settings.min_gc}-{settings.max_gc}",
                f"  Minimum binding gap: {settings.min_binding_gap}",
                f"  Preferred pair distance: {settings.preferred_pair_min}-{settings.preferred_pair_max}",
            ]
        )
    if result.missing_regions:
        lines.append("Missing regions: " + ", ".join(format_interval(region) for region in result.missing_regions))
    lines.append("Coverage distribution:")
    for row in coverage_distribution(result).to_dict("records"):
        lines.append(f"  {row['Coverage']}: {row['Bases']} bp")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    main()
