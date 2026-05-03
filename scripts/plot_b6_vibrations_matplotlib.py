#!/usr/bin/env python3
"""Create Matplotlib plots for the selected B6 vibrational analysis."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

FIGURE_ORDER = [
    (
        "Figure_10_B6_vibrational_frequencies",
        "Figure 10. B6 Vibrational Frequencies",
        "Bar chart of all 12 non-zero normal-mode frequencies from B6_all_vibrational_frequencies.csv.",
    ),
    (
        "Figure_11_B6_max_amplitude_by_mode",
        "Figure 11. Maximum Amplitude By Mode",
        "Maximum normalized normal-mode displacement amplitude and dominant atom label for each mode.",
    ),
    (
        "Figure_12_B6_frequency_vs_amplitude",
        "Figure 12. Frequency Vs Maximum Amplitude",
        "Scatter plot for checking whether higher-frequency modes show systematically different displacement amplitudes.",
    ),
    (
        "Figure_13_B6_atom_participation_heatmap",
        "Figure 13. Atom Participation Heatmap",
        "Heatmap of atom-by-mode displacement amplitudes from B6_normal_mode_amplitudes.csv.",
    ),
    (
        "Figure_14_B6_dominant_atom_by_mode",
        "Figure 14. Dominant Atom By Mode",
        "Dominant atom index for every normal mode, based on the largest normalized displacement amplitude.",
    ),
    (
        "Figure_15_B6_frequency_distribution",
        "Figure 15. Frequency Distribution",
        "Histogram showing how B6 vibrational modes are distributed across low, medium, and high frequency ranges.",
    ),
    (
        "Figure_16_final_relative_energies_labeled",
        "Figure 16. Final Relative Energies",
        "Final PBE0-D4/def2-TZVP candidate energies by rank, with multiplicity labels and the best_B6 marker.",
    ),
    (
        "Figure_17_screening_success_rate",
        "Figure 17. Screening Success Rate",
        "Successful versus failed/not-normal R2SCAN-3C screening calculations.",
    ),
    (
        "B6_vibrational_spectrum_lines",
        "Supplement. Frequency Line Spectrum",
        "Line-spectrum representation of the same non-zero vibrational frequencies.",
    ),
]


def html_escape(text: object) -> str:
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def to_float(text: str) -> float:
    return float(text)


def truthy(text: str) -> bool:
    return str(text).strip().lower() == "true"


def save_figure(fig: plt.Figure, output_dir: Path, stem: str) -> List[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = [output_dir / f"{stem}.png", output_dir / f"{stem}.svg"]
    fig.savefig(paths[0], dpi=220, bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight")
    plt.close(fig)
    return paths


def plot_frequency_bars(freq_rows: Sequence[Dict[str, str]], output_dir: Path) -> List[Path]:
    modes = [int(row["mode_number"]) for row in freq_rows]
    freqs = [to_float(row["frequency_cm-1"]) for row in freq_rows]

    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    bars = ax.bar(modes, freqs, color="#4f8a8b", edgecolor="#1f4142", linewidth=0.8)
    ax.set_title("B6 vibrational frequencies")
    ax.set_xlabel("Normal mode number")
    ax.set_ylabel("Frequency, cm$^{-1}$")
    ax.set_xticks(modes)
    ax.grid(axis="y", color="#d9dddd", linewidth=0.8)
    ax.set_axisbelow(True)

    for bar, freq in zip(bars, freqs):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(freqs) * 0.018,
            f"{freq:.0f}",
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=90,
        )

    return save_figure(fig, output_dir, "Figure_10_B6_vibrational_frequencies")


def plot_spectrum_lines(freq_rows: Sequence[Dict[str, str]], output_dir: Path) -> List[Path]:
    freqs = [to_float(row["frequency_cm-1"]) for row in freq_rows]
    modes = [int(row["mode_number"]) for row in freq_rows]

    fig, ax = plt.subplots(figsize=(10.5, 3.8))
    for mode, freq in zip(modes, freqs):
        ax.vlines(freq, 0.0, 1.0, color="#4f8a8b", linewidth=2.2)
        ax.text(freq, 1.04, f"{mode}", ha="center", va="bottom", fontsize=8)

    ax.set_title("B6 vibrational spectrum, line representation")
    ax.set_xlabel("Frequency, cm$^{-1}$")
    ax.set_yticks([])
    ax.set_ylim(0, 1.18)
    ax.set_xlim(0, max(freqs) * 1.08)
    ax.grid(axis="x", color="#d9dddd", linewidth=0.8)
    ax.set_axisbelow(True)

    return save_figure(fig, output_dir, "B6_vibrational_spectrum_lines")


def amplitude_matrix(amplitude_rows: Sequence[Dict[str, str]]) -> Tuple[List[int], List[int], List[List[float]]]:
    modes = sorted({int(row["mode_number"]) for row in amplitude_rows})
    atoms = sorted({int(row["atom_index"]) for row in amplitude_rows})
    values = {(int(row["atom_index"]), int(row["mode_number"])): to_float(row["amplitude"]) for row in amplitude_rows}
    matrix = [[values.get((atom, mode), 0.0) for mode in modes] for atom in atoms]
    return atoms, modes, matrix


def plot_amplitude_heatmap(amplitude_rows: Sequence[Dict[str, str]], output_dir: Path) -> List[Path]:
    atoms, modes, matrix = amplitude_matrix(amplitude_rows)

    fig, ax = plt.subplots(figsize=(10.8, 4.8))
    image = ax.imshow(matrix, aspect="auto", cmap="viridis")
    ax.set_title("B6 normal-mode displacement amplitudes")
    ax.set_xlabel("Normal mode number")
    ax.set_ylabel("Atom")
    ax.set_xticks(range(len(modes)), [str(mode) for mode in modes])
    ax.set_yticks(range(len(atoms)), [f"B{atom}" for atom in atoms])

    for y, atom_values in enumerate(matrix):
        for x, value in enumerate(atom_values):
            ax.text(x, y, f"{value:.2f}", ha="center", va="center", fontsize=7, color="white" if value > 0.36 else "#1f2727")

    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label("Normalized amplitude")

    return save_figure(fig, output_dir, "Figure_13_B6_atom_participation_heatmap")


def plot_max_amplitude(mode_rows: Sequence[Dict[str, str]], output_dir: Path) -> List[Path]:
    modes = [int(row["mode_number"]) for row in mode_rows]
    freqs = [to_float(row["frequency_cm-1"]) for row in mode_rows]
    max_amp = [to_float(row["max_amplitude"]) for row in mode_rows]
    dominant_atoms = [row["dominant_atom_index"] for row in mode_rows]

    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    ax.plot(modes, max_amp, marker="o", color="#4f8a8b", linewidth=2)
    ax.set_title("Maximum atom displacement amplitude by normal mode")
    ax.set_xlabel("Normal mode number")
    ax.set_ylabel("Maximum normalized amplitude")
    ax.set_xticks(modes)
    ax.grid(color="#d9dddd", linewidth=0.8)
    ax.set_axisbelow(True)

    for mode, freq, amp, atom in zip(modes, freqs, max_amp, dominant_atoms):
        ax.text(mode, amp + 0.012, f"B{atom}\n{freq:.0f}", ha="center", va="bottom", fontsize=8)

    return save_figure(fig, output_dir, "Figure_11_B6_max_amplitude_by_mode")


def plot_frequency_vs_amplitude(mode_rows: Sequence[Dict[str, str]], output_dir: Path) -> List[Path]:
    freqs = [to_float(row["frequency_cm-1"]) for row in mode_rows]
    max_amp = [to_float(row["max_amplitude"]) for row in mode_rows]
    modes = [int(row["mode_number"]) for row in mode_rows]

    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    ax.scatter(freqs, max_amp, s=70, color="#4f8a8b", edgecolor="#1f4142", linewidth=0.8)
    ax.set_title("B6 frequency vs maximum normal-mode amplitude")
    ax.set_xlabel("Frequency, cm$^{-1}$")
    ax.set_ylabel("Maximum normalized amplitude")
    ax.grid(color="#d9dddd", linewidth=0.8)
    ax.set_axisbelow(True)

    for mode, freq, amp in zip(modes, freqs, max_amp):
        ax.text(freq, amp + 0.008, str(mode), ha="center", va="bottom", fontsize=8)

    return save_figure(fig, output_dir, "Figure_12_B6_frequency_vs_amplitude")


def plot_dominant_atom_by_mode(mode_rows: Sequence[Dict[str, str]], output_dir: Path) -> List[Path]:
    modes = [int(row["mode_number"]) for row in mode_rows]
    dominant_atoms = [int(row["dominant_atom_index"]) for row in mode_rows]
    freqs = [to_float(row["frequency_cm-1"]) for row in mode_rows]

    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    ax.step(modes, dominant_atoms, where="mid", color="#4f8a8b", linewidth=1.6)
    ax.scatter(modes, dominant_atoms, s=70, color="#4f8a8b", edgecolor="#1f4142", linewidth=0.8)
    ax.set_title("Dominant atom in each B6 normal mode")
    ax.set_xlabel("Normal mode number")
    ax.set_ylabel("Dominant atom index")
    ax.set_xticks(modes)
    ax.set_yticks(sorted(set(dominant_atoms)))
    ax.grid(color="#d9dddd", linewidth=0.8)
    ax.set_axisbelow(True)

    for mode, atom, freq in zip(modes, dominant_atoms, freqs):
        ax.text(mode, atom + 0.08, f"{freq:.0f}", ha="center", va="bottom", fontsize=8)

    return save_figure(fig, output_dir, "Figure_14_B6_dominant_atom_by_mode")


def plot_frequency_distribution(freq_rows: Sequence[Dict[str, str]], output_dir: Path) -> List[Path]:
    freqs = [to_float(row["frequency_cm-1"]) for row in freq_rows]
    bins = [0, 400, 800, 1200, 1600]
    labels = ["0-400", "400-800", "800-1200", "1200-1600"]

    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    counts, _, patches = ax.hist(freqs, bins=bins, color="#4f8a8b", edgecolor="#1f4142", linewidth=0.8)
    ax.set_title("Distribution of B6 vibrational frequencies")
    ax.set_xlabel("Frequency range, cm$^{-1}$")
    ax.set_ylabel("Number of modes")
    ax.set_xticks([(bins[i] + bins[i + 1]) / 2 for i in range(len(labels))], labels)
    ax.grid(axis="y", color="#d9dddd", linewidth=0.8)
    ax.set_axisbelow(True)

    for patch, count in zip(patches, counts):
        ax.text(patch.get_x() + patch.get_width() / 2, count + 0.06, f"{int(count)}", ha="center", va="bottom", fontsize=9)

    return save_figure(fig, output_dir, "Figure_15_B6_frequency_distribution")


def plot_final_relative_energies_labeled(final_rows: Sequence[Dict[str, str]], output_dir: Path) -> List[Path]:
    rows = [row for row in final_rows if row.get("relative_energy_ev")]
    ranks = list(range(1, len(rows) + 1))
    energies = [to_float(row["relative_energy_ev"]) for row in rows]
    multiplicities = [row.get("multiplicity", "") for row in rows]

    fig, ax = plt.subplots(figsize=(10.0, 4.8))
    ax.plot(ranks, energies, marker="o", color="#4f8a8b", linewidth=1.8)
    ax.scatter([1], [energies[0]], s=110, color="#b3432f", edgecolor="#6e2419", linewidth=1.0, zorder=4, label="best_B6")
    ax.set_title("Final candidates: relative energies with best_B6 marker")
    ax.set_xlabel("Final candidate rank")
    ax.set_ylabel("Relative energy, eV")
    ax.set_xticks(ranks)
    ax.grid(color="#d9dddd", linewidth=0.8)
    ax.set_axisbelow(True)

    for rank, energy, mult in zip(ranks, energies, multiplicities):
        ax.text(rank, energy + max(energies + [1e-6]) * 0.04, f"m={mult}", ha="center", va="bottom", fontsize=8)
    ax.legend(loc="upper left")

    return save_figure(fig, output_dir, "Figure_16_final_relative_energies_labeled")


def plot_screening_success_rate(screening_rows: Sequence[Dict[str, str]], output_dir: Path) -> List[Path]:
    total = len(screening_rows)
    successful = sum(1 for row in screening_rows if truthy(row.get("normal_termination", "")))
    failed = total - successful

    fig, (ax_bar, ax_pie) = plt.subplots(1, 2, figsize=(10.0, 4.6), gridspec_kw={"width_ratios": [1.05, 1.0]})
    labels = ["successful", "failed/not normal"]
    values = [successful, failed]
    colors = ["#4f8a8b", "#b75b4a"]

    bars = ax_bar.bar(labels, values, color=colors, edgecolor="#333", linewidth=0.8)
    ax_bar.set_title("Screening calculation status")
    ax_bar.set_ylabel("Number of calculations")
    ax_bar.grid(axis="y", color="#d9dddd", linewidth=0.8)
    ax_bar.set_axisbelow(True)
    for bar, value in zip(bars, values):
        ax_bar.text(bar.get_x() + bar.get_width() / 2, value + max(values) * 0.025, str(value), ha="center", va="bottom", fontsize=10)

    ax_pie.pie(values, labels=labels, colors=colors, autopct=lambda pct: f"{pct:.1f}%", startangle=90)
    ax_pie.set_title(f"Total screening .out = {total}")

    return save_figure(fig, output_dir, "Figure_17_screening_success_rate")


def write_manifest(output_dir: Path, files: Sequence[Tuple[Path, str]]) -> None:
    manifest = output_dir / "plots_manifest.csv"
    with manifest.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "description"])
        for path, description in files:
            writer.writerow([path.name, description])


def write_gallery_page(output_dir: Path) -> Path:
    page = output_dir / "index.html"
    cards: List[str] = []
    for stem, title, description in FIGURE_ORDER:
        svg = f"{stem}.svg"
        png = f"{stem}.png"
        if not (output_dir / svg).exists():
            continue
        cards.append(
            f"""
      <article class="plot">
        <header>
          <h2>{html_escape(title)}</h2>
          <p>{html_escape(description)}</p>
          <div class="links">
            <a href="{html_escape(svg)}">SVG</a>
            <a href="{html_escape(png)}">PNG</a>
          </div>
        </header>
        <a class="image-link" href="{html_escape(svg)}">
          <img src="{html_escape(svg)}" alt="{html_escape(title)}">
        </a>
      </article>"""
        )

    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>B6 Matplotlib Graph Gallery</title>
  <style>
    :root {{
      color-scheme: light;
      font-family: Arial, sans-serif;
      background: #f5f7f7;
      color: #1f2727;
    }}
    body {{
      margin: 0;
      background: #f5f7f7;
    }}
    .topbar {{
      background: #ffffff;
      border-bottom: 1px solid #d8dfdf;
      padding: 18px 24px 16px;
      position: sticky;
      top: 0;
      z-index: 2;
    }}
    h1 {{
      margin: 0 0 8px;
      font-size: 24px;
      line-height: 1.2;
    }}
    .summary {{
      margin: 0;
      color: #586464;
      font-size: 14px;
      line-height: 1.45;
    }}
    .nav {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-top: 12px;
    }}
    a {{
      color: #23676a;
      text-decoration: none;
    }}
    a:hover {{
      text-decoration: underline;
    }}
    .nav a, .links a {{
      border: 1px solid #b9c9c9;
      border-radius: 6px;
      padding: 6px 9px;
      background: #f9fbfb;
      font-size: 13px;
    }}
    main {{
      max-width: 1260px;
      margin: 0 auto;
      padding: 18px;
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(430px, 1fr));
      gap: 16px;
    }}
    .plot {{
      background: #ffffff;
      border: 1px solid #d7dfdf;
      border-radius: 8px;
      overflow: hidden;
    }}
    .plot header {{
      padding: 13px 14px 10px;
      border-bottom: 1px solid #e4e9e9;
    }}
    .plot h2 {{
      margin: 0 0 6px;
      font-size: 16px;
      line-height: 1.25;
    }}
    .plot p {{
      margin: 0;
      color: #5c6666;
      font-size: 13px;
      line-height: 1.45;
    }}
    .links {{
      display: flex;
      gap: 8px;
      margin-top: 10px;
    }}
    .image-link {{
      display: block;
      background: #ffffff;
    }}
    img {{
      width: 100%;
      height: auto;
      display: block;
    }}
    @media (max-width: 520px) {{
      main {{
        grid-template-columns: 1fr;
        padding: 10px;
      }}
      .topbar {{
        padding: 14px;
      }}
    }}
  </style>
</head>
<body>
  <section class="topbar">
    <h1>B6 Matplotlib Graph Gallery</h1>
    <p class="summary">Vibrational, normal-mode, final-energy, and screening-status plots generated from the B6 ORCA campaign CSV files.</p>
    <nav class="nav">
      <a href="../README.md">Vibration README</a>
      <a href="../../../../README.md">Main README</a>
      <a href="plots_manifest.csv">Plot Manifest</a>
    </nav>
  </section>
  <main>
{''.join(cards)}
  </main>
</body>
</html>
"""
    page.write_text(content, encoding="utf-8")
    return page


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot B6 vibrational analysis with Matplotlib.")
    parser.add_argument("--project-dir", default=".")
    parser.add_argument("--input-dir", default="results/vibrations/B6")
    parser.add_argument("--output-dir", default="results/vibrations/B6/matplotlib_plots")
    parser.add_argument("--final-csv", default="results/final_results.csv")
    parser.add_argument("--screening-csv", default="results/results.csv")
    args = parser.parse_args()

    project_dir = Path(args.project_dir).resolve()
    input_dir = Path(args.input_dir)
    if not input_dir.is_absolute():
        input_dir = project_dir / input_dir
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = project_dir / output_dir
    final_csv = Path(args.final_csv)
    if not final_csv.is_absolute():
        final_csv = project_dir / final_csv
    screening_csv = Path(args.screening_csv)
    if not screening_csv.is_absolute():
        screening_csv = project_dir / screening_csv

    freq_rows = read_csv_rows(input_dir / "B6_all_vibrational_frequencies.csv")
    mode_rows = read_csv_rows(input_dir / "B6_mode_summary.csv")
    amplitude_rows = read_csv_rows(input_dir / "B6_normal_mode_amplitudes.csv")
    final_rows = read_csv_rows(final_csv)
    screening_rows = read_csv_rows(screening_csv)

    manifest_items: List[Tuple[Path, str]] = []
    for path in plot_frequency_bars(freq_rows, output_dir):
        manifest_items.append((path, "Figure 10. Bar chart of all 12 non-zero B6 vibrational frequencies."))
    for path in plot_max_amplitude(mode_rows, output_dir):
        manifest_items.append((path, "Figure 11. Maximum normal-mode displacement amplitude by mode."))
    for path in plot_frequency_vs_amplitude(mode_rows, output_dir):
        manifest_items.append((path, "Figure 12. Scatter plot of frequency versus maximum normal-mode amplitude."))
    for path in plot_amplitude_heatmap(amplitude_rows, output_dir):
        manifest_items.append((path, "Figure 13. Atom participation heatmap by normal mode."))
    for path in plot_dominant_atom_by_mode(mode_rows, output_dir):
        manifest_items.append((path, "Figure 14. Dominant atom index by normal mode."))
    for path in plot_frequency_distribution(freq_rows, output_dir):
        manifest_items.append((path, "Figure 15. Histogram of B6 vibrational frequency distribution."))
    for path in plot_final_relative_energies_labeled(final_rows, output_dir):
        manifest_items.append((path, "Figure 16. Final relative energies with best_B6 marker and multiplicity labels."))
    for path in plot_screening_success_rate(screening_rows, output_dir):
        manifest_items.append((path, "Figure 17. Screening success/fail summary."))

    # Extra line-spectrum view kept as a supplementary Matplotlib plot.
    for path in plot_spectrum_lines(freq_rows, output_dir):
        manifest_items.append((path, "Supplementary line-spectrum representation of B6 vibrational frequencies."))

    write_manifest(output_dir, manifest_items)
    gallery_page = write_gallery_page(output_dir)
    print(f"Wrote Matplotlib plots to: {output_dir.resolve()}")
    print(f"Wrote manifest: {(output_dir / 'plots_manifest.csv').resolve()}")
    print(f"Wrote gallery page: {gallery_page.resolve()}")


if __name__ == "__main__":
    main()
