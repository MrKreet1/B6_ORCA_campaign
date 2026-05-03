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


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def to_float(text: str) -> float:
    return float(text)


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

    return save_figure(fig, output_dir, "B6_vibrational_frequencies_bar")


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

    return save_figure(fig, output_dir, "B6_normal_mode_amplitudes_heatmap")


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

    return save_figure(fig, output_dir, "B6_max_amplitude_by_mode")


def write_manifest(output_dir: Path, files: Sequence[Tuple[Path, str]]) -> None:
    manifest = output_dir / "plots_manifest.csv"
    with manifest.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "description"])
        for path, description in files:
            writer.writerow([path.name, description])


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot B6 vibrational analysis with Matplotlib.")
    parser.add_argument("--input-dir", default="results/vibrations/B6")
    parser.add_argument("--output-dir", default="results/vibrations/B6/matplotlib_plots")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    freq_rows = read_csv_rows(input_dir / "B6_all_vibrational_frequencies.csv")
    mode_rows = read_csv_rows(input_dir / "B6_mode_summary.csv")
    amplitude_rows = read_csv_rows(input_dir / "B6_normal_mode_amplitudes.csv")

    manifest_items: List[Tuple[Path, str]] = []
    for path in plot_frequency_bars(freq_rows, output_dir):
        manifest_items.append((path, "Bar plot of all non-zero B6 vibrational frequencies."))
    for path in plot_spectrum_lines(freq_rows, output_dir):
        manifest_items.append((path, "Line-spectrum representation of B6 vibrational frequencies."))
    for path in plot_amplitude_heatmap(amplitude_rows, output_dir):
        manifest_items.append((path, "Heatmap of normalized atom displacement amplitudes by mode."))
    for path in plot_max_amplitude(mode_rows, output_dir):
        manifest_items.append((path, "Maximum normal-mode displacement amplitude and dominant atom by mode."))

    write_manifest(output_dir, manifest_items)
    print(f"Wrote Matplotlib plots to: {output_dir.resolve()}")
    print(f"Wrote manifest: {(output_dir / 'plots_manifest.csv').resolve()}")


if __name__ == "__main__":
    main()
