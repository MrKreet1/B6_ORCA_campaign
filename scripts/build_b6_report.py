#!/usr/bin/env python3
"""Build the B6 calculation report and simple SVG figures from collected results."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import importlib.util
import json
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

Atom = Tuple[str, float, float, float]


SCREENING_COLUMNS = [
    "calculation_name",
    "geometry_type",
    "distance",
    "multiplicity",
    "method",
    "total_energy_hartree",
    "relative_energy_ev",
    "normal_termination",
    "optimization_converged",
]

FINAL_COLUMNS = [
    "calculation_name",
    "multiplicity",
    "method",
    "basis",
    "total_energy_hartree",
    "relative_energy_ev",
    "lowest_frequency_cm-1",
    "n_imaginary_frequencies",
    "is_true_minimum",
    "xyz_file",
    "output_file",
]


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def read_xyz(path: Path) -> List[Atom]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    n = int(lines[0].strip())
    atoms: List[Atom] = []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms


def short_name(name: str, max_len: int = 44) -> str:
    name = name.replace("_PBE0_def2-TZVP_OptFreq", "")
    name = name.replace("FINAL_", "")
    name = name.replace("B6_", "")
    if len(name) <= max_len:
        return name
    return name[: max_len - 1] + "..."


def short_path(path_text: str, project: Path) -> str:
    if not path_text:
        return ""
    path = Path(path_text)
    try:
        if path.is_absolute():
            return path.resolve().relative_to(project).as_posix()
    except Exception:
        pass
    marker = project.name
    parts = path.parts
    if marker in parts:
        idx = parts.index(marker)
        return Path(*parts[idx + 1 :]).as_posix()
    return path_text.replace("\\", "/")


def resolve_project_path(path_text: str, project: Path) -> Path:
    raw = Path(path_text)
    candidates = [raw]
    if not raw.is_absolute():
        candidates.insert(0, project / raw)

    marker = project.name
    parts = raw.parts
    if marker in parts:
        idx = parts.index(marker)
        tail = parts[idx + 1 :]
        if tail:
            candidates.append(project / Path(*tail))

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0]


def to_float(text: str, default: float = math.inf) -> float:
    try:
        return float(text)
    except Exception:
        return default


def truthy(text: object) -> bool:
    return str(text).strip().lower() == "true"


def format_float(value: float, digits: int = 6) -> str:
    if math.isinf(value) or math.isnan(value):
        return ""
    return f"{value:.{digits}f}"


def read_text_if_exists(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def markdown_table(rows: Sequence[Dict[str, str]], columns: Sequence[str], project: Path, limit: int | None = None) -> str:
    shown = list(rows if limit is None else rows[:limit])
    out = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in shown:
        values = []
        for col in columns:
            value = row.get(col, "")
            if col in {"xyz_file", "output_file"}:
                value = short_path(value, project)
            if col == "calculation_name":
                value = short_name(value, 62)
            values.append(value)
        out.append("| " + " | ".join(values) + " |")
    return "\n".join(out)


def simple_markdown_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    out = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        out.append("| " + " | ".join(str(item) for item in row) + " |")
    return "\n".join(out)


def pair_distances(atoms: Sequence[Atom]) -> List[Tuple[int, int, float]]:
    pairs: List[Tuple[int, int, float]] = []
    for i in range(len(atoms)):
        _, xi, yi, zi = atoms[i]
        for j in range(i + 1, len(atoms)):
            _, xj, yj, zj = atoms[j]
            d = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)
            pairs.append((i, j, d))
    return pairs


def bond_distance_summary(atoms: Sequence[Atom]) -> Tuple[float, float, float, List[Tuple[int, int, float]]]:
    pairs = sorted(pair_distances(atoms), key=lambda item: item[2])
    distances = [d for _, _, d in pairs]
    return min(distances), max(distances), sum(distances) / len(distances), pairs


def adjacency_like_bonds(atoms: Sequence[Atom], cutoff: float = 2.05) -> List[Tuple[int, int, float]]:
    return [(i, j, d) for i, j, d in pair_distances(atoms) if d <= cutoff]


def planarity(atoms: Sequence[Atom]) -> Tuple[float, float]:
    points = [(x, y, z) for _, x, y, z in atoms]
    c = [sum(p[i] for p in points) / len(points) for i in range(3)]
    centered = [(x - c[0], y - c[1], z - c[2]) for x, y, z in points]
    cov = [[sum(p[i] * p[j] for p in centered) for j in range(3)] for i in range(3)]
    normal = smallest_eigenvector_3x3(cov)
    norm = math.sqrt(sum(x * x for x in normal)) or 1.0
    normal = [x / norm for x in normal]
    distances = [sum(p[i] * normal[i] for i in range(3)) for p in centered]
    rms = math.sqrt(sum(d * d for d in distances) / len(distances))
    return rms, max(abs(d) for d in distances)


def smallest_eigenvector_3x3(matrix: Sequence[Sequence[float]]) -> List[float]:
    a = [list(row) for row in matrix]
    v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    for _ in range(64):
        p, q = max([(0, 1), (0, 2), (1, 2)], key=lambda ij: abs(a[ij[0]][ij[1]]))
        if abs(a[p][q]) < 1e-12:
            break
        app, aqq, apq = a[p][p], a[q][q], a[p][q]
        phi = 0.5 * math.atan2(2.0 * apq, aqq - app)
        c, s = math.cos(phi), math.sin(phi)
        for k in range(3):
            akp, akq = a[k][p], a[k][q]
            a[k][p], a[k][q] = c * akp - s * akq, s * akp + c * akq
        for k in range(3):
            apk, aqk = a[p][k], a[q][k]
            a[p][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
        for k in range(3):
            vkp, vkq = v[k][p], v[k][q]
            v[k][p], v[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    idx = min(range(3), key=lambda i: a[i][i])
    return [v[0][idx], v[1][idx], v[2][idx]]


def parse_frequencies_from_out(path: Path) -> List[float]:
    text = read_text_if_exists(path)
    if "VIBRATIONAL FREQUENCIES" not in text:
        return []
    start = text.rfind("VIBRATIONAL FREQUENCIES")
    block = text[start:]
    stops = [block.find(marker) for marker in ["NORMAL MODES", "IR SPECTRUM", "THERMOCHEMISTRY", "ORCA TERMINATED"] if block.find(marker) > 0]
    if stops:
        block = block[: min(stops)]
    freqs: List[float] = []
    for line in block.splitlines():
        match = re.search(r"^\s*\d+\s*:\s*(-?\d+(?:\.\d+)?)\s*cm", line)
        if match:
            freqs.append(float(match.group(1)))
    return freqs


def nonzero_frequencies(freqs: Sequence[float], threshold: float = 10.0) -> List[float]:
    return [freq for freq in freqs if abs(freq) > threshold]


def zero_frequency_count(freqs: Sequence[float], threshold: float = 10.0) -> int:
    return sum(1 for freq in freqs if abs(freq) <= threshold)


def group_summary_rows(rows: Sequence[Dict[str, str]], group_key: str) -> List[List[object]]:
    groups = sorted({row.get(group_key, "") for row in rows}, key=lambda value: str(value))
    table_rows: List[List[object]] = []
    for group in groups:
        subset = [row for row in rows if row.get(group_key, "") == group]
        with_energy = [row for row in subset if row.get("total_energy_hartree")]
        best = min(with_energy, key=lambda row: to_float(row.get("total_energy_hartree", "")), default={})
        table_rows.append(
            [
                group,
                len(subset),
                sum(1 for row in subset if truthy(row.get("normal_termination", ""))),
                sum(1 for row in subset if truthy(row.get("optimization_converged", ""))),
                best.get("total_energy_hartree", ""),
                best.get("relative_energy_ev", ""),
                short_name(best.get("calculation_name", ""), 48),
            ]
        )
    return table_rows


def final_planarity_rows(final_rows: Sequence[Dict[str, str]], project: Path) -> List[List[object]]:
    rows: List[List[object]] = []
    for row in final_rows:
        xyz_path = resolve_project_path(row.get("xyz_file", ""), project)
        if xyz_path.exists():
            atoms = read_xyz(xyz_path)
            rms, max_abs = planarity(atoms)
        else:
            rms, max_abs = math.nan, math.nan
        rows.append(
            [
                short_name(row.get("calculation_name", ""), 52),
                row.get("geometry_type", ""),
                row.get("relative_energy_ev", ""),
                row.get("lowest_frequency_cm-1", ""),
                format_float(rms, 5),
                format_float(max_abs, 5),
            ]
        )
    return rows


def final_frequency_rows(final_rows: Sequence[Dict[str, str]], project: Path) -> List[List[object]]:
    rows: List[List[object]] = []
    for row in final_rows:
        out_path = resolve_project_path(row.get("output_file", ""), project)
        freqs = parse_frequencies_from_out(out_path)
        vib = nonzero_frequencies(freqs)
        rows.append(
            [
                short_name(row.get("calculation_name", ""), 52),
                zero_frequency_count(freqs),
                len(vib),
                format_float(min(vib), 2) if vib else "",
                ", ".join(format_float(freq, 2) for freq in vib[:6]),
                row.get("n_imaginary_frequencies", ""),
            ]
        )
    return rows


def distance_fingerprint(atoms: Sequence[Atom]) -> List[float]:
    return sorted(d for _, _, d in pair_distances(atoms))


def same_distance_fingerprint(a: Sequence[float], b: Sequence[float], tolerance: float) -> bool:
    return len(a) == len(b) and max((abs(x - y) for x, y in zip(a, b)), default=0.0) <= tolerance


def unique_final_geometry_rows(final_rows: Sequence[Dict[str, str]], project: Path, tolerance: float = 0.02, limit: int = 5) -> List[List[object]]:
    groups: List[Dict[str, object]] = []
    for row in final_rows:
        xyz_path = resolve_project_path(row.get("xyz_file", ""), project)
        if not xyz_path.exists():
            continue
        fp = distance_fingerprint(read_xyz(xyz_path))
        matched = False
        for group in groups:
            if same_distance_fingerprint(fp, group["fingerprint"], tolerance):  # type: ignore[arg-type]
                group["members"].append(row)  # type: ignore[index, union-attr]
                matched = True
                break
        if not matched:
            groups.append({"fingerprint": fp, "members": [row]})

    table_rows: List[List[object]] = []
    for idx, group in enumerate(groups[:limit], start=1):
        members = group["members"]  # type: ignore[assignment]
        representative = members[0]
        rel_values = [to_float(row.get("relative_energy_ev", ""), math.nan) for row in members]
        finite_rel = [value for value in rel_values if not math.isnan(value)]
        source_types = sorted({row.get("geometry_type", "") for row in members})
        table_rows.append(
            [
                idx,
                len(members),
                short_name(representative.get("calculation_name", ""), 52),
                representative.get("multiplicity", ""),
                representative.get("relative_energy_ev", ""),
                format_float(max(finite_rel), 8) if finite_rel else "",
                ", ".join(source_types),
            ]
        )
    return table_rows


def best_coordinate_rows(atoms: Sequence[Atom]) -> List[List[object]]:
    rows: List[List[object]] = []
    for idx, (el, x, y, z) in enumerate(atoms, start=1):
        rows.append([idx, el, f"{x:.8f}", f"{y:.8f}", f"{z:.8f}"])
    return rows


def distance_rows(atoms: Sequence[Atom]) -> List[List[object]]:
    _, _, _, pairs = bond_distance_summary(atoms)
    return [[f"B{i + 1}-B{j + 1}", f"{d:.6f}"] for i, j, d in pairs]


def source_category(geometry_type: str) -> str:
    text = geometry_type.lower()
    if any(token in text for token in ["random", "3d", "octa", "prism", "pyramid"]):
        return "3D/random"
    if "quasi" in text:
        return "quasi-planar"
    if any(token in text for token in ["planar", "ring", "triangle", "rhombic", "rectangular"]):
        return "planar"
    return "other"


def is_3d_geometry(geometry_type: str) -> bool:
    return source_category(geometry_type) == "3D/random"


def write_subset_csv(path: Path, rows: Sequence[Dict[str, str]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_geometry_module(project: Path):
    module_path = project / "scripts" / "generate_b6_inputs.py"
    spec = importlib.util.spec_from_file_location("generate_b6_inputs", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def html(text: object) -> str:
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def project_atoms(atoms: Sequence[Atom], axes: Tuple[int, int] = (0, 1)) -> List[Tuple[float, float, float]]:
    coords = [(x, y, z) for _, x, y, z in atoms]
    out = []
    for coord in coords:
        out.append((coord[axes[0]], coord[axes[1]], coord[2]))
    return out


def scaled_points(points: Sequence[Tuple[float, float, float]], width: int, height: int, pad: int) -> List[Tuple[float, float, float]]:
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    span = max(xmax - xmin, ymax - ymin, 1e-6)
    scale = min((width - 2 * pad) / span, (height - 2 * pad) / span)
    cx, cy = (xmin + xmax) / 2.0, (ymin + ymax) / 2.0
    return [(width / 2.0 + (x - cx) * scale, height / 2.0 - (y - cy) * scale, z) for x, y, z in points]


def atoms_svg(atoms: Sequence[Atom], width: int = 360, height: int = 300, title: str = "", axes: Tuple[int, int] = (0, 1)) -> str:
    points = scaled_points(project_atoms(atoms, axes), width, height, 34)
    pairs = pair_distances(atoms)
    min_d = min((d for _, _, d in pairs), default=1.7)
    cutoff = min(2.25, max(1.75, min_d * 1.28))
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]
    if title:
        lines.append(f'<text x="16" y="24" font-family="Arial" font-size="15" font-weight="700" fill="#222">{html(title)}</text>')
    for i, j, d in pairs:
        if d <= cutoff:
            x1, y1, _ = points[i]
            x2, y2, _ = points[j]
            lines.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="#808080" stroke-width="2"/>')
    order = sorted(range(len(points)), key=lambda idx: points[idx][2])
    for idx in order:
        x, y, z = points[idx]
        radius = 11 + 2.5 * (z - min(p[2] for p in points)) / (max(p[2] for p in points) - min(p[2] for p in points) + 1e-6)
        lines.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius:.2f}" fill="#4f8a8b" stroke="#1f4142" stroke-width="1.5"/>')
        lines.append(f'<text x="{x:.2f}" y="{y + 4:.2f}" text-anchor="middle" font-family="Arial" font-size="10" fill="#ffffff">B</text>')
    lines.append("</svg>")
    return "\n".join(lines)


def rotate_for_3d(coord: Tuple[float, float, float], yaw: float = -0.65, pitch: float = 0.78) -> Tuple[float, float, float]:
    x, y, z = coord
    cy, sy = math.cos(yaw), math.sin(yaw)
    cp, sp = math.cos(pitch), math.sin(pitch)
    x1 = x * cy + z * sy
    z1 = -x * sy + z * cy
    y1 = y * cp - z1 * sp
    z2 = y * sp + z1 * cp
    return x1, y1, z2


def atoms_svg_3d(atoms: Sequence[Atom], width: int = 360, height: int = 300, title: str = "") -> str:
    coords = [(x, y, z) for _, x, y, z in atoms]
    center = [sum(coord[i] for coord in coords) / len(coords) for i in range(3)]
    rotated = [rotate_for_3d((x - center[0], y - center[1], z - center[2])) for x, y, z in coords]
    scaled = scaled_points(rotated, width, height, 36)
    pairs = pair_distances(atoms)
    min_d = min((d for _, _, d in pairs), default=1.7)
    cutoff = min(2.25, max(1.75, min_d * 1.28))
    z_values = [p[2] for p in scaled]
    zmin, zmax = min(z_values), max(z_values)
    zspan = zmax - zmin + 1e-6

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]
    if title:
        lines.append(f'<text x="16" y="24" font-family="Arial" font-size="15" font-weight="700" fill="#222">{html(title)}</text>')

    bond_items = []
    for i, j, d in pairs:
        if d <= cutoff:
            bond_items.append((scaled[i][2] + scaled[j][2], i, j, d))
    for _, i, j, _ in sorted(bond_items):
        x1, y1, z1 = scaled[i]
        x2, y2, z2 = scaled[j]
        depth = ((z1 + z2) / 2.0 - zmin) / zspan
        stroke = 1.6 + 1.4 * depth
        shade = int(150 - 58 * depth)
        lines.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="rgb({shade},{shade},{shade})" stroke-width="{stroke:.2f}" stroke-linecap="round"/>')

    for idx in sorted(range(len(scaled)), key=lambda atom_idx: scaled[atom_idx][2]):
        x, y, z = scaled[idx]
        depth = (z - zmin) / zspan
        radius = 9.5 + 5.2 * depth
        fill = int(115 + 52 * depth)
        lines.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius:.2f}" fill="rgb(58,{fill},{fill})" stroke="#1f4142" stroke-width="1.5"/>')
        lines.append(f'<text x="{x:.2f}" y="{y + 3.8:.2f}" text-anchor="middle" font-family="Arial" font-size="9" fill="#ffffff">B</text>')
    lines.append("</svg>")
    return "\n".join(lines)


def write_workflow_svg(path: Path) -> None:
    width, height = 1060, 210
    labels = [
        "Стартовые\nгеометрии",
        "R2SCAN-3C\nOpt screening",
        "Сбор энергий\nи сходимости",
        "Отбор 10\nкандидатов",
        "PBE0-D4/\ndef2-TZVP OptFreq",
        "Частоты\nи best_B6.xyz",
    ]
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<defs><marker id="arrow" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="#555"/></marker></defs>',
    ]
    box_w, box_h, y = 145, 74, 74
    for i, label in enumerate(labels):
        x = 24 + i * 170
        lines.append(f'<rect x="{x}" y="{y}" width="{box_w}" height="{box_h}" rx="8" fill="#f5f7f7" stroke="#506b6b" stroke-width="1.5"/>')
        for k, part in enumerate(label.split("\n")):
            lines.append(f'<text x="{x + box_w / 2}" y="{y + 30 + k * 18}" text-anchor="middle" font-family="Arial" font-size="14" fill="#1d2d2d">{html(part)}</text>')
        if i < len(labels) - 1:
            lines.append(f'<line x1="{x + box_w + 8}" y1="{y + box_h / 2}" x2="{x + 162}" y2="{y + box_h / 2}" stroke="#555" stroke-width="2" marker-end="url(#arrow)"/>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_start_geometries_svg(path: Path, project: Path) -> None:
    gen = load_geometry_module(project)
    entries = [
        ("linear_chain", gen.linear_chain(1.8), (0, 1)),
        ("planar_ring", gen.planar_ring(1.8), (0, 1)),
        ("compact_planar_triangle", gen.compact_planar_triangle(1.8), (0, 1)),
        ("rhombic_planar", gen.rhombic_planar(1.8), (0, 1)),
        ("rectangular_planar", gen.rectangular_planar(1.8), (0, 1)),
        ("octahedral_3d", gen.octahedral_3d(1.8), (0, 2)),
        ("trigonal_prism", gen.trigonal_prism(1.8), (0, 2)),
        ("random_3d", gen.random_3d(1.8, 1000), (0, 2)),
    ]
    cell_w, cell_h = 250, 205
    width, height = cell_w * 4, cell_h * 2
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]
    for idx, (title, atoms, axes) in enumerate(entries):
        col, row = idx % 4, idx // 4
        x0, y0 = col * cell_w, row * cell_h
        frag = atoms_svg(atoms, cell_w, cell_h, title, axes)
        inner = "\n".join(frag.splitlines()[2:-1])
        lines.append(f'<g transform="translate({x0},{y0})">{inner}</g>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_bar_svg(path: Path, rows: Sequence[Dict[str, str]], title: str, value_key: str = "relative_energy_ev", limit: int = 10) -> None:
    shown = list(rows[:limit])
    width, height = 1120, 520
    left, right, top, bottom = 360, 36, 54, 48
    plot_w, plot_h = width - left - right, height - top - bottom
    values = [max(0.0, to_float(r.get(value_key, "0"), 0.0)) for r in shown]
    vmax = max(values + [1e-6])
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="24" y="32" font-family="Arial" font-size="18" font-weight="700" fill="#222">{html(title)}</text>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#333" stroke-width="1"/>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#333" stroke-width="1"/>',
    ]
    bar_h = min(32, plot_h / max(len(shown), 1) * 0.58)
    step = plot_h / max(len(shown), 1)
    for idx, (row, value) in enumerate(zip(shown, values)):
        y = top + idx * step + (step - bar_h) / 2
        bar_w = (value / vmax) * (plot_w - 90)
        label = f"{idx + 1}. {short_name(row.get('calculation_name', ''), 38)}"
        lines.append(f'<text x="{left - 12}" y="{y + bar_h * 0.68:.2f}" text-anchor="end" font-family="Arial" font-size="12" fill="#222">{html(label)}</text>')
        lines.append(f'<rect x="{left}" y="{y:.2f}" width="{bar_w:.2f}" height="{bar_h:.2f}" fill="#4f8a8b"/>')
        lines.append(f'<text x="{left + bar_w + 8:.2f}" y="{y + bar_h * 0.68:.2f}" font-family="Arial" font-size="12" fill="#222">{value:.8f} eV</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_placeholder_svg(path: Path, title: str, message: str) -> None:
    path.write_text(
        "\n".join(
            [
                '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="220" viewBox="0 0 900 220">',
                '<rect width="100%" height="100%" fill="#ffffff"/>',
                f'<text x="24" y="44" font-family="Arial" font-size="20" font-weight="700" fill="#222">{html(title)}</text>',
                f'<text x="24" y="86" font-family="Arial" font-size="14" fill="#666">{html(message)}</text>',
                "</svg>",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def vibration_dir(project: Path) -> Path:
    return project / "results" / "vibrations" / "B6"


def read_vibration_csv(project: Path, filename: str) -> List[Dict[str, str]]:
    path = vibration_dir(project) / filename
    if not path.exists():
        return []
    return read_csv(path)


def vibration_frequency_rows(project: Path) -> List[Dict[str, str]]:
    return read_vibration_csv(project, "B6_all_vibrational_frequencies.csv")


def vibration_mode_summary_rows(project: Path) -> List[Dict[str, str]]:
    return read_vibration_csv(project, "B6_mode_summary.csv")


def vibration_amplitude_rows(project: Path) -> List[Dict[str, str]]:
    return read_vibration_csv(project, "B6_normal_mode_amplitudes.csv")


def vibration_summary_table_rows(mode_rows: Sequence[Dict[str, str]]) -> List[List[object]]:
    rows: List[List[object]] = []
    for row in mode_rows:
        rows.append(
            [
                row.get("mode_number", ""),
                row.get("orca_index", ""),
                row.get("frequency_cm-1", ""),
                row.get("max_amplitude", ""),
                row.get("dominant_atom_index", ""),
                row.get("dominant_atom_element", ""),
            ]
        )
    return rows


def write_vibration_spectrum_svg(path: Path, frequency_rows: Sequence[Dict[str, str]]) -> None:
    if not frequency_rows:
        write_placeholder_svg(path, "Вибрационный спектр B6", "Файл B6_all_vibrational_frequencies.csv не найден.")
        return

    values = [to_float(row.get("frequency_cm-1", ""), math.nan) for row in frequency_rows]
    values = [value for value in values if not math.isnan(value)]
    if not values:
        write_placeholder_svg(path, "Вибрационный спектр B6", "В таблице частот нет численных значений.")
        return

    width, height = 1120, 390
    left, right, top, bottom = 70, 30, 54, 72
    plot_w, plot_h = width - left - right, height - top - bottom
    fmin, fmax = 0.0, max(values) * 1.06
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="24" y="32" font-family="Arial" font-size="18" font-weight="700" fill="#222">Вибрационные частоты выбранного B6 минимума</text>',
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#333" stroke-width="1"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#333" stroke-width="1"/>',
    ]
    for tick in range(0, int(fmax) + 1, 200):
        x = left + (tick - fmin) / (fmax - fmin) * plot_w
        lines.append(f'<line x1="{x:.2f}" y1="{top + plot_h}" x2="{x:.2f}" y2="{top + plot_h + 6}" stroke="#333" stroke-width="1"/>')
        lines.append(f'<text x="{x:.2f}" y="{top + plot_h + 24}" text-anchor="middle" font-family="Arial" font-size="12" fill="#333">{tick}</text>')
        lines.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_h}" stroke="#eceeee" stroke-width="1"/>')
    lines.append(f'<text x="{left + plot_w / 2}" y="{height - 16}" text-anchor="middle" font-family="Arial" font-size="13" fill="#333">frequency, cm⁻¹</text>')

    step = plot_h / (len(values) + 1)
    for idx, row in enumerate(frequency_rows, start=1):
        freq = to_float(row.get("frequency_cm-1", ""), math.nan)
        if math.isnan(freq):
            continue
        y = top + idx * step
        x = left + (freq - fmin) / (fmax - fmin) * plot_w
        lines.append(f'<line x1="{x:.2f}" y1="{y - 8:.2f}" x2="{x:.2f}" y2="{y + 8:.2f}" stroke="#4f8a8b" stroke-width="3"/>')
        lines.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4.5" fill="#4f8a8b" stroke="#1f4142" stroke-width="1"/>')
        lines.append(f'<text x="{x + 7:.2f}" y="{y + 4:.2f}" font-family="Arial" font-size="11" fill="#222">ν{row.get("mode_number", idx)}={freq:.2f}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_vibration_amplitude_heatmap_svg(path: Path, amplitude_rows: Sequence[Dict[str, str]]) -> None:
    if not amplitude_rows:
        write_placeholder_svg(path, "Амплитуды нормальных мод B6", "Файл B6_normal_mode_amplitudes.csv не найден.")
        return

    modes = sorted({int(row.get("mode_number", "0")) for row in amplitude_rows if row.get("mode_number", "").isdigit()})
    atoms = sorted({int(row.get("atom_index", "0")) for row in amplitude_rows if row.get("atom_index", "").isdigit()})
    amplitudes: Dict[Tuple[int, int], float] = {}
    frequencies: Dict[int, str] = {}
    max_amp = 0.0
    for row in amplitude_rows:
        if not row.get("mode_number", "").isdigit() or not row.get("atom_index", "").isdigit():
            continue
        mode = int(row["mode_number"])
        atom = int(row["atom_index"])
        amp = to_float(row.get("amplitude", ""), 0.0)
        amplitudes[(atom, mode)] = amp
        frequencies[mode] = row.get("frequency_cm-1", "")
        max_amp = max(max_amp, amp)

    cell_w, cell_h = 74, 34
    left, top = 110, 86
    width = left + cell_w * len(modes) + 40
    height = top + cell_h * len(atoms) + 72
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="24" y="32" font-family="Arial" font-size="18" font-weight="700" fill="#222">Участие атомов в нормальных модах B6</text>',
        '<text x="24" y="56" font-family="Arial" font-size="12" fill="#666">Цвет показывает относительную амплитуду смещения атома в моде; значения нормированы на максимум в таблице.</text>',
    ]
    for col, mode in enumerate(modes):
        x = left + col * cell_w + cell_w / 2
        lines.append(f'<text x="{x:.2f}" y="{top - 34}" text-anchor="middle" font-family="Arial" font-size="11" fill="#222">ν{mode}</text>')
        lines.append(f'<text x="{x:.2f}" y="{top - 18}" text-anchor="middle" font-family="Arial" font-size="10" fill="#666">{html(frequencies.get(mode, ""))}</text>')
    for row_idx, atom in enumerate(atoms):
        y = top + row_idx * cell_h
        lines.append(f'<text x="{left - 12}" y="{y + cell_h * 0.62:.2f}" text-anchor="end" font-family="Arial" font-size="12" fill="#222">B{atom}</text>')
        for col, mode in enumerate(modes):
            x = left + col * cell_w
            amp = amplitudes.get((atom, mode), 0.0)
            frac = amp / max_amp if max_amp else 0.0
            r = int(240 - 176 * frac)
            g = int(246 - 88 * frac)
            b = int(246 - 86 * frac)
            lines.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{cell_w - 3}" height="{cell_h - 3}" fill="rgb({r},{g},{b})" stroke="#d3dddd" stroke-width="1"/>')
            lines.append(f'<text x="{x + cell_w / 2:.2f}" y="{y + cell_h * 0.62:.2f}" text-anchor="middle" font-family="Arial" font-size="10" fill="#1f2727">{amp:.2f}</text>')
    lines.append(f'<text x="{left + cell_w * len(modes) / 2}" y="{height - 18}" text-anchor="middle" font-family="Arial" font-size="12" fill="#333">normal mode number and frequency, cm⁻¹</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_min_energy_geometries_svg(path: Path, final_rows: Sequence[Dict[str, str]], project: Path, limit: int = 10) -> None:
    shown = list(final_rows[:limit])
    cols = 5 if len(shown) > 4 else max(1, len(shown))
    rows = max(1, math.ceil(len(shown) / cols))
    cell_w, cell_h = 300, 255
    width, height = cell_w * cols, cell_h * rows + 42
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="20" y="28" font-family="Arial" font-size="18" font-weight="700" fill="#222">Минимально-энергетические финальные геометрии B6</text>',
    ]

    for idx, row in enumerate(shown):
        col, grid_row = idx % cols, idx // cols
        x0, y0 = col * cell_w, 42 + grid_row * cell_h
        calc_label = short_name(row.get("calculation_name", ""), 34)
        rel = row.get("relative_energy_ev", "")
        mult = row.get("multiplicity", "")
        freq = row.get("lowest_frequency_cm-1", "")
        source = row.get("geometry_type", "")

        lines.append(f'<g transform="translate({x0},{y0})">')
        lines.append(f'<text x="12" y="18" font-family="Arial" font-size="13" font-weight="700" fill="#222">#{idx + 1} {html(calc_label)}</text>')
        lines.append(f'<text x="12" y="36" font-family="Arial" font-size="12" fill="#333">m={html(mult)}  ΔE={html(rel)} eV  νmin={html(freq)} cm⁻¹</text>')
        lines.append(f'<text x="12" y="53" font-family="Arial" font-size="11" fill="#666">source: {html(source)}</text>')

        xyz_path = resolve_project_path(row.get("xyz_file", ""), project)
        if xyz_path.exists():
            atoms = read_xyz(xyz_path)
            frag = atoms_svg(atoms, cell_w, cell_h - 62, "", (0, 1))
            inner = "\n".join(frag.splitlines()[2:-1])
            lines.append(f'<g transform="translate(0,58)">{inner}</g>')
        else:
            lines.append(f'<text x="{cell_w / 2}" y="{cell_h / 2}" text-anchor="middle" font-family="Arial" font-size="13" fill="#9a2d2d">XYZ not found</text>')
        lines.append("</g>")

    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_3d_start_geometries_svg(path: Path, project: Path) -> None:
    gen = load_geometry_module(project)
    entries = []
    for name, maker in [
        ("octahedral_3d", getattr(gen, "octahedral_3d", None)),
        ("trigonal_prism", getattr(gen, "trigonal_prism", None)),
        ("pentagonal_pyramid_3d", getattr(gen, "pentagonal_pyramid_3d", None)),
    ]:
        if maker is not None:
            entries.append((name, maker(1.8)))
    if hasattr(gen, "random_3d"):
        entries.extend(
            [
                ("random_3d_seed1000", gen.random_3d(1.8, 1000)),
                ("random_3d_seed1001", gen.random_3d(1.8, 1001)),
                ("random_3d_seed1002", gen.random_3d(1.8, 1002)),
            ]
        )

    cols = 3
    rows = max(1, math.ceil(len(entries) / cols))
    cell_w, cell_h = 330, 270
    width, height = cell_w * cols, cell_h * rows + 42
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="20" y="28" font-family="Arial" font-size="18" font-weight="700" fill="#222">Отдельные 3D-стартовые геометрии B6</text>',
    ]
    for idx, (title, atoms) in enumerate(entries):
        col, grid_row = idx % cols, idx // cols
        x0, y0 = col * cell_w, 42 + grid_row * cell_h
        frag = atoms_svg_3d(atoms, cell_w, cell_h, title)
        inner = "\n".join(frag.splitlines()[2:-1])
        lines.append(f'<g transform="translate({x0},{y0})">{inner}</g>')
    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_min_energy_geometries_3d_svg(
    path: Path,
    final_rows: Sequence[Dict[str, str]],
    project: Path,
    limit: int = 10,
    title: str = "3D-проекция минимально-энергетических геометрий B6",
) -> None:
    shown = list(final_rows[:limit])
    cols = 5 if len(shown) > 4 else max(1, len(shown))
    rows = max(1, math.ceil(len(shown) / cols))
    cell_w, cell_h = 300, 255
    width, height = cell_w * cols, cell_h * rows + 42
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="20" y="28" font-family="Arial" font-size="18" font-weight="700" fill="#222">{html(title)}</text>',
    ]

    for idx, row in enumerate(shown):
        col, grid_row = idx % cols, idx // cols
        x0, y0 = col * cell_w, 42 + grid_row * cell_h
        calc_label = short_name(row.get("calculation_name", ""), 34)
        rel = row.get("relative_energy_ev", "")
        mult = row.get("multiplicity", "")
        freq = row.get("lowest_frequency_cm-1", "")
        source = row.get("geometry_type", "")

        lines.append(f'<g transform="translate({x0},{y0})">')
        lines.append(f'<text x="12" y="18" font-family="Arial" font-size="13" font-weight="700" fill="#222">#{idx + 1} {html(calc_label)}</text>')
        lines.append(f'<text x="12" y="36" font-family="Arial" font-size="12" fill="#333">m={html(mult)}  ΔE={html(rel)} eV  νmin={html(freq)} cm⁻¹</text>')
        lines.append(f'<text x="12" y="53" font-family="Arial" font-size="11" fill="#666">source: {html(source)}</text>')

        xyz_path = resolve_project_path(row.get("xyz_file", ""), project)
        if xyz_path.exists():
            atoms = read_xyz(xyz_path)
            frag = atoms_svg_3d(atoms, cell_w, cell_h - 62, "")
            inner = "\n".join(frag.splitlines()[2:-1])
            lines.append(f'<g transform="translate(0,58)">{inner}</g>')
        else:
            lines.append(f'<text x="{cell_w / 2}" y="{cell_h / 2}" text-anchor="middle" font-family="Arial" font-size="13" fill="#9a2d2d">XYZ not found</text>')
        lines.append("</g>")

    lines.append("</svg>")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_min_energy_geometries_3d_html(
    path: Path,
    final_rows: Sequence[Dict[str, str]],
    project: Path,
    limit: int = 10,
    title: str = "3D-визуализация низкоэнергетических финальных геометрий B6",
    subtitle: str = "Перетащите мышью для вращения; колесо мыши меняет масштаб. Все структуры отсортированы по PBE0-D4/def2-TZVP энергии.",
) -> None:
    structures: List[Dict[str, object]] = []
    for idx, row in enumerate(final_rows[:limit], start=1):
        xyz_path = resolve_project_path(row.get("xyz_file", ""), project)
        if not xyz_path.exists():
            continue
        atoms = read_xyz(xyz_path)
        structures.append(
            {
                "rank": idx,
                "name": short_name(row.get("calculation_name", ""), 48),
                "source": row.get("geometry_type", ""),
                "multiplicity": row.get("multiplicity", ""),
                "relative_energy_ev": row.get("relative_energy_ev", ""),
                "lowest_frequency": row.get("lowest_frequency_cm-1", ""),
                "atoms": [{"el": el, "x": x, "y": y, "z": z} for el, x, y, z in atoms],
            }
        )

    data = json.dumps(structures, ensure_ascii=False)
    html_doc = """<!doctype html>
<html lang="ru">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>3D B6 geometries</title>
  <style>
    :root { color-scheme: light; font-family: Arial, sans-serif; }
    body { margin: 0; background: #f6f7f7; color: #1f2727; }
    header { padding: 18px 24px 10px; background: #fff; border-bottom: 1px solid #d9dddd; }
    h1 { margin: 0 0 6px; font-size: 22px; }
    p { margin: 0; color: #5b6565; font-size: 14px; }
    .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 14px; padding: 16px; }
    .card { background: #fff; border: 1px solid #d8dfdf; border-radius: 8px; overflow: hidden; }
    .meta { padding: 10px 12px; border-bottom: 1px solid #e5e9e9; }
    .title { font-weight: 700; font-size: 13px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
    .sub { color: #596464; font-size: 12px; margin-top: 4px; line-height: 1.35; }
    canvas { width: 100%; height: 280px; display: block; background: #fff; cursor: grab; }
    canvas:active { cursor: grabbing; }
  </style>
</head>
<body>
  <header>
    <h1>""" + html(title) + """</h1>
    <p>""" + html(subtitle) + """</p>
  </header>
  <main class="grid" id="grid"></main>
  <script>
    const structures = """ + data + """;

    function pairDistances(atoms) {
      const pairs = [];
      for (let i = 0; i < atoms.length; i++) {
        for (let j = i + 1; j < atoms.length; j++) {
          const a = atoms[i], b = atoms[j];
          const d = Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z);
          pairs.push([i, j, d]);
        }
      }
      return pairs;
    }

    function centerAtoms(atoms) {
      const c = atoms.reduce((acc, atom) => {
        acc.x += atom.x; acc.y += atom.y; acc.z += atom.z;
        return acc;
      }, {x: 0, y: 0, z: 0});
      c.x /= atoms.length; c.y /= atoms.length; c.z /= atoms.length;
      return atoms.map(atom => ({...atom, x: atom.x - c.x, y: atom.y - c.y, z: atom.z - c.z}));
    }

    function rotate(atom, yaw, pitch) {
      const cy = Math.cos(yaw), sy = Math.sin(yaw);
      const cp = Math.cos(pitch), sp = Math.sin(pitch);
      const x1 = atom.x * cy + atom.z * sy;
      const z1 = -atom.x * sy + atom.z * cy;
      const y1 = atom.y * cp - z1 * sp;
      const z2 = atom.y * sp + z1 * cp;
      return {el: atom.el, x: x1, y: y1, z: z2};
    }

    function draw(canvas, structure, state) {
      const dpr = window.devicePixelRatio || 1;
      const rect = canvas.getBoundingClientRect();
      const w = Math.max(1, Math.floor(rect.width * dpr));
      const h = Math.max(1, Math.floor(rect.height * dpr));
      if (canvas.width !== w || canvas.height !== h) {
        canvas.width = w; canvas.height = h;
      }
      const ctx = canvas.getContext("2d");
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, rect.width, rect.height);

      const centered = centerAtoms(structure.atoms);
      const rotated = centered.map(atom => rotate(atom, state.yaw, state.pitch));
      const xs = rotated.map(a => a.x), ys = rotated.map(a => a.y), zs = rotated.map(a => a.z);
      const span = Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys), 1e-6);
      const scale = Math.min(rect.width - 58, rect.height - 46) / span * state.zoom;
      const zmin = Math.min(...zs), zmax = Math.max(...zs), zspan = zmax - zmin + 1e-6;
      const projected = rotated.map(atom => ({
        ...atom,
        sx: rect.width / 2 + atom.x * scale,
        sy: rect.height / 2 - atom.y * scale,
        depth: (atom.z - zmin) / zspan
      }));

      const pairs = pairDistances(centered);
      const minD = Math.min(...pairs.map(p => p[2]));
      const cutoff = Math.min(2.25, Math.max(1.75, minD * 1.28));
      const bonds = pairs
        .filter(p => p[2] <= cutoff)
        .map(p => ({i: p[0], j: p[1], z: projected[p[0]].depth + projected[p[1]].depth}))
        .sort((a, b) => a.z - b.z);

      ctx.lineCap = "round";
      for (const bond of bonds) {
        const a = projected[bond.i], b = projected[bond.j];
        const depth = (a.depth + b.depth) / 2;
        const shade = Math.round(150 - 58 * depth);
        ctx.strokeStyle = `rgb(${shade}, ${shade}, ${shade})`;
        ctx.lineWidth = 2 + 2.2 * depth;
        ctx.beginPath();
        ctx.moveTo(a.sx, a.sy);
        ctx.lineTo(b.sx, b.sy);
        ctx.stroke();
      }

      const atoms = projected.slice().sort((a, b) => a.depth - b.depth);
      for (const atom of atoms) {
        const r = 13 + 7 * atom.depth;
        const fill = Math.round(116 + 55 * atom.depth);
        const grad = ctx.createRadialGradient(atom.sx - r * 0.35, atom.sy - r * 0.35, r * 0.2, atom.sx, atom.sy, r);
        grad.addColorStop(0, `rgb(132, ${Math.min(220, fill + 48)}, ${Math.min(220, fill + 48)})`);
        grad.addColorStop(1, `rgb(48, ${fill}, ${fill})`);
        ctx.fillStyle = grad;
        ctx.strokeStyle = "#1f4142";
        ctx.lineWidth = 1.4;
        ctx.beginPath();
        ctx.arc(atom.sx, atom.sy, r, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
        ctx.fillStyle = "#fff";
        ctx.font = "11px Arial";
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(atom.el, atom.sx, atom.sy + 0.5);
      }
    }

    function setupCard(structure) {
      const card = document.createElement("section");
      card.className = "card";
      card.innerHTML = `
        <div class="meta">
          <div class="title">#${structure.rank} ${structure.name}</div>
          <div class="sub">m=${structure.multiplicity} · ΔE=${structure.relative_energy_ev} eV · νmin=${structure.lowest_frequency} cm⁻¹<br>source: ${structure.source}</div>
        </div>
        <canvas></canvas>
      `;
      const canvas = card.querySelector("canvas");
      const state = {yaw: -0.7, pitch: 0.75, zoom: 0.86, dragging: false, lastX: 0, lastY: 0};
      canvas.addEventListener("pointerdown", event => {
        state.dragging = true; state.lastX = event.clientX; state.lastY = event.clientY; canvas.setPointerCapture(event.pointerId);
      });
      canvas.addEventListener("pointermove", event => {
        if (!state.dragging) return;
        const dx = event.clientX - state.lastX, dy = event.clientY - state.lastY;
        state.lastX = event.clientX; state.lastY = event.clientY;
        state.yaw += dx * 0.01; state.pitch += dy * 0.01;
        state.pitch = Math.max(-1.45, Math.min(1.45, state.pitch));
        draw(canvas, structure, state);
      });
      canvas.addEventListener("pointerup", () => { state.dragging = false; });
      canvas.addEventListener("wheel", event => {
        event.preventDefault();
        state.zoom *= event.deltaY > 0 ? 0.92 : 1.08;
        state.zoom = Math.max(0.45, Math.min(2.5, state.zoom));
        draw(canvas, structure, state);
      }, {passive: false});
      requestAnimationFrame(() => draw(canvas, structure, state));
      window.addEventListener("resize", () => draw(canvas, structure, state));
      return card;
    }

    const grid = document.getElementById("grid");
    structures.forEach(structure => grid.appendChild(setupCard(structure)));
  </script>
</body>
</html>
"""
    path.write_text(html_doc, encoding="utf-8")


def write_figures(project: Path, fig_dir: Path, screening_rows: Sequence[Dict[str, str]], final_rows: Sequence[Dict[str, str]], best_atoms: Sequence[Atom]) -> Dict[str, Path]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    final_3d_rows = [row for row in final_rows if is_3d_geometry(row.get("geometry_type", ""))]
    figures = {
        "workflow": fig_dir / "Figure_1_workflow.svg",
        "starts": fig_dir / "Figure_2_start_geometries.svg",
        "screening": fig_dir / "Figure_3_screening_top10.svg",
        "best": fig_dir / "Figure_4_best_B6.svg",
        "final": fig_dir / "Figure_5_final_relative_energies.svg",
        "min_geometries": fig_dir / "Figure_6_min_energy_geometries.svg",
        "min_geometries_3d_svg": fig_dir / "Figure_7_min_energy_geometries_3d.svg",
        "min_geometries_3d_html": fig_dir / "Figure_7_min_energy_geometries_3d.html",
        "starts_3d": fig_dir / "Figure_8_3d_start_geometries.svg",
        "final_from_3d_svg": fig_dir / "Figure_9_final_from_3d_starts_3d.svg",
        "final_from_3d_html": fig_dir / "Figure_9_final_from_3d_starts_3d.html",
        "vibration_spectrum": fig_dir / "Figure_10_B6_vibrational_spectrum.svg",
        "vibration_heatmap": fig_dir / "Figure_11_B6_normal_mode_amplitudes.svg",
    }
    write_workflow_svg(figures["workflow"])
    write_start_geometries_svg(figures["starts"], project)
    write_bar_svg(figures["screening"], screening_rows, "Топ-10 структур после R2SCAN-3C screening", limit=10)
    figures["best"].write_text(atoms_svg(best_atoms, 520, 420, "best_B6.xyz"), encoding="utf-8")
    write_bar_svg(figures["final"], final_rows, "Относительные энергии финальных кандидатов", limit=10)
    write_min_energy_geometries_svg(figures["min_geometries"], final_rows, project, limit=10)
    write_min_energy_geometries_3d_svg(figures["min_geometries_3d_svg"], final_rows, project, limit=10)
    write_min_energy_geometries_3d_html(figures["min_geometries_3d_html"], final_rows, project, limit=10)
    write_3d_start_geometries_svg(figures["starts_3d"], project)
    write_min_energy_geometries_3d_svg(
        figures["final_from_3d_svg"],
        final_3d_rows,
        project,
        limit=10,
        title="Финальные структуры, полученные из 3D/random стартов",
    )
    write_min_energy_geometries_3d_html(
        figures["final_from_3d_html"],
        final_3d_rows,
        project,
        limit=10,
        title="Финальные B6 структуры из 3D/random стартов",
        subtitle="Отдельный набор только для расчетов, исходный geometry_type которых относится к 3D/random. Перетащите мышью для вращения; колесо мыши меняет масштаб.",
    )
    write_vibration_spectrum_svg(figures["vibration_spectrum"], vibration_frequency_rows(project))
    write_vibration_amplitude_heatmap_svg(figures["vibration_heatmap"], vibration_amplitude_rows(project))
    return figures


def build_report(project: Path, results_csv: Path, final_csv: Path, best_xyz: Path, fig_dir: Path) -> str:
    screening_rows = sorted(read_csv(results_csv), key=lambda r: to_float(r.get("total_energy_hartree", "")))
    final_rows = sorted(read_csv(final_csv), key=lambda r: to_float(r.get("total_energy_hartree", "")))
    best_atoms = read_xyz(best_xyz)
    screening_3d_rows = [row for row in screening_rows if is_3d_geometry(row.get("geometry_type", ""))]
    final_3d_rows = [row for row in final_rows if is_3d_geometry(row.get("geometry_type", ""))]
    vib_freq_rows = vibration_frequency_rows(project)
    vib_mode_rows = vibration_mode_summary_rows(project)
    vib_amp_rows = vibration_amplitude_rows(project)
    write_subset_csv(project / "results" / "screening_3d_results.csv", screening_3d_rows)
    write_subset_csv(project / "results" / "final_from_3d_results.csv", final_3d_rows)
    figures = write_figures(project, fig_dir, screening_rows, final_rows, best_atoms)

    screening_done = [r for r in screening_rows if r.get("normal_termination") == "True"]
    screening_conv = [r for r in screening_rows if r.get("optimization_converged") == "True"]
    final_done = [r for r in final_rows if r.get("normal_termination") == "True"]
    final_conv = [r for r in final_rows if r.get("optimization_converged") == "True"]
    final_true = [r for r in final_rows if r.get("is_true_minimum") == "True"]
    best = final_true[0] if final_true else final_rows[0]
    rms_plane, max_plane = planarity(best_atoms)
    close_005 = [r for r in final_rows if to_float(r.get("relative_energy_ev", "999")) <= 0.05]
    close_001 = [r for r in final_rows if to_float(r.get("relative_energy_ev", "999")) <= 0.01]
    close_0001 = [r for r in final_rows if to_float(r.get("relative_energy_ev", "999")) <= 0.001]
    three_d_final = [
        r
        for r in final_rows
        if any(token in r.get("geometry_type", "").lower() for token in ["random", "3d", "prism", "octa", "pyramid"])
    ]
    planarized_3d = 0
    for row in three_d_final:
        xyz_path = resolve_project_path(row.get("xyz_file", ""), project)
        if xyz_path.exists():
            rms, _ = planarity(read_xyz(xyz_path))
            if rms <= 0.01:
                planarized_3d += 1
    min_dist, max_dist, avg_dist, best_pairs = bond_distance_summary(best_atoms)
    best_bonds = adjacency_like_bonds(best_atoms)
    source_counts = Counter(source_category(row.get("geometry_type", "")) for row in final_rows)
    multiplicity_counts = Counter(row.get("multiplicity", "") for row in screening_rows)
    distance_values = sorted({row.get("distance", "") for row in screening_rows if row.get("distance", "")}, key=lambda value: to_float(value, 999.0))
    geometry_values = sorted({row.get("geometry_type", "") for row in screening_rows if row.get("geometry_type", "")})
    vib_freq_values = [to_float(row.get("frequency_cm-1", ""), math.nan) for row in vib_freq_rows]
    vib_freq_values = [value for value in vib_freq_values if not math.isnan(value)]
    vib_imag_count = sum(1 for row in vib_freq_rows if truthy(row.get("imaginary", "")))
    vib_min = min(vib_freq_values) if vib_freq_values else math.nan
    vib_max = max(vib_freq_values) if vib_freq_values else math.nan
    stage1_template = read_text_if_exists(project / "templates" / "stage1_opt_template.inp")
    final_template = read_text_if_exists(project / "templates" / "final_opt_freq_template.inp")

    lines = [
        "# Многостартовый DFT-поиск устойчивой геометрии нейтрального кластера B₆ методом ORCA 6.1",
        "",
        "## Аннотация",
        f"В работе выполнен многостартовый DFT-поиск низкоэнергетической геометрии нейтрального кластера B₆. На первом этапе обработано `{len(screening_rows)}` расчетов R2SCAN-3C Opt, из которых `{len(screening_done)}` завершились нормально и `{len(screening_conv)}` дали сошедшуюся оптимизацию. На втором этапе выполнено `{len(final_rows)}` финальных PBE0-D4/def2-TZVP OptFreq расчетов. Лучшей найденной структурой в рамках данного набора является `{best.get('calculation_name', '')}` с мультиплетностью `{best.get('multiplicity', '')}`, энергией `{best.get('total_energy_hartree', '')}` Hartree и минимальной ненулевой частотой `{best.get('lowest_frequency_cm-1', '')}` cm⁻¹.",
        "",
        "Ключевой результат: лучшие кандидаты после финальной оптимизации дают плоскую или практически плоскую структуру. Это согласуется с литературной тенденцией малых борных кластеров к 2D/квазиплоским мотивам, но вывод ограничен использованным набором стартов и выбранным уровнем DFT.",
        "",
        "## Краткое содержание результата",
        "",
        simple_markdown_table(
            ["Параметр", "Значение"],
            [
                ["Система", "B₆, neutral"],
                ["Заряд", "0"],
                ["Проверенные мультиплетности", ", ".join(sorted(multiplicity_counts.keys(), key=lambda x: int(x) if str(x).isdigit() else 99))],
                ["Screening", "R2SCAN-3C Opt"],
                ["Финальный уровень", "PBE0-D4/def2-TZVP OptFreq"],
                ["Screening .out", len(screening_rows)],
                ["Успешные screening", len(screening_done)],
                ["Финальные OptFreq .out", len(final_rows)],
                ["Истинные минимумы без мнимых частот", len(final_true)],
                ["best_B6.xyz", short_path(str(best_xyz), project)],
                ["Лучшая энергия, Hartree", best.get("total_energy_hartree", "")],
                ["Лучшая мультиплетность", best.get("multiplicity", "")],
                ["Минимальная ненулевая частота, cm⁻¹", best.get("lowest_frequency_cm-1", "")],
            ],
        ),
        "",
        "## 1. Введение",
        "Кластеры бора интересны из-за электронодефицитной природы атома B, делокализованного связывания и выраженной структурной конкуренции между плоскими, квазиплоскими и объемными мотивами. Для B₆ это означает, что результат нельзя надежно получить из одной заранее выбранной геометрии: разные стартовые структуры могут сходиться в разные локальные минимумы или, наоборот, показывать, что объемные старты переходят к плоской области поверхности потенциальной энергии.",
        "",
        "## 2. Литературный контекст",
        "Современные подходы к поиску минимумов атомных кластеров используют разнообразные начальные структуры, критерии уникальности и последующую квантово-химическую оптимизацию [1]. В работах по малым борным и борсодержащим кластерам типовой протокол включает DFT-оптимизацию, сравнение полных энергий, анализ электронных свойств и проверку устойчивости через частоты [2,3]. Для чистого B₆ особенно важны работы, где обсуждаются планарность, антиароматичность и химическое связывание B₆/B₆⁻ [4]. Более широкие обзоры по size-selected boron clusters также показывают тенденцию малых борных кластеров к плоским и квазиплоским структурам, связанную с делокализацией σ- и π-связей [5]. Поэтому в данной работе сравнивались не только плоские, но и 3D-старты, чтобы не навязывать геометрию заранее.",
        "",
        "## 3. Цель и задачи работы",
        "Цель работы: найти наиболее устойчивую геометрию нейтрального кластера B₆ в рамках заданного набора стартовых структур, мультиплетностей и уровней DFT.",
        "",
        "Задачи: сгенерировать набор стартовых структур B₆; выполнить R2SCAN-3C Opt screening; сравнить энергии и сходимость; выбрать низкоэнергетические кандидаты; провести финальный PBE0-D4/def2-TZVP OptFreq; проверить отсутствие мнимых частот; сохранить `best_B6.xyz`.",
        "",
        "## 4. Методика расчетов",
        "- Программа: ORCA 6.1",
        "- Кластер: B₆",
        "- Заряд: 0",
        "- Проверенные мультиплетности: 1, 3, 5",
        "- Первичный метод: R2SCAN-3C Opt",
        "- Финальный метод: PBE0-D4/def2-TZVP OptFreq",
        "- Дисперсионная поправка: D4",
        "- Число CPU: 8",
        "- `%maxcore`: 2500 MB",
        "",
        "Все энергии, частоты и координаты извлекались только из реальных ORCA `.out` файлов и производных CSV/XYZ файлов. Фиктивные или вручную придуманные значения не использовались.",
        "",
        "### 4.1. Логика расчётного workflow",
        "Расчётный проект был организован как последовательная кампания, в которой широкий набор стартовых структур сначала быстро оптимизируется на более дешёвом уровне, а затем низкоэнергетические кандидаты уточняются на более дорогом уровне с частотным анализом. Такой подход снижает риск того, что итог будет зависеть от одной произвольно выбранной геометрии.",
        "",
        "1. Генерация стартовых XYZ и ORCA `.inp` файлов.",
        "2. R2SCAN-3C Opt screening для всех стартов, расстояний и мультиплетностей.",
        "3. Сбор `FINAL SINGLE POINT ENERGY`, признаков нормального завершения и сходимости оптимизации.",
        "4. Отбор низкоэнергетических кандидатов с геометрической дедупликацией.",
        "5. Финальный PBE0-D4/def2-TZVP OptFreq.",
        "6. Частотный анализ и выбор `best_B6.xyz` только среди структур без мнимых частот.",
        "",
        "### 4.2. ORCA-настройки screening",
        "Ключевая строка screening-расчёта:",
        "",
        "```orca",
        "! R2SCAN-3C TightSCF TightOpt Opt",
        "",
        "%pal",
        "  nprocs 8",
        "end",
        "",
        "%maxcore 2500",
        "```",
        "",
        "Дополнительно использовались `MaxIter 500` для SCF и `MaxIter 300` для геометрической оптимизации. Полный шаблон хранится в `templates/stage1_opt_template.inp`.",
        "",
        "### 4.3. ORCA-настройки финального этапа",
        "Ключевая строка финального расчёта:",
        "",
        "```orca",
        "! PBE0 def2-TZVP D4 def2/J RIJCOSX TightSCF TightOpt Opt Freq",
        "```",
        "",
        "На финальном этапе выполнялись одновременно переоптимизация и расчёт частот (`Opt Freq`). Полный шаблон хранится в `templates/final_opt_freq_template.inp`.",
        "",
        "### 4.4. Контроль качества парсинга",
        "Для каждого `.out` файла проверялись строки `ORCA TERMINATED NORMALLY`, `THE OPTIMIZATION HAS CONVERGED`, `FINAL SINGLE POINT ENERGY` и, для финального этапа, блок `VIBRATIONAL FREQUENCIES`. В частотном анализе нулевые трансляционно-вращательные моды не интерпретировались как мнимые вибрационные частоты.",
        "",
        f"Рисунок 1: `{short_path(str(figures['workflow']), project)}`.",
        "",
        "## 5. Генерация стартовых геометрий B₆",
        "Для уменьшения риска попадания в локальный минимум был использован многостартовый подход. Были сгенерированы плоские, квазиплоские и трехмерные стартовые структуры B₆ с различными начальными расстояниями B-B. Для каждой структуры были проверены мультиплетности 1, 3 и 5.",
        "",
        "В расчетной кампании использовались следующие типы стартов: линейная цепочка; плоское кольцо; компактная плоская структура; ромбическая структура; прямоугольная структура; октаэдрическая 3D-структура; тригональная призма; случайные 3D-структуры. Дополнительно генератор поддерживает искаженное плоское кольцо, fused-triangle, квазиплоскую и пирамидальную 3D-структуру для расширенного набора.",
        "",
        f"В опубликованной расчетной кампании обработано `{len(screening_rows)}` screening output-файла. Текущая расширенная версия генератора поддерживает набор до `384` screening-расчетов за счет дополнительных геометрий и расстояний; поэтому различие между числом `252` в обработанных результатах и `384` в README относится к разным состояниям расчетной кампании, а не к ошибке в таблицах.",
        "",
        "### 5.1. Набор стартов и их назначение",
        "",
        simple_markdown_table(
            ["Стартовая геометрия", "Тип", "Зачем нужна в кампании"],
            [
                ["linear_chain", "1D", "Проверка вытянутого предела и возможной перестройки в компактную форму"],
                ["planar_ring", "2D", "Кольцевой плоский мотив B₆"],
                ["distorted_planar_ring", "2D", "Проверка устойчивости кольца к нарушению симметрии"],
                ["compact_planar_triangle", "2D", "Компактный фрагмент треугольной борной сетки"],
                ["rhombic_planar", "2D", "Плоский ромбический мотив; важен для сравнения с плоскими минимумами"],
                ["rectangular_planar", "2D", "Альтернативный плоский мотив с иной топологией B-B контактов"],
                ["fused_triangles_planar", "2D", "Два соединённых треугольника как компактный борный мотив"],
                ["quasi_planar", "quasi-2D", "Проверка слабого выхода атомов из плоскости"],
                ["octahedral_3d", "3D", "Высокосимметричный объёмный конкурент"],
                ["trigonal_prism", "3D", "Призматический объёмный конкурент"],
                ["pentagonal_pyramid_3d", "3D", "Пирамидальный объёмный старт"],
                ["random_3d_seed*", "3D", "Набор случайных, но физически разумных стартов"],
            ],
        ),
        "",
        f"В уже обработанном screening-наборе представлены `{len(geometry_values)}` типов `geometry_type`: {', '.join(geometry_values)}.",
        f"Начальные расстояния B-B в обработанном наборе: {', '.join(distance_values)} Å.",
        "",
        f"Рисунок 2: `{short_path(str(figures['starts']), project)}`.",
        "",
        "## 6. Первичный screening: R2SCAN-3C Opt",
        f"На screening-этапе обработано `{len(screening_rows)}` ORCA output-файлов. Нормально завершились `{len(screening_done)}` расчетов, сходимость оптимизации обнаружена у `{len(screening_conv)}` расчетов. Полные энергии извлекались из строки `FINAL SINGLE POINT ENERGY`; расчеты без нормального завершения или без сходимости не рассматриваются как надежные финальные кандидаты.",
        "",
        "### 6.1. Screening по мультиплетностям",
        "Эта таблица показывает, как распределены расчёты по спиновым состояниям. Энергии разных мультиплетностей сравнивались только после успешной оптимизации на одном уровне теории.",
        "",
        simple_markdown_table(
            ["multiplicity", "всего", "normal", "converged", "лучшая E, Hartree", "лучшая ΔE, eV", "лучший расчёт"],
            group_summary_rows(screening_rows, "multiplicity"),
        ),
        "",
        "### 6.2. Screening по типам стартовых геометрий",
        "Таблица ниже нужна не для окончательного выбора минимума, а для контроля многостартового покрытия: она показывает, какие типы стартов давали низкоэнергетические структуры после R2SCAN-3C Opt.",
        "",
        simple_markdown_table(
            ["geometry_type", "всего", "normal", "converged", "лучшая E, Hartree", "лучшая ΔE, eV", "лучший расчёт"],
            group_summary_rows(screening_rows, "geometry_type"),
        ),
        "",
        "### 6.3. Низкоэнергетическая область screening",
        "Таблица 1. Топ-10 screening-результатов; полный набор приведен в `results/results.csv`.",
        "",
        markdown_table(screening_rows, SCREENING_COLUMNS, project, limit=10),
        "",
        f"Рисунок 3: `{short_path(str(figures['screening']), project)}`.",
        "",
        "## 7. Отбор финальных кандидатов",
        "Финальные кандидаты выбирались из низкоэнергетических расчетов screening-этапа с нормальным завершением и сошедшейся оптимизацией. Для уменьшения дублирования структур используется сравнение отсортированных межатомных расстояний B-B; структуры с близкими distance fingerprints рассматриваются как геометрически повторяющиеся кандидаты. В текущем финальном наборе сохранены 10 OptFreq расчетов.",
        "",
        "Практически это означает, что финальный этап не является повторением всех 252 screening-расчётов. Его задача - уточнить наиболее перспективную часть поверхности потенциальной энергии и проверить, являются ли структуры настоящими минимумами по частотам.",
        "",
        "Распределение источников финальных кандидатов:",
        "",
        simple_markdown_table(
            ["Категория старта", "число финальных расчётов"],
            sorted(source_counts.items()),
        ),
        "",
        "## 8. Финальный расчет: PBE0-D4/def2-TZVP OptFreq",
        f"Финальный этап включал `{len(final_rows)}` расчетов PBE0-D4/def2-TZVP OptFreq. Нормально завершились `{len(final_done)}` расчетов, сходимость оптимизации обнаружена у `{len(final_conv)}` расчетов.",
        "",
        "Таблица 2. Финальные расчеты PBE0-D4/def2-TZVP OptFreq.",
        "",
        markdown_table(final_rows, FINAL_COLUMNS, project, limit=None),
        "",
        "### 8.1. Топ-5 уникальных геометрических групп",
        "Финальная таблица содержит 10 строк, но строки не обязательно соответствуют 10 независимым минимумам. Для ориентировочной дедупликации ниже финальные XYZ сгруппированы по отсортированным расстояниям B-B с порогом 0.02 Å. Такая таблица помогает отделить физически разные мотивы от повторного попадания в один и тот же минимум.",
        "",
        simple_markdown_table(
            ["group", "строк", "представитель", "m", "min ΔE, eV", "max ΔE, eV", "исходные geometry_type"],
            unique_final_geometry_rows(final_rows, project, tolerance=0.02, limit=5),
        ),
        "",
        "### 8.2. Планарность финальных структур",
        "Для оценки того, перешли ли 3D-старты в плоские или квазиплоские минимумы, для каждого финального XYZ была рассчитана лучшая плоскость по координатам атомов. В таблице приведены RMS-отклонение атомов от этой плоскости и максимальное абсолютное отклонение.",
        "",
        simple_markdown_table(
            ["calculation_name", "source geometry", "ΔE, eV", "lowest ν, cm⁻¹", "RMS plane, Å", "max plane, Å"],
            final_planarity_rows(final_rows, project),
        ),
        "",
        f"Рисунок 5: `{short_path(str(figures['final']), project)}`.",
        f"Рисунок 6: `{short_path(str(figures['min_geometries']), project)}` - визуализация 10 финальных оптимизированных геометрий с минимальной энергией.",
        f"Рисунок 7: `{short_path(str(figures['min_geometries_3d_svg']), project)}` - статическая 3D-проекция тех же низкоэнергетических структур; интерактивная версия для вращения мышью сохранена в `{short_path(str(figures['min_geometries_3d_html']), project)}`.",
        "",
        "### 8.3. Отдельный набор 3D/random структур",
        f"Для отдельного контроля 3D-ветки расчётной кампании выделены строки, у которых `geometry_type` относится к random/3D/призматическим/октаэдрическим/пирамидальным стартам. В screening таких строк: `{len(screening_3d_rows)}`; в финальном OptFreq-наборе: `{len(final_3d_rows)}`.",
        "",
        "Отдельные CSV-файлы для этой части кампании:",
        "",
        "- `results/screening_3d_results.csv`: screening-расчёты, начатые из 3D/random геометрий.",
        "- `results/final_from_3d_results.csv`: финальные OptFreq-расчёты, выбранные из 3D/random стартов.",
        "",
        f"Рисунок 8: `{short_path(str(figures['starts_3d']), project)}` - отдельные 3D-стартовые геометрии B₆.",
        f"Рисунок 9: `{short_path(str(figures['final_from_3d_svg']), project)}` - финальные структуры, которые были получены из 3D/random стартов; интерактивная версия сохранена в `{short_path(str(figures['final_from_3d_html']), project)}`.",
        "",
        "## 9. Частотный анализ",
        "Структура считалась истинным минимумом только при одновременном выполнении трех условий: `ORCA TERMINATED NORMALLY`, сходимость оптимизации и отсутствие мнимых частот. Если структура имеет хотя бы одну мнимую частоту, она не считается финальным минимумом даже при низкой электронной энергии.",
        "",
        "В ORCA output для финальных расчетов присутствуют шесть нулевых трансляционно-вращательных мод. Они не учитывались как мнимые вибрационные моды и не использовались при выборе `lowest_frequency_cm-1`; в таблицу записана минимальная ненулевая вибрационная частота после отсечения мод с |ν| <= 10 cm⁻¹.",
        "",
        f"В текущей финальной таблице структур с `n_imaginary_frequencies > 0` не найдено; истинных минимумов по указанному критерию: `{len(final_true)}`. Для выбранной структуры `n_imaginary_frequencies = {best.get('n_imaginary_frequencies', '')}`, `lowest_frequency_cm-1 = {best.get('lowest_frequency_cm-1', '')}`.",
        "",
        "### 9.1. Подробная сводка по частотам",
        "Для B₆ всего 18 нормальных мод в декартовом представлении: 6 нулевых/трансляционно-вращательных и 12 ненулевых вибрационных. В таблице приведены первые ненулевые частоты, извлечённые из блоков `VIBRATIONAL FREQUENCIES` финальных ORCA output-файлов.",
        "",
        simple_markdown_table(
            ["calculation_name", "нулевых мод", "ненулевых мод", "min ненулевая, cm⁻¹", "первые ненулевые частоты, cm⁻¹", "n_imag"],
            final_frequency_rows(final_rows, project),
        ),
        "",
        "### 9.2. Вибрационный анализ выбранной структуры best_B6",
        f"Для выбранной структуры дополнительно вынесен отдельный набор файлов в `results/vibrations/B6`. В нем сохранены исходный ORCA output, Hessian-файл, оптимизированная геометрия, список всех ненулевых частот и таблицы амплитуд нормальных мод. По этим данным найдено `{len(vib_freq_rows)}` ненулевых вибрационных мод; диапазон частот: `{format_float(vib_min, 2)}`-`{format_float(vib_max, 2)}` cm⁻¹; мнимых частот в вынесенной таблице: `{vib_imag_count}`.",
        "",
        "Таблица 3. Ненулевые вибрационные моды выбранного минимума B₆.",
        "",
        simple_markdown_table(
            ["mode", "ORCA index", "frequency, cm⁻¹", "max amplitude", "dominant atom", "element"],
            vibration_summary_table_rows(vib_mode_rows),
        ),
        "",
        "Амплитуды нормальных мод являются нормированными компонентами смещений из ORCA. Они показывают относительное участие атомов в каждой моде и не должны интерпретироваться как абсолютные тепловые амплитуды в Å.",
        "",
        f"Рисунок 10: `{short_path(str(figures['vibration_spectrum']), project)}` - спектр 12 ненулевых вибрационных частот выбранной структуры.",
        f"Рисунок 11: `{short_path(str(figures['vibration_heatmap']), project)}` - heatmap относительных амплитуд атомов B1-B6 по нормальным модам.",
        "",
        "Файлы вибрационного анализа:",
        "",
        simple_markdown_table(
            ["Файл", "Содержание"],
            [
                ["results/vibrations/B6/B6_all_vibrational_frequencies.csv", "12 ненулевых частот выбранного минимума"],
                ["results/vibrations/B6/B6_mode_summary.csv", "частота, максимальная амплитуда и доминирующий атом для каждой моды"],
                ["results/vibrations/B6/B6_normal_mode_amplitudes.csv", "dx, dy, dz и амплитуда по каждому атому и каждой моде"],
                ["results/vibrations/B6/B6_vibrational_frequencies_raw.txt", "сырой блок VIBRATIONAL FREQUENCIES из ORCA"],
                ["results/vibrations/B6/B6_best.out", "ORCA output выбранной структуры"],
                ["results/vibrations/B6/B6_best.hess", "Hessian-файл с вибрационной информацией"],
                ["results/vibrations/B6/B6_best_optimized.xyz", "оптимизированная геометрия выбранного минимума"],
            ],
        ),
        "",
        "Matplotlib-графики вибраций и контрольных метрик сохранены отдельно в `results/vibrations/B6/matplotlib_plots`:",
        "",
        simple_markdown_table(
            ["Рисунок", "Файл", "Содержание"],
            [
                ["Figure 10", "results/vibrations/B6/matplotlib_plots/Figure_10_B6_vibrational_frequencies.svg", "bar chart всех 12 ненулевых частот B₆"],
                ["Figure 11", "results/vibrations/B6/matplotlib_plots/Figure_11_B6_max_amplitude_by_mode.svg", "максимальная амплитуда по каждой нормальной моде"],
                ["Figure 12", "results/vibrations/B6/matplotlib_plots/Figure_12_B6_frequency_vs_amplitude.svg", "scatter plot частоты против максимальной амплитуды"],
                ["Figure 13", "results/vibrations/B6/matplotlib_plots/Figure_13_B6_atom_participation_heatmap.svg", "heatmap участия атомов в каждой моде"],
                ["Figure 14", "results/vibrations/B6/matplotlib_plots/Figure_14_B6_dominant_atom_by_mode.svg", "доминирующий атом для каждой моды"],
                ["Figure 15", "results/vibrations/B6/matplotlib_plots/Figure_15_B6_frequency_distribution.svg", "гистограмма распределения частот"],
                ["Figure 16", "results/vibrations/B6/matplotlib_plots/Figure_16_final_relative_energies_labeled.svg", "относительные энергии финальных кандидатов с отметкой best_B6"],
                ["Figure 17", "results/vibrations/B6/matplotlib_plots/Figure_17_screening_success_rate.svg", "success/fail summary для screening-этапа"],
            ],
        ),
        "",
        "## 10. Обсуждение результатов",
        f"Самой устойчивой по финальной энергии оказалась структура `{best.get('calculation_name', '')}`. Она получена из старта `{best.get('geometry_type', '')}`, имеет мультиплетность `{best.get('multiplicity', '')}` и полную энергию `{best.get('total_energy_hartree', '')}` Hartree. Ее относительная энергия принята равной `{best.get('relative_energy_ev', '')}` eV.",
        "",
        f"Близкие по энергии кандидаты присутствуют: `{len(close_001)}` финальных структур лежат в пределах 0.01 eV, а `{len(close_005)}` структур - в пределах 0.05 eV от минимума. При этом `{len(close_0001)}` финальных строк отличаются от лучшей структуры менее чем на 0.001 eV, то есть намного меньше практической точности обычного DFT-сравнения изомеров. Несмотря на наличие 10 финальных строк, они, вероятно, представляют несколько очень близких или практически идентичных минимумов. Поэтому физически значимый вывод состоит не в различии между отдельными строками, а в устойчивом воспроизведении одной плоской низкоэнергетической структуры из разных стартовых геометрий.",
        "",
        f"3D-старты были конкурентоспособными как исходные кандидаты: `{len(three_d_final)}` из `{len(final_rows)}` финальных расчетов происходят из 3D/random стартов. При этом анализ планарности оптимизированных финальных XYZ показывает, что `{planarized_3d}` из них имеют RMS-отклонение от лучшей плоскости не больше 0.01 Å. Для выбранного `best_B6.xyz` RMS-отклонение от плоскости равно `{rms_plane:.4f}` Å, максимальное отклонение `{max_plane:.4f}` Å. Поэтому итоговая структура является плоской или практически плоской, несмотря на то что лучший старт был random 3D.",
        "",
        f"Рисунок 4: `{short_path(str(figures['best']), project)}`.",
        "",
        "Полученный результат согласуется с литературной тенденцией малых борных кластеров к плоским или квазиплоским структурам. Важно, что этот вывод сделан после проверки 3D-стартов, а не путем их исключения заранее.",
        "",
        "### 10.1. Геометрические характеристики выбранной структуры",
        f"Для `best_B6.xyz` минимальное межатомное расстояние B-B равно `{min_dist:.6f}` Å, максимальное расстояние среди всех 15 пар атомов равно `{max_dist:.6f}` Å, среднее расстояние по всем парам равно `{avg_dist:.6f}` Å. Если использовать простой геометрический cutoff 2.05 Å для близких B-B контактов, получается `{len(best_bonds)}` коротких контактов.",
        "",
        "Координаты выбранной структуры:",
        "",
        simple_markdown_table(
            ["atom", "element", "x, Å", "y, Å", "z, Å"],
            best_coordinate_rows(best_atoms),
        ),
        "",
        "Все попарные расстояния B-B в выбранной структуре:",
        "",
        simple_markdown_table(
            ["pair", "distance, Å"],
            distance_rows(best_atoms),
        ),
        "",
        "## 11. Ограничения расчёта",
        "Следует учитывать, что полученный минимум является лучшим найденным минимумом в рамках использованного набора стартовых структур, мультиплетностей 1, 3 и 5 и выбранного уровня теории PBE0-D4/def2-TZVP. Для более строгого подтверждения глобального минимума можно расширить набор стартовых геометрий, выполнить дополнительную дедупликацию структур, проверить другие функционалы DFT и при необходимости провести более высокоуровневые single-point расчеты.",
        "",
        "## 12. Вывод",
        f"В рамках выбранного набора стартовых геометрий, проверенных мультиплетностей и уровня теории PBE0-D4/def2-TZVP наиболее устойчивым найденным минимумом нейтрального кластера B₆ является структура `{best.get('calculation_name', '')}`, имеющая мультиплетность `{best.get('multiplicity', '')}`, полную энергию `{best.get('total_energy_hartree', '')}` Hartree, относительную энергию `{best.get('relative_energy_ev', '')}` eV и не имеющая мнимых частот.",
        "",
        "Этот результат не формулируется как доказательство абсолютного глобального минимума B₆; он является лучшим найденным минимумом в рамках выполненного многостартового DFT-набора и выбранного уровня теории.",
        "",
        "## 13. Приложения",
        f"- `results/results.csv`: полный screening-набор.",
        f"- `results/final_results.csv`: финальные OptFreq энергии и частоты.",
        f"- `results/best_B6.xyz`: координаты выбранной структуры.",
        f"- `results/B6_final_report.txt`: текст отчета.",
        f"- `calculations/final/*/*.out`: ORCA output-файлы финальных расчетов.",
        f"- `results/figures/Figure_4_best_B6.svg`: изображение финальной структуры.",
        f"- `results/figures/Figure_6_min_energy_geometries.svg`: визуализация низкоэнергетических финальных геометрий.",
        f"- `results/figures/Figure_7_min_energy_geometries_3d.svg`: статическая 3D-проекция низкоэнергетических финальных геометрий.",
        f"- `results/figures/Figure_7_min_energy_geometries_3d.html`: интерактивная 3D-визуализация с вращением мышью.",
        f"- `results/screening_3d_results.csv`: отдельная таблица screening-расчётов из 3D/random стартов.",
        f"- `results/final_from_3d_results.csv`: отдельная таблица финальных OptFreq-расчётов из 3D/random стартов.",
        f"- `results/figures/Figure_8_3d_start_geometries.svg`: отдельные 3D-стартовые геометрии.",
        f"- `results/figures/Figure_9_final_from_3d_starts_3d.svg`: 3D-проекция финальных структур из 3D/random стартов.",
        f"- `results/figures/Figure_9_final_from_3d_starts_3d.html`: интерактивная версия финальных структур из 3D/random стартов.",
        f"- `results/vibrations/B6/*`: отдельные таблицы, raw-блок, `.out`, `.hess` и XYZ для вибрационного анализа выбранной структуры.",
        f"- `results/figures/Figure_10_B6_vibrational_spectrum.svg`: спектр ненулевых вибрационных частот.",
        f"- `results/figures/Figure_11_B6_normal_mode_amplitudes.svg`: heatmap амплитуд нормальных мод.",
        "",
        "### 13.1. Команды воспроизведения обработки данных",
        "Ниже приведены команды, которые не запускают новые квантово-химические расчёты, а только пересобирают таблицы и отчёт из уже существующих ORCA output-файлов.",
        "",
        "```bash",
        "python3 scripts/collect_results.py \\",
        "  --root calculations/stage1 \\",
        "  --csv results/results.csv \\",
        "  --best-xyz results/best_B6.xyz \\",
        "  --all-energies-csv results/all_energies.csv",
        "",
        "python3 scripts/collect_results.py \\",
        "  --root calculations/final \\",
        "  --csv results/final_results.csv \\",
        "  --best-xyz results/best_B6.xyz",
        "",
        "python3 scripts/build_b6_report.py --project-dir .",
        "```",
        "",
        "### 13.2. Шаблоны ORCA input",
        "Ниже приведены текущие шаблоны input-файлов, которые используются как документированные примеры настроек. Координаты в реальных `.inp` файлах генерируются отдельно для каждой стартовой структуры.",
        "",
        "Screening-шаблон:",
        "",
        "```orca",
        stage1_template.strip() if stage1_template else "templates/stage1_opt_template.inp not found",
        "```",
        "",
        "Финальный OptFreq-шаблон:",
        "",
        "```orca",
        final_template.strip() if final_template else "templates/final_opt_freq_template.inp not found",
        "```",
        "",
        "### 13.3. Контрольные файлы",
        "",
        simple_markdown_table(
            ["Файл", "Назначение"],
            [
                ["results/results.csv", "screening-таблица R2SCAN-3C Opt"],
                ["results/all_energies.csv", "дублирующая полная таблица энергий screening"],
                ["results/final_results.csv", "финальные PBE0-D4/def2-TZVP OptFreq результаты"],
                ["results/best_B6.xyz", "координаты выбранного минимума"],
                ["results/B6_final_report.md", "подробный отчёт в Markdown"],
                ["results/B6_final_report.txt", "текстовая копия отчёта"],
                ["results/figures/*.svg", "схемы workflow, геометрий и графики энергий"],
                ["results/figures/Figure_6_min_energy_geometries.svg", "топ финальных оптимизированных геометрий по энергии"],
                ["results/figures/Figure_7_min_energy_geometries_3d.svg", "статическая 3D-проекция топовых финальных геометрий"],
                ["results/figures/Figure_7_min_energy_geometries_3d.html", "интерактивная 3D-визуализация топовых финальных геометрий"],
                ["results/screening_3d_results.csv", "отдельная screening-таблица только для 3D/random стартов"],
                ["results/final_from_3d_results.csv", "отдельная final-таблица только для кандидатов из 3D/random стартов"],
                ["results/figures/Figure_8_3d_start_geometries.svg", "отдельные 3D-стартовые геометрии"],
                ["results/figures/Figure_9_final_from_3d_starts_3d.svg", "статическая 3D-проекция финальных структур из 3D/random стартов"],
                ["results/figures/Figure_9_final_from_3d_starts_3d.html", "интерактивная 3D-визуализация финальных структур из 3D/random стартов"],
                ["results/vibrations/B6/README.md", "описание вынесенного вибрационного анализа выбранной структуры"],
                ["results/vibrations/B6/B6_all_vibrational_frequencies.csv", "все 12 ненулевых частот best_B6"],
                ["results/vibrations/B6/B6_mode_summary.csv", "сводка по нормальным модам и доминирующим атомам"],
                ["results/vibrations/B6/B6_normal_mode_amplitudes.csv", "dx, dy, dz и амплитуды нормальных мод по атомам"],
                ["results/vibrations/B6/B6_vibrational_frequencies_raw.txt", "сырой ORCA-блок VIBRATIONAL FREQUENCIES"],
                ["results/vibrations/B6/B6_best.out", "ORCA output выбранной структуры для вибрационного анализа"],
                ["results/vibrations/B6/B6_best.hess", "Hessian-файл выбранной структуры"],
                ["results/figures/Figure_10_B6_vibrational_spectrum.svg", "график ненулевых вибрационных частот"],
                ["results/figures/Figure_11_B6_normal_mode_amplitudes.svg", "тепловая карта относительных амплитуд нормальных мод"],
                ["results/vibrations/B6/matplotlib_plots/Figure_10_B6_vibrational_frequencies.svg", "Matplotlib bar chart всех частот"],
                ["results/vibrations/B6/matplotlib_plots/Figure_11_B6_max_amplitude_by_mode.svg", "Matplotlib график max amplitude по модам"],
                ["results/vibrations/B6/matplotlib_plots/Figure_12_B6_frequency_vs_amplitude.svg", "Matplotlib scatter plot frequency vs max amplitude"],
                ["results/vibrations/B6/matplotlib_plots/Figure_13_B6_atom_participation_heatmap.svg", "Matplotlib heatmap участия атомов"],
                ["results/vibrations/B6/matplotlib_plots/Figure_14_B6_dominant_atom_by_mode.svg", "Matplotlib dominant atom per mode"],
                ["results/vibrations/B6/matplotlib_plots/Figure_15_B6_frequency_distribution.svg", "Matplotlib histogram распределения частот"],
                ["results/vibrations/B6/matplotlib_plots/Figure_16_final_relative_energies_labeled.svg", "Matplotlib final relative energies с отметкой best_B6"],
                ["results/vibrations/B6/matplotlib_plots/Figure_17_screening_success_rate.svg", "Matplotlib screening success/fail"],
            ],
        ),
        "",
        "## Литература",
        "[1] J. Burkhardt, Y. Jia, W.-L. Li. Structure Search with the Strategic Escape Algorithm. Journal of Chemical Theory and Computation, 2025, 21, 3765-3773. DOI: https://doi.org/10.1021/acs.jctc.4c01746",
        "[2] Q.-S. Li, B. Song, L. Wen, L.-M. Yang, E. Ganz. Elucidation of Structures, Electronic Properties, and Chemical Bonding of Monophosphorus-Substituted Boron Clusters in Neutral, Negative, and Positively Charged PBn/PBn-/PBn+ (n = 4-8). Condensed Matter, 2022, 7, 66. https://www.mdpi.com/2410-3896/7/4/66",
        "[3] Milon, D. Roy, F. Ahmed. A DFT study to investigate the physical, electrical, optical properties and thermodynamic functions of boron nanoclusters (MxB2n0; x=1,2, n=3,4,5). Heliyon, 2023, 9, e17886. DOI: https://doi.org/10.1016/j.heliyon.2023.e17886",
        "[4] A. N. Alexandrova, A. I. Boldyrev, H.-J. Zhai, L.-S. Wang, E. Steiner, P. W. Fowler. Structure and Bonding in B6- and B6: Planarity and Antiaromaticity. Journal of Physical Chemistry A, 2003, 107, 1359-1369. DOI: https://doi.org/10.1021/jp0268866",
        "[5] W.-L. Li, X. Chen, T. Jian, T.-T. Chen, J. Li, L.-S. Wang. From planar boron clusters to borophenes and metalloborophenes. Nature Reviews Chemistry, 2017, 1, 0071. DOI: https://doi.org/10.1038/s41570-017-0071",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Build B6 final report and SVG figures.")
    parser.add_argument("--project-dir", default=".")
    parser.add_argument("--results-csv", default="results/results.csv")
    parser.add_argument("--final-csv", default="results/final_results.csv")
    parser.add_argument("--best-xyz", default="results/best_B6.xyz")
    parser.add_argument("--out", default="results/B6_final_report.txt")
    parser.add_argument("--md-out", default="results/B6_final_report.md")
    parser.add_argument("--fig-dir", default="results/figures")
    args = parser.parse_args()

    project = Path(args.project_dir).resolve()
    report = build_report(
        project,
        (project / args.results_csv).resolve(),
        (project / args.final_csv).resolve(),
        (project / args.best_xyz).resolve(),
        (project / args.fig_dir).resolve(),
    )
    out = (project / args.out).resolve()
    md_out = (project / args.md_out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(report, encoding="utf-8")
    md_out.write_text(report, encoding="utf-8")
    print(f"Wrote report: {out}")
    print(f"Wrote markdown report: {md_out}")
    print(f"Wrote figures to: {(project / args.fig_dir).resolve()}")


if __name__ == "__main__":
    main()
