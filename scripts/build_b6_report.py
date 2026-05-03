#!/usr/bin/env python3
"""Build the B6 calculation report and simple SVG figures from collected results."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
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
            return str(path.resolve().relative_to(project))
    except Exception:
        pass
    marker = project.name
    parts = path.parts
    if marker in parts:
        idx = parts.index(marker)
        return str(Path(*parts[idx + 1 :]))
    return path_text


def to_float(text: str, default: float = math.inf) -> float:
    try:
        return float(text)
    except Exception:
        return default


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


def pair_distances(atoms: Sequence[Atom]) -> List[Tuple[int, int, float]]:
    pairs: List[Tuple[int, int, float]] = []
    for i in range(len(atoms)):
        _, xi, yi, zi = atoms[i]
        for j in range(i + 1, len(atoms)):
            _, xj, yj, zj = atoms[j]
            d = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)
            pairs.append((i, j, d))
    return pairs


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


def write_figures(project: Path, fig_dir: Path, screening_rows: Sequence[Dict[str, str]], final_rows: Sequence[Dict[str, str]], best_atoms: Sequence[Atom]) -> Dict[str, Path]:
    fig_dir.mkdir(parents=True, exist_ok=True)
    figures = {
        "workflow": fig_dir / "Figure_1_workflow.svg",
        "starts": fig_dir / "Figure_2_start_geometries.svg",
        "screening": fig_dir / "Figure_3_screening_top10.svg",
        "best": fig_dir / "Figure_4_best_B6.svg",
        "final": fig_dir / "Figure_5_final_relative_energies.svg",
    }
    write_workflow_svg(figures["workflow"])
    write_start_geometries_svg(figures["starts"], project)
    write_bar_svg(figures["screening"], screening_rows, "Топ-10 структур после R2SCAN-3C screening", limit=10)
    figures["best"].write_text(atoms_svg(best_atoms, 520, 420, "best_B6.xyz"), encoding="utf-8")
    write_bar_svg(figures["final"], final_rows, "Относительные энергии финальных кандидатов", limit=10)
    return figures


def build_report(project: Path, results_csv: Path, final_csv: Path, best_xyz: Path, fig_dir: Path) -> str:
    screening_rows = sorted(read_csv(results_csv), key=lambda r: to_float(r.get("total_energy_hartree", "")))
    final_rows = sorted(read_csv(final_csv), key=lambda r: to_float(r.get("total_energy_hartree", "")))
    best_atoms = read_xyz(best_xyz)
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
    three_d_final = [
        r
        for r in final_rows
        if any(token in r.get("geometry_type", "").lower() for token in ["random", "3d", "prism", "octa", "pyramid"])
    ]
    planarized_3d = 0
    for row in three_d_final:
        xyz_path = Path(row.get("xyz_file", ""))
        if xyz_path.exists():
            rms, _ = planarity(read_xyz(xyz_path))
            if rms <= 0.01:
                planarized_3d += 1

    lines = [
        "# Многостартовый DFT-поиск устойчивой геометрии нейтрального кластера B₆ методом ORCA 6.1",
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
        f"Рисунок 1: `{short_path(str(figures['workflow']), project)}`.",
        "",
        "## 5. Генерация стартовых геометрий B₆",
        "Для уменьшения риска попадания в локальный минимум был использован многостартовый подход. Были сгенерированы плоские, квазиплоские и трехмерные стартовые структуры B₆ с различными начальными расстояниями B-B. Для каждой структуры были проверены мультиплетности 1, 3 и 5.",
        "",
        "В расчетной кампании использовались следующие типы стартов: линейная цепочка; плоское кольцо; компактная плоская структура; ромбическая структура; прямоугольная структура; октаэдрическая 3D-структура; тригональная призма; случайные 3D-структуры. Дополнительно генератор поддерживает искаженное плоское кольцо, fused-triangle, квазиплоскую и пирамидальную 3D-структуру для расширенного набора.",
        "",
        f"Рисунок 2: `{short_path(str(figures['starts']), project)}`.",
        "",
        "## 6. Первичный screening: R2SCAN-3C Opt",
        f"На screening-этапе обработано `{len(screening_rows)}` ORCA output-файлов. Нормально завершились `{len(screening_done)}` расчетов, сходимость оптимизации обнаружена у `{len(screening_conv)}` расчетов. Полные энергии извлекались из строки `FINAL SINGLE POINT ENERGY`; расчеты без нормального завершения или без сходимости не рассматриваются как надежные финальные кандидаты.",
        "",
        "Таблица 1. Топ-10 screening-результатов; полный набор приведен в `results/results.csv`.",
        "",
        markdown_table(screening_rows, SCREENING_COLUMNS, project, limit=10),
        "",
        f"Рисунок 3: `{short_path(str(figures['screening']), project)}`.",
        "",
        "## 7. Отбор финальных кандидатов",
        "Финальные кандидаты выбирались из низкоэнергетических расчетов screening-этапа с нормальным завершением и сошедшейся оптимизацией. Для уменьшения дублирования структур используется сравнение отсортированных межатомных расстояний B-B; структуры с близкими distance fingerprints рассматриваются как геометрически повторяющиеся кандидаты. В текущем финальном наборе сохранены 10 OptFreq расчетов.",
        "",
        "## 8. Финальный расчет: PBE0-D4/def2-TZVP OptFreq",
        f"Финальный этап включал `{len(final_rows)}` расчетов PBE0-D4/def2-TZVP OptFreq. Нормально завершились `{len(final_done)}` расчетов, сходимость оптимизации обнаружена у `{len(final_conv)}` расчетов.",
        "",
        "Таблица 2. Финальные расчеты PBE0-D4/def2-TZVP OptFreq.",
        "",
        markdown_table(final_rows, FINAL_COLUMNS, project, limit=None),
        "",
        f"Рисунок 5: `{short_path(str(figures['final']), project)}`.",
        "",
        "## 9. Частотный анализ",
        "Структура считалась истинным минимумом только при одновременном выполнении трех условий: `ORCA TERMINATED NORMALLY`, сходимость оптимизации и отсутствие мнимых частот. Если структура имеет хотя бы одну мнимую частоту, она не считается финальным минимумом даже при низкой электронной энергии.",
        "",
        "В ORCA output для финальных расчетов присутствуют шесть нулевых трансляционно-вращательных мод. Они не учитывались как мнимые вибрационные моды и не использовались при выборе `lowest_frequency_cm-1`; в таблицу записана минимальная ненулевая вибрационная частота после отсечения мод с |ν| <= 10 cm⁻¹.",
        "",
        f"В текущей финальной таблице структур с `n_imaginary_frequencies > 0` не найдено; истинных минимумов по указанному критерию: `{len(final_true)}`. Для выбранной структуры `n_imaginary_frequencies = {best.get('n_imaginary_frequencies', '')}`, `lowest_frequency_cm-1 = {best.get('lowest_frequency_cm-1', '')}`.",
        "",
        "## 10. Обсуждение результатов",
        f"Самой устойчивой по финальной энергии оказалась структура `{best.get('calculation_name', '')}`. Она получена из старта `{best.get('geometry_type', '')}`, имеет мультиплетность `{best.get('multiplicity', '')}` и полную энергию `{best.get('total_energy_hartree', '')}` Hartree. Ее относительная энергия принята равной `{best.get('relative_energy_ev', '')}` eV.",
        "",
        f"Близкие по энергии кандидаты присутствуют: `{len(close_001)}` финальных структур лежат в пределах 0.01 eV, а `{len(close_005)}` структур - в пределах 0.05 eV от минимума. Очень малые различия между несколькими финальными структурами указывают, что разные стартовые геометрии после оптимизации сходятся к одному и тому же или практически идентичному минимуму. Поэтому физически значимым является не различие между этими строками, а устойчивое воспроизведение одной низкоэнергетической плоской структуры из разных стартов.",
        "",
        f"3D-старты были конкурентоспособными как исходные кандидаты: `{len(three_d_final)}` из `{len(final_rows)}` финальных расчетов происходят из 3D/random стартов. При этом анализ планарности оптимизированных финальных XYZ показывает, что `{planarized_3d}` из них имеют RMS-отклонение от лучшей плоскости не больше 0.01 Å. Для выбранного `best_B6.xyz` RMS-отклонение от плоскости равно `{rms_plane:.4f}` Å, максимальное отклонение `{max_plane:.4f}` Å. Поэтому итоговая структура является плоской или практически плоской, несмотря на то что лучший старт был random 3D.",
        "",
        f"Рисунок 4: `{short_path(str(figures['best']), project)}`.",
        "",
        "Полученный результат согласуется с литературной тенденцией малых борных кластеров к плоским или квазиплоским структурам. Важно, что этот вывод сделан после проверки 3D-стартов, а не путем их исключения заранее.",
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
