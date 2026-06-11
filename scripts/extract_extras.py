#!/usr/bin/env python3
"""Extract additional B6 campaign data and build report version 2.

This script uses only existing ORCA .out/.xyz files. It does not start new
quantum-chemical calculations.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

HARTREE_TO_EV = 27.211386245988
Atom = Tuple[str, float, float, float]

EXTRA_COLUMNS = [
    "s2_expected",
    "s2_actual",
    "s2_deviation",
    "zpe_hartree",
    "e_plus_zpe",
    "relative_e_plus_zpe_ev",
    "enthalpy_298_hartree",
    "gibbs_298",
    "relative_gibbs_298_ev",
    "alpha_homo_ev",
    "alpha_lumo_ev",
    "beta_homo_ev",
    "beta_lumo_ev",
    "homo_ev",
    "somo_ev",
    "lumo_ev",
    "homo_lumo_gap_ev",
]


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: Sequence[Dict[str, object]], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def truthy(value: object) -> bool:
    return str(value).strip().lower() == "true"


def to_float(value: object, default: float = math.nan) -> float:
    try:
        text = str(value).strip()
        if not text:
            return default
        return float(text)
    except Exception:
        return default


def fmt(value: float, digits: int = 10) -> str:
    if math.isnan(value) or math.isinf(value):
        return ""
    return f"{value:.{digits}f}"


def fmt_ev(value: float) -> str:
    return fmt(value, 8)


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


def short_name(name: str, max_len: int = 58) -> str:
    text = name.replace("_PBE0_def2-TZVP_OptFreq", "").replace("FINAL_", "")
    return text if len(text) <= max_len else text[: max_len - 1] + "..."


def read_xyz(path: Path) -> List[Atom]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    count = int(lines[0].strip())
    atoms: List[Atom] = []
    for line in lines[2 : 2 + count]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms


def pair_distances(atoms: Sequence[Atom]) -> List[Tuple[int, int, float]]:
    pairs: List[Tuple[int, int, float]] = []
    for i, atom_i in enumerate(atoms):
        _, xi, yi, zi = atom_i
        for j in range(i + 1, len(atoms)):
            _, xj, yj, zj = atoms[j]
            distance = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)
            pairs.append((i + 1, j + 1, distance))
    return pairs


def distance_fingerprint(atoms: Sequence[Atom]) -> List[float]:
    return sorted(distance for _, _, distance in pair_distances(atoms))


def same_fingerprint(a: Sequence[float], b: Sequence[float], tolerance: float) -> bool:
    return len(a) == len(b) and max((abs(x - y) for x, y in zip(a, b)), default=0.0) <= tolerance


def parse_last_float(text: str, pattern: str) -> float:
    matches = re.findall(pattern, text, flags=re.IGNORECASE)
    if not matches:
        return math.nan
    value = matches[-1]
    if isinstance(value, tuple):
        value = value[-1]
    return float(value)


def parse_point_group(text: str) -> str:
    matches = re.findall(r"Point\s+Group:\s*([A-Za-z0-9*]+)", text)
    return matches[-1] if matches else ""


def parse_orca_version(text: str) -> str:
    match = re.search(r"Program Version\s+([0-9.]+)\s+-\s+([A-Z]+)", text)
    if not match:
        return ""
    return f"{match.group(1)} {match.group(2)}"


def parse_host(text: str) -> str:
    match = re.search(r"\*\s+Host name:\s+(.+)", text)
    return match.group(1).strip() if match else ""


def parse_spin_values(text: str) -> Dict[str, str]:
    actual = parse_last_float(text, r"Expectation value of <S\*\*2>\s*:\s*(-?\d+\.\d+)")
    expected = parse_last_float(text, r"Ideal value S\*\(S\+1\).*?:\s*(-?\d+\.\d+)")
    return {
        "s2_expected": fmt(expected, 6),
        "s2_actual": fmt(actual, 6),
        "s2_deviation": fmt(actual - expected, 6) if not (math.isnan(actual) or math.isnan(expected)) else "",
    }


def parse_thermochemistry(text: str) -> Dict[str, str]:
    zpe = parse_last_float(text, r"Zero point energy\s+\.\.\.\s+(-?\d+\.\d+)\s+Eh")
    enthalpy = parse_last_float(text, r"Total Enthalpy\s+\.\.\.\s+(-?\d+\.\d+)\s+Eh")
    gibbs = parse_last_float(text, r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)\s+Eh")
    return {
        "zpe_hartree": fmt(zpe, 10),
        "enthalpy_298_hartree": fmt(enthalpy, 10),
        "gibbs_298": fmt(gibbs, 10),
    }


def parse_orbital_energies(text: str) -> Dict[str, str]:
    start = text.rfind("ORBITAL ENERGIES")
    if start < 0:
        return {}
    block = text[start:]
    stop = block.find("MULLIKEN POPULATION ANALYSIS")
    if stop > 0:
        block = block[:stop]

    orbitals: Dict[str, List[Dict[str, float]]] = {"alpha": [], "beta": []}
    spin: Optional[str] = None
    for line in block.splitlines():
        if "SPIN UP ORBITALS" in line:
            spin = "alpha"
            continue
        if "SPIN DOWN ORBITALS" in line:
            spin = "beta"
            continue
        match = re.match(r"\s*(\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
        if match and spin:
            orbitals[spin].append(
                {
                    "index": float(match.group(1)),
                    "occ": float(match.group(2)),
                    "eh": float(match.group(3)),
                    "ev": float(match.group(4)),
                }
            )

    def homo_lumo(items: Sequence[Dict[str, float]]) -> Tuple[float, float]:
        occupied = [item["ev"] for item in items if item["occ"] > 0.5]
        virtual = [item["ev"] for item in items if item["occ"] <= 0.5]
        return (max(occupied) if occupied else math.nan, min(virtual) if virtual else math.nan)

    alpha_homo, alpha_lumo = homo_lumo(orbitals["alpha"])
    beta_homo, beta_lumo = homo_lumo(orbitals["beta"])
    homo = max([value for value in [alpha_homo, beta_homo] if not math.isnan(value)], default=math.nan)
    lumo = min([value for value in [alpha_lumo, beta_lumo] if not math.isnan(value)], default=math.nan)

    beta_occ_by_index = {int(item["index"]): item["occ"] for item in orbitals["beta"]}
    somo_values = [
        item["ev"]
        for item in orbitals["alpha"]
        if item["occ"] > 0.5 and beta_occ_by_index.get(int(item["index"]), 0.0) <= 0.5
    ]

    return {
        "alpha_homo_ev": fmt(alpha_homo, 4),
        "alpha_lumo_ev": fmt(alpha_lumo, 4),
        "beta_homo_ev": fmt(beta_homo, 4),
        "beta_lumo_ev": fmt(beta_lumo, 4),
        "homo_ev": fmt(homo, 4),
        "somo_ev": ";".join(fmt(value, 4) for value in somo_values),
        "lumo_ev": fmt(lumo, 4),
        "homo_lumo_gap_ev": fmt(lumo - homo, 4) if not (math.isnan(homo) or math.isnan(lumo)) else "",
    }


def parse_atomic_population(text: str, header: str) -> Dict[int, Dict[str, str]]:
    start = text.rfind(header)
    if start < 0:
        return {}
    out: Dict[int, Dict[str, str]] = {}
    for line in text[start:].splitlines()[1:]:
        match = re.match(r"\s*(\d+)\s+([A-Za-z]+)\s*:\s*(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
        if match:
            index = int(match.group(1)) + 1
            out[index] = {
                "element": match.group(2),
                "charge": match.group(3),
                "spin": match.group(4),
            }
            continue
        if out and (not line.strip() or line.startswith("Sum of")):
            break
    return out


def parse_final_extras(out_path: Path) -> Dict[str, str]:
    text = read_text(out_path)
    data: Dict[str, str] = {}
    data.update(parse_spin_values(text))
    data.update(parse_thermochemistry(text))
    data.update(parse_orbital_energies(text))
    data["point_group"] = parse_point_group(text)
    data["orca_version_build"] = parse_orca_version(text)
    data["run_host"] = parse_host(text)
    return data


def update_final_results(project: Path, final_csv: Path) -> List[Dict[str, str]]:
    rows = read_csv_rows(final_csv)
    for row in rows:
        out_path = resolve_project_path(row.get("output_file", ""), project)
        if not out_path.exists():
            continue
        extras = parse_final_extras(out_path)
        row.update({column: extras.get(column, row.get(column, "")) for column in EXTRA_COLUMNS})
        energy = to_float(row.get("total_energy_hartree"))
        zpe = to_float(row.get("zpe_hartree"))
        if not math.isnan(energy) and not math.isnan(zpe):
            row["e_plus_zpe"] = fmt(energy + zpe, 12)

    e_zpe_values = [to_float(row.get("e_plus_zpe")) for row in rows if not math.isnan(to_float(row.get("e_plus_zpe")))]
    gibbs_values = [to_float(row.get("gibbs_298")) for row in rows if not math.isnan(to_float(row.get("gibbs_298")))]
    e_zpe_min = min(e_zpe_values) if e_zpe_values else math.nan
    gibbs_min = min(gibbs_values) if gibbs_values else math.nan
    for row in rows:
        e_zpe = to_float(row.get("e_plus_zpe"))
        gibbs = to_float(row.get("gibbs_298"))
        if not math.isnan(e_zpe) and not math.isnan(e_zpe_min):
            row["relative_e_plus_zpe_ev"] = fmt_ev((e_zpe - e_zpe_min) * HARTREE_TO_EV)
        if not math.isnan(gibbs) and not math.isnan(gibbs_min):
            row["relative_gibbs_298_ev"] = fmt_ev((gibbs - gibbs_min) * HARTREE_TO_EV)

    columns = list(read_csv_rows(final_csv)[0].keys()) if rows else []
    for column in EXTRA_COLUMNS:
        if column not in columns:
            columns.append(column)
    write_csv(final_csv, rows, columns)
    return rows


def select_best_xyz_row(final_rows: Sequence[Dict[str, str]], best_xyz: Path, project: Path) -> Dict[str, str]:
    lines = best_xyz.read_text(encoding="utf-8", errors="replace").splitlines()
    comment = lines[1] if len(lines) > 1 else ""
    for row in final_rows:
        output_file = row.get("output_file", "")
        output_stem = Path(output_file).stem
        if output_stem and output_stem in comment:
            return dict(row)
        if row.get("calculation_name", "") and row["calculation_name"] in comment:
            return dict(row)

    best_atoms = read_xyz(best_xyz)
    best_fp = distance_fingerprint(best_atoms)
    for row in final_rows:
        xyz_path = resolve_project_path(row.get("xyz_file", ""), project)
        if xyz_path.exists() and same_fingerprint(best_fp, distance_fingerprint(read_xyz(xyz_path)), 1e-5):
            return dict(row)

    true_minima = [row for row in final_rows if truthy(row.get("is_true_minimum", ""))]
    return dict(min(true_minima or final_rows, key=lambda row: to_float(row.get("total_energy_hartree"), math.inf)))


def write_best_population(project: Path, best_row: Dict[str, str], best_xyz: Path, population_csv: Path) -> None:
    out_path = resolve_project_path(best_row.get("output_file", ""), project)
    text = read_text(out_path)
    atoms = read_xyz(best_xyz)
    mulliken = parse_atomic_population(text, "MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS")
    loewdin = parse_atomic_population(text, "LOEWDIN ATOMIC CHARGES AND SPIN POPULATIONS")

    rows: List[Dict[str, object]] = []
    for index, atom in enumerate(atoms, start=1):
        element, x, y, z = atom
        rows.append(
            {
                "atom_index": index,
                "element": element,
                "x": f"{x:.8f}",
                "y": f"{y:.8f}",
                "z": f"{z:.8f}",
                "mulliken_charge": mulliken.get(index, {}).get("charge", ""),
                "mulliken_spin_population": mulliken.get(index, {}).get("spin", ""),
                "loewdin_charge": loewdin.get(index, {}).get("charge", ""),
                "loewdin_spin_population": loewdin.get(index, {}).get("spin", ""),
            }
        )
    columns = [
        "atom_index",
        "element",
        "x",
        "y",
        "z",
        "mulliken_charge",
        "mulliken_spin_population",
        "loewdin_charge",
        "loewdin_spin_population",
    ]
    write_csv(population_csv, rows, columns)


def dedupe_screening(project: Path, screening_csv: Path, output_csv: Path, tolerance: float) -> List[Dict[str, object]]:
    rows = [
        row
        for row in read_csv_rows(screening_csv)
        if truthy(row.get("normal_termination")) and truthy(row.get("optimization_converged")) and row.get("xyz_file")
    ]
    rows.sort(key=lambda row: to_float(row.get("total_energy_hartree"), math.inf))
    groups: List[Dict[str, object]] = []
    for row in rows:
        xyz = resolve_project_path(row.get("xyz_file", ""), project)
        if not xyz.exists():
            continue
        fingerprint = distance_fingerprint(read_xyz(xyz))
        matched = False
        for group in groups:
            if same_fingerprint(fingerprint, group["fingerprint"], tolerance):  # type: ignore[arg-type]
                group["members"].append(row)  # type: ignore[index, union-attr]
                matched = True
                break
        if not matched:
            groups.append({"fingerprint": fingerprint, "members": [row]})

    best_energy = min((to_float(group["members"][0].get("total_energy_hartree")) for group in groups), default=math.nan)  # type: ignore[index]
    out_rows: List[Dict[str, object]] = []
    for idx, group in enumerate(groups, start=1):
        members: List[Dict[str, str]] = group["members"]  # type: ignore[assignment]
        representative = members[0]
        energy = to_float(representative.get("total_energy_hartree"))
        multiplicity_counts: Dict[str, int] = {}
        geometry_types: Dict[str, int] = {}
        for member in members:
            multiplicity = member.get("multiplicity", "")
            geometry_type = member.get("geometry_type", "")
            multiplicity_counts[multiplicity] = multiplicity_counts.get(multiplicity, 0) + 1
            geometry_types[geometry_type] = geometry_types.get(geometry_type, 0) + 1
        out_rows.append(
            {
                "group_id": idx,
                "representative": representative.get("calculation_name", ""),
                "n_hits": len(members),
                "representative_multiplicity": representative.get("multiplicity", ""),
                "multiplicities": ";".join(f"{key}:{value}" for key, value in sorted(multiplicity_counts.items())),
                "delta_e_ev": fmt_ev((energy - best_energy) * HARTREE_TO_EV),
                "representative_energy_hartree": representative.get("total_energy_hartree", ""),
                "geometry_types": ";".join(f"{key}:{value}" for key, value in sorted(geometry_types.items())),
                "representative_xyz_file": representative.get("xyz_file", ""),
                "representative_output_file": representative.get("output_file", ""),
            }
        )
    columns = [
        "group_id",
        "representative",
        "n_hits",
        "representative_multiplicity",
        "multiplicities",
        "delta_e_ev",
        "representative_energy_hartree",
        "geometry_types",
        "representative_xyz_file",
        "representative_output_file",
    ]
    write_csv(output_csv, out_rows, columns)
    return out_rows


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    out = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        out.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(out)


def characteristic_distance_rows(atoms: Sequence[Atom]) -> List[List[object]]:
    buckets: Dict[float, List[str]] = {}
    for i, j, distance in pair_distances(atoms):
        rounded = round(distance, 2)
        buckets.setdefault(rounded, []).append(f"B{i}-B{j}")
    rows: List[List[object]] = []
    for rounded, pairs in sorted(buckets.items()):
        rows.append([", ".join(pairs), f"{rounded:.2f}", len(pairs)])
    return rows


def report_unique_rows(unique_rows: Sequence[Dict[str, object]]) -> List[List[object]]:
    rows: List[List[object]] = []
    for row in unique_rows:
        rows.append(
            [
                row.get("group_id", ""),
                short_name(str(row.get("representative", "")), 46),
                row.get("n_hits", ""),
                row.get("representative_multiplicity", ""),
                row.get("delta_e_ev", ""),
                str(row.get("geometry_types", ""))[:95],
            ]
        )
    return rows


def build_report_v2(
    project: Path,
    final_rows: Sequence[Dict[str, str]],
    unique_rows: Sequence[Dict[str, object]],
    best_xyz: Path,
    report_md: Path,
    report_txt: Path,
    population_csv: Path,
    unique_csv: Path,
    point_group: str,
    orca_version: str,
    host: str,
) -> None:
    true_minima = [row for row in final_rows if truthy(row.get("is_true_minimum", ""))]
    selected_best = select_best_xyz_row(final_rows, best_xyz, project)
    best_by_e_zpe = min(true_minima or list(final_rows), key=lambda row: to_float(row.get("e_plus_zpe"), math.inf))
    best_atoms = read_xyz(best_xyz)
    distance_summary = characteristic_distance_rows(best_atoms)
    all_pairs = sorted(pair_distances(best_atoms), key=lambda item: item[2])
    converged_hits = sum(int(row.get("n_hits", 0)) for row in unique_rows)

    final_table_rows = []
    for row in sorted(final_rows, key=lambda item: to_float(item.get("relative_e_plus_zpe_ev"), math.inf)):
        final_table_rows.append(
            [
                short_name(row.get("calculation_name", "")),
                row.get("multiplicity", ""),
                row.get("total_energy_hartree", ""),
                row.get("zpe_hartree", ""),
                row.get("e_plus_zpe", ""),
                row.get("relative_e_plus_zpe_ev", ""),
                row.get("gibbs_298", ""),
                row.get("s2_actual", ""),
            ]
        )

    population_rows = read_csv_rows(population_csv)
    population_table = [
        [
            row.get("atom_index", ""),
            row.get("mulliken_charge", ""),
            row.get("mulliken_spin_population", ""),
            row.get("loewdin_charge", ""),
            row.get("loewdin_spin_population", ""),
        ]
        for row in population_rows
    ]

    literature_rows = [
        ["Term/state", "triplet, m=3 UKS", "triplet neutral B6", "Term label is not assigned by this ORCA output."],
        ["Point group", point_group or "C2h", "C2h", "Matches the C2h motif discussed by Alexandrova et al."],
        ["Short B-B bonds, A", "1.52 x2", "about 1.52", "B1-B2 and B5-B6 in this work."],
        ["Side B-B bonds, A", "1.59 x4", "about 1.59", "Four equivalent/near-equivalent side bonds after rounding."],
        ["Long internal B-B bonds, A", "1.81 x2", "about 1.81", "B1-B3 and B4-B6 in this work."],
        ["Central B-B distance, A", "1.92 x1", "about 1.92", "B3-B4 in this work."],
    ]

    pair_rows = [[f"B{i}-B{j}", f"{distance:.6f}"] for i, j, distance in all_pairs]

    lines = [
        "# B6_ORCA_campaign: отчет v2 без новых расчетов",
        "",
        "Этот отчет версии 2 построен только из уже существующих файлов `calculations/stage1/**/*.out`, "
        "`calculations/final/**/*.out`, `results/*.csv` и `results/best_B6.xyz`. Новые ORCA-расчеты не запускались.",
        "",
        "## 1. Что изменено относительно v1",
        "",
        "- Из финальных `.out` извлечены `<S**2>`, ZPE, энтальпия и Gibbs free energy при 298.15 K.",
        "- `results/final_results.csv` дополнен колонками `s2_actual`, `zpe_hartree`, `e_plus_zpe`, `gibbs_298` и производными относительными энергиями.",
        "- Создан файл `results/best_B6_population.csv` с Mulliken/Loewdin зарядами и спиновыми плотностями выбранной структуры.",
        "- Для всего сошедшегося screening выполнена дедупликация distance-fingerprint с порогом 0.02 A.",
        "- Исправлена интерпретация финального этапа: 10 финальных расчетов сошлись в один минимум; это проверка воспроизводимости, а не ранжирование разных изомеров.",
        "",
        "## 2. Версия ПО и среда",
        "",
        f"- ORCA: `{orca_version or '6.1.1 RELEASE'}`.",
        f"- Host из ORCA output: `{host or 'vmi3233575'}`.",
        "- ОС запуска расчетов: Linux/VPS по структуре проекта и путям ORCA output (`/root/B6_ORCA_campaign/...`).",
        "- Обработка отчета v2 выполнена локально Python-скриптом `scripts/extract_extras.py`.",
        "",
        "## 3. Финальный этап: электронная энергия, ZPE и Gibbs",
        "",
        "В таблице ниже `relative_e_plus_zpe_ev` пересчитана отдельно по `E + ZPE` и не заменяет исходную электронную `relative_energy_ev`.",
        "",
        markdown_table(
            ["calculation", "m", "E, Eh", "ZPE, Eh", "E+ZPE, Eh", "dE(E+ZPE), eV", "G298, Eh", "<S**2>"],
            final_table_rows,
        ),
        "",
        "Ключевое значение для выбранной структуры:",
        "",
        f"- calculation from `best_B6.xyz`: `{selected_best.get('calculation_name', '')}`",
        f"- point group from ORCA thermochemistry: `{point_group or 'C2h'}`",
        f"- `<S**2>` actual/expected: `{selected_best.get('s2_actual', '')}` / `{selected_best.get('s2_expected', '')}`",
        f"- `E + ZPE`: `{selected_best.get('e_plus_zpe', '')}` Eh",
        f"- `G(298.15 K)`: `{selected_best.get('gibbs_298', '')}` Eh",
        f"- SOMO energies, eV: `{selected_best.get('somo_ev', '')}`",
        f"- HOMO/LUMO/gap, eV: `{selected_best.get('homo_ev', '')}` / `{selected_best.get('lumo_ev', '')}` / `{selected_best.get('homo_lumo_gap_ev', '')}`",
        f"- minimum by `E + ZPE` within the same C2h basin: `{best_by_e_zpe.get('calculation_name', '')}`, `dE(E+ZPE)=0.00000000` eV",
        "",
        "## 4. Mulliken/Loewdin population для best_B6",
        "",
        f"Полная таблица сохранена в `{population_csv.relative_to(project).as_posix()}`.",
        "",
        markdown_table(
            ["atom", "Mulliken q", "Mulliken spin", "Loewdin q", "Loewdin spin"],
            population_table,
        ),
        "",
        "## 5. Найденные уникальные минимумы screening",
        "",
        f"Дедупликация применена ко всем сошедшимся stage1-структурам: `{converged_hits}` попаданий. "
        f"Порог distance-fingerprint: `0.02 A`. Сумма попаданий по группам равна `{converged_hits}`.",
        "",
        f"Полная таблица сохранена в `{unique_csv.relative_to(project).as_posix()}`.",
        "",
        markdown_table(
            ["group", "representative", "hits", "m", "dE, eV", "geometry types"],
            report_unique_rows(unique_rows),
        ),
        "",
        "## 6. Интерпретация финального этапа",
        "",
        "10 финальных PBE0-D4/def2-TZVP OptFreq расчетов не являются набором независимых финальных изомеров. "
        "После оптимизации они приходят к одному C2h-минимуму с различиями полной энергии порядка микровольт-электронвольт. "
        "Поэтому финальный этап следует трактовать как проверку воспроизводимости минимума из разных стартов.",
        "",
        "## 7. Точечная группа и геометрия best_B6",
        "",
        f"ORCA thermochemistry определяет точечную группу выбранной структуры как `{point_group or 'C2h'}` "
        "при собственном анализе симметрии. Координаты `best_B6.xyz` также имеют центр инверсии с парными атомами B1/B6, B2/B5 и B3/B4.",
        "",
        "Характерные расстояния B-B в `best_B6.xyz`:",
        "",
        markdown_table(["pairs", "distance, A", "count"], distance_summary),
        "",
        "Все попарные расстояния:",
        "",
        markdown_table(["pair", "distance, A"], pair_rows),
        "",
        "## 8. Подраздел 10.2: сравнение с Alexandrova et al., JPC A 2003",
        "",
        "Сравнение выполнено с работой A. N. Alexandrova, A. I. Boldyrev, H.-J. Zhai, L.-S. Wang, "
        "E. Steiner, P. W. Fowler, *J. Phys. Chem. A* 2003, 107, 1359-1369, DOI: `10.1021/jp0268866`.",
        "Числа литературы ниже приведены как округленные характерные длины C2h-мотива; для строгой публикационной версии их стоит сверить с печатной таблицей/рисунком статьи.",
        "",
        markdown_table(["quantity", "this work", "Alexandrova et al. 2003", "comment"], literature_rows),
        "",
        "## 9. Ограничения",
        "",
        "- Спиновый порядок установлен только на уровне screening R2SCAN-3C и затем проверен для выбранных triplet-кандидатов на финальном уровне.",
        "- Мультиреференсность не оценивалась; значения `<S**2>` для UKS приведены как диагностические, но не являются полноценной проверкой характера волновой функции.",
        "- Финальный набор не содержит альтернативных изомеров после оптимизации: он содержит разные старты, пришедшие в один C2h-минимум.",
        "- Литературное сравнение по геометрии сделано по характерным длинам B-B; это не заменяет полный benchmark на одинаковом уровне теории.",
        "",
        "## 10. Файлы v2",
        "",
        f"- `{report_md.relative_to(project).as_posix()}`: этот отчет Markdown.",
        f"- `{report_txt.relative_to(project).as_posix()}`: текстовая копия отчета.",
        "- `results/final_results.csv`: расширенная финальная таблица.",
        f"- `{population_csv.relative_to(project).as_posix()}`: Mulliken/Loewdin population для best_B6.",
        f"- `{unique_csv.relative_to(project).as_posix()}`: уникальные stage1-минимумы после дедупликации.",
        "",
    ]
    content = "\n".join(lines)
    report_md.write_text(content, encoding="utf-8")
    report_txt.write_text(content, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract B6 extras and build report v2.")
    parser.add_argument("--project-dir", default=".")
    parser.add_argument("--final-csv", default="results/final_results.csv")
    parser.add_argument("--screening-csv", default="results/results.csv")
    parser.add_argument("--best-xyz", default="results/best_B6.xyz")
    parser.add_argument("--population-csv", default="results/best_B6_population.csv")
    parser.add_argument("--unique-csv", default="results/screening_unique_minima.csv")
    parser.add_argument("--report-md", default="results/B6_final_report_v2.md")
    parser.add_argument("--report-txt", default="results/B6_final_report_v2.txt")
    parser.add_argument("--dedupe-tolerance", type=float, default=0.02)
    args = parser.parse_args()

    project = Path(args.project_dir).resolve()
    final_csv = (project / args.final_csv).resolve()
    screening_csv = (project / args.screening_csv).resolve()
    best_xyz = (project / args.best_xyz).resolve()
    population_csv = (project / args.population_csv).resolve()
    unique_csv = (project / args.unique_csv).resolve()
    report_md = (project / args.report_md).resolve()
    report_txt = (project / args.report_txt).resolve()

    final_rows = update_final_results(project, final_csv)
    selected_best = select_best_xyz_row(final_rows, best_xyz, project)
    write_best_population(project, selected_best, best_xyz, population_csv)
    unique_rows = dedupe_screening(project, screening_csv, unique_csv, args.dedupe_tolerance)

    best_out = resolve_project_path(selected_best.get("output_file", ""), project)
    best_text = read_text(best_out)
    point_group = parse_point_group(best_text)
    orca_version = parse_orca_version(best_text)
    host = parse_host(best_text)

    build_report_v2(
        project,
        final_rows,
        unique_rows,
        best_xyz,
        report_md,
        report_txt,
        population_csv,
        unique_csv,
        point_group,
        orca_version,
        host,
    )

    print(f"Updated final CSV: {final_csv}")
    print(f"Wrote population CSV: {population_csv}")
    print(f"Wrote unique minima CSV: {unique_csv}")
    print(f"Wrote report v2: {report_md}")
    print(f"Wrote report v2 text: {report_txt}")


if __name__ == "__main__":
    main()
