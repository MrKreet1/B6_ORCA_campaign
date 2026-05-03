#!/usr/bin/env python3
"""
collect_results.py

Обходит папки ORCA-расчётов, читает реальные .out-файлы и собирает results.csv.
Ничего не придумывает: энергии, частоты и координаты берутся только из файлов ORCA.

Критерий true minimum:
- ORCA TERMINATED NORMALLY;
- OPTIMIZATION HAS CONVERGED;
- есть частотный блок для финального расчёта;
- нет отрицательных vibrational frequencies.

Для stage1 без Freq has_imaginary_frequencies остаётся пустым/NA, а is_true_minimum=False.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

HARTREE_TO_EV = 27.211386245988
FREQUENCY_ZERO_THRESHOLD_CM = 10.0

COLUMNS = [
    "calculation_name",
    "geometry_type",
    "distance",
    "charge",
    "multiplicity",
    "method",
    "basis",
    "total_energy_hartree",
    "relative_energy_ev",
    "normal_termination",
    "optimization_converged",
    "has_imaginary_frequencies",
    "n_imaginary_frequencies",
    "lowest_frequency_cm-1",
    "is_true_minimum",
    "xyz_file",
    "output_file",
]


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def find_last_energy(text: str) -> Optional[float]:
    # ORCA обычно пишет: FINAL SINGLE POINT ENERGY     -XXX.XXXXXXXX
    values = re.findall(r"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+(-?\d+\.\d+(?:[Ee][+-]?\d+)?)", text)
    if not values:
        return None
    return float(values[-1])


def normal_termination(text: str) -> bool:
    return "ORCA TERMINATED NORMALLY" in text


def opt_converged(text: str) -> bool:
    patterns = [
        "THE OPTIMIZATION HAS CONVERGED",
        "OPTIMIZATION HAS CONVERGED",
        "Geometry convergence",
    ]
    if any(p.lower() in text.lower() for p in patterns[:2]):
        return True
    # Не считаем простой блок Geometry convergence доказательством сходимости.
    return False


def parse_metadata(calc_dir: Path) -> Dict[str, object]:
    meta = calc_dir / "metadata.json"
    if meta.exists():
        try:
            return json.loads(meta.read_text(encoding="utf-8"))
        except Exception:
            pass
    return {}


def infer_from_name(name: str) -> Dict[str, object]:
    info: Dict[str, object] = {}
    m = re.search(r"B6_(.+)_d([0-9.]+)_q(-?\d+)_m(\d+)_", name)
    if m:
        info["geometry_type"] = m.group(1)
        try:
            info["distance"] = float(m.group(2))
        except ValueError:
            info["distance"] = ""
        info["charge"] = int(m.group(3))
        info["multiplicity"] = int(m.group(4))
    return info


def parse_frequencies(text: str) -> List[float]:
    """Извлекает частоты из блока VIBRATIONAL FREQUENCIES ORCA.

    Типичные строки:
      0:       0.00 cm**-1
      6:     123.45 cm**-1
      7:     -50.12 cm**-1 ***imaginary mode***
    """
    if "VIBRATIONAL FREQUENCIES" not in text:
        return []

    start = text.rfind("VIBRATIONAL FREQUENCIES")
    block = text[start:]
    # Останавливаемся до следующего крупного раздела, если он есть.
    stop_markers = [
        "NORMAL MODES",
        "IR SPECTRUM",
        "THERMOCHEMISTRY",
        "ORCA TERMINATED",
    ]
    stops = [block.find(marker) for marker in stop_markers if block.find(marker) > 0]
    if stops:
        block = block[: min(stops)]

    freqs: List[float] = []
    for line in block.splitlines():
        # Требуем cm**-1, чтобы не ловить другие числа.
        if "cm**-1" not in line and "cm-1" not in line:
            continue
        m = re.search(r"^\s*\d+\s*:\s*(-?\d+(?:\.\d+)?)\s*cm", line)
        if m:
            freqs.append(float(m.group(1)))
    return freqs


def frequency_status(freqs: List[float]) -> Tuple[str, str, str]:
    """Возвращает has_imaginary, count, lowest_frequency.

    Пустая строка означает, что Freq не был найден, например для stage1 Opt.
    Нулевые/почти нулевые трансляционно-вращательные моды не учитываются
    при выборе lowest_frequency_cm-1 и при подсчете мнимых вибрационных мод.
    """
    if not freqs:
        return "", "", ""
    vibrational = [f for f in freqs if abs(f) > FREQUENCY_ZERO_THRESHOLD_CM]
    if not vibrational:
        return "", "", ""
    negative = [f for f in vibrational if f < -FREQUENCY_ZERO_THRESHOLD_CM]
    lowest = min(vibrational)
    has_imag = bool(negative)
    return str(has_imag), str(len(negative)), f"{lowest:.6f}"


def parse_last_cartesian_block(text: str) -> Optional[List[Tuple[str, float, float, float]]]:
    """Берёт последний блок CARTESIAN COORDINATES (ANGSTROEM) из .out."""
    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines) if "CARTESIAN COORDINATES (ANGSTROEM)" in line]
    if not starts:
        return None

    for start in reversed(starts):
        atoms: List[Tuple[str, float, float, float]] = []
        i = start + 1
        # Пропускаем пустые строки и разделители.
        while i < len(lines) and (not lines[i].strip() or set(lines[i].strip()) <= {"-"}):
            i += 1
        while i < len(lines):
            line = lines[i].strip()
            if not line or set(line) <= {"-"}:
                break
            parts = line.split()
            if len(parts) >= 4 and re.fullmatch(r"[A-Za-z]{1,3}", parts[0]):
                try:
                    atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
                    i += 1
                    continue
                except ValueError:
                    pass
            break
        if len(atoms) == 6:
            return atoms
    return None


def write_xyz(path: Path, atoms: List[Tuple[str, float, float, float]], comment: str) -> None:
    lines = [f"{len(atoms)}", comment]
    for el, x, y, z in atoms:
        lines.append(f"{el:2s} {x:16.8f} {y:16.8f} {z:16.8f}")
    content = "\n".join(lines) + "\n"
    if path.exists() and path.read_text(encoding="utf-8", errors="replace") == content:
        return
    path.write_text(content, encoding="utf-8")


def write_rows_csv(path: Path, rows: List[Dict[str, object]], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in columns})


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def portable_path(path_value: object, project_dir: Path) -> str:
    """Return a repository-relative POSIX path when possible."""
    if not path_value:
        return ""
    raw = Path(str(path_value))
    try:
        resolved = raw.resolve() if raw.is_absolute() else (project_dir / raw).resolve()
        return resolved.relative_to(project_dir).as_posix()
    except Exception:
        pass

    marker = project_dir.name
    parts = raw.parts
    if marker in parts:
        idx = parts.index(marker)
        tail = parts[idx + 1 :]
        if tail:
            return Path(*tail).as_posix()
    return str(path_value).replace("\\", "/")


def resolve_project_path(path_text: str, project_dir: Path) -> Path:
    raw = Path(path_text)
    candidates = [raw]
    if not raw.is_absolute():
        candidates.insert(0, project_dir / raw)

    marker = project_dir.name
    parts = raw.parts
    if marker in parts:
        idx = parts.index(marker)
        tail = parts[idx + 1 :]
        if tail:
            candidates.append(project_dir / Path(*tail))

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[0]


def source_family(geometry_type: object) -> str:
    text = str(geometry_type).lower()
    if any(token in text for token in ["random", "3d", "octa", "prism", "pyramid"]):
        return "3D/random start"
    if any(token in text for token in ["planar", "ring", "triangle", "rhombic", "rectangular", "quasi"]):
        return "planar/quasi-planar start"
    return "other start"


def write_final_report(report_path: Path, rows: List[Dict[str, object]], output_csv: Path, best_xyz: Path) -> None:
    screening_csv = report_path.parent / "results.csv"
    screening_rows = read_csv_rows(screening_csv)
    screening_completed = [r for r in screening_rows if r.get("normal_termination") == "True"]
    screening_converged = [r for r in screening_rows if r.get("optimization_converged") == "True"]
    completed = [r for r in rows if r.get("normal_termination") == "True"]
    converged = [r for r in rows if r.get("optimization_converged") == "True"]
    true_minima = [r for r in rows if r.get("is_true_minimum") == "True"]
    best = true_minima[0] if true_minima else None
    final_families: Dict[str, int] = {}
    for row in rows:
        family = source_family(row.get("geometry_type", ""))
        final_families[family] = final_families.get(family, 0) + 1

    top_rows = []
    for row in rows[:10]:
        top_rows.append(
            " | ".join(
                [
                    str(row.get("calculation_name", "")),
                    f"m={row.get('multiplicity', '')}",
                    f"E={row.get('total_energy_hartree', '')}",
                    f"dE={row.get('relative_energy_ev', '')} eV",
                    f"imag={row.get('n_imaginary_frequencies', '')}",
                    f"true_min={row.get('is_true_minimum', '')}",
                ]
            )
        )

    lines = [
        "Финальный отчет по расчетному исследованию B6",
        "",
        "1. Введение",
        "Малые кластеры бора чувствительны к стартовой геометрии и спиновому состоянию, поэтому поиск минимума B6 требует многостартового подхода. В расчете проверялись плоские, квазиплоские и 3D-старты, а финальный выбор выполнялся только после частотного анализа.",
        "",
        "2. Цель и задачи",
        "Цель: найти наиболее устойчивую геометрию нейтрального B6 в рамках заданного набора DFT-расчетов. Задачи: сгенерировать старты, провести R2SCAN-3C Opt screening, выбрать низкоэнергетические кандидаты, выполнить PBE0-D4/def2-TZVP Opt Freq и исключить структуры с мнимыми частотами.",
        "",
        "3. Методика",
        "- ПО: ORCA 6.1",
        "- Заряд: 0",
        "- Мультиплетности: 1, 3, 5",
        "- Screening: R2SCAN-3C Opt",
        "- Финальный уровень: PBE0-D4/def2-TZVP Opt Freq",
        "- Ресурсы в input-файлах: nprocs 8, maxcore 2500 MB",
        "- Критерий минимума: normal_termination=True, optimization_converged=True, n_imaginary_frequencies=0",
        "",
        "4. Генерация стартовых геометрий",
        "Генератор покрывает линейные, кольцевые, искаженные кольцевые, компактные плоские, ромбические, прямоугольные, fused-triangle, квазиплоские, октаэдрические, призматические, пирамидальные и random 3D-старты. Основной диапазон расстояний B-B: 1.5, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 3.5 Angstrom.",
        "",
        "5. Screening R2SCAN-3C Opt",
        f"- Таблица screening: {screening_csv if screening_rows else 'не найдена рядом с отчетом'}",
        f"- Строк в screening-таблице: {len(screening_rows)}",
        f"- Нормально завершено: {len(screening_completed)}",
        f"- Сошедшаяся оптимизация: {len(screening_converged)}",
        "Энергии screening извлекались из FINAL SINGLE POINT ENERGY только в реальных ORCA output-файлах.",
        "",
        "6. Отбор финальных кандидатов",
        "Финальные input-файлы готовятся из низкоэнергетических сошедшихся структур. Дедупликация выполняется по отсортированному набору межатомных расстояний B-B с настраиваемым порогом.",
        "",
        "7. Финальные PBE0-D4/def2-TZVP Opt Freq расчеты",
        f"- Таблица final: {output_csv}",
        f"- Всего финальных output-файлов в таблице: {len(rows)}",
        f"- ORCA TERMINATED NORMALLY: {len(completed)}",
        f"- THE OPTIMIZATION HAS CONVERGED: {len(converged)}",
        f"- Истинных минимумов без мнимых частот: {len(true_minima)}",
        f"- Источники финальных кандидатов: {', '.join(f'{k}: {v}' for k, v in sorted(final_families.items()))}",
        "",
        "8. Частотный анализ",
        "Структуры с n_imaginary_frequencies > 0 исключаются из финального выбора, даже если их электронная энергия ниже. В текущей финальной таблице критерий истинного минимума берется из колонок has_imaginary_frequencies, n_imaginary_frequencies и is_true_minimum.",
        "",
        "9. Энергии финальных кандидатов",
        "calculation_name | multiplicity | total_energy_hartree | relative_energy_ev | n_imaginary_frequencies | is_true_minimum",
        *top_rows,
        "",
    ]

    if best:
        lines.extend(
            [
                "10. Выбор финального минимума",
                f"- calculation_name: {best.get('calculation_name', '')}",
                f"- geometry_type: {best.get('geometry_type', '')}",
                f"- multiplicity: {best.get('multiplicity', '')}",
                f"- method/basis: {best.get('method', '')}/{best.get('basis', '')}",
                f"- total_energy_hartree: {best.get('total_energy_hartree', '')}",
                f"- relative_energy_ev: {best.get('relative_energy_ev', '')}",
                f"- lowest_frequency_cm-1: {best.get('lowest_frequency_cm-1', '')}",
                f"- n_imaginary_frequencies: {best.get('n_imaginary_frequencies', '')}",
                f"- xyz_file: {best.get('xyz_file', '')}",
                f"- output_file: {best.get('output_file', '')}",
                f"- best_B6.xyz: {best_xyz}",
                "",
                "Финальная структура выбрана как структура с минимальной полной энергией среди расчетов, которые завершились нормально, имели сошедшуюся оптимизацию и не содержали мнимых частот.",
            ]
        )
    else:
        lines.extend(
            [
                "Финальная структура не выбрана: среди обработанных расчетов нет истинного минимума по заданным критериям.",
                "Нужно проверить несошедшиеся расчеты, мнимые частоты и при необходимости перезапустить Opt Freq.",
            ]
        )

    lines.extend(
        [
            "",
            "11. Обсуждение плоских и 3D-стартов",
            "В финальном наборе присутствуют кандидаты, происходящие как из planar/quasi-planar, так и из 3D/random стартов. Если 3D-старты после оптимизации дают ту же низкоэнергетическую область, это поддерживает вывод о необходимости многостартовой проверки, а не выбора одной заранее заданной геометрии.",
            "",
            "12. Ограничения работы",
            "Вывод справедлив для реально выполненного набора стартов, мультиплетностей и уровней теории. Новые старты, более плотная дедупликация, другие функционалы или учет дополнительных поправок могут изменить относительный порядок близких изомеров.",
            "",
            "13. Приложения",
            f"- final_results.csv: {output_csv}",
            f"- best_B6.xyz: {best_xyz}",
            f"- results.csv: {screening_csv if screening_rows else 'не найден'}",
            "",
            "Все численные энергии, частоты и координаты в этом отчете взяты из ORCA output-файлов и сгенерированных CSV/XYZ файлов. Фиктивные значения не использовались.",
        ]
    )

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote final report: {report_path}")


def collect(
    root: Path,
    output_csv: Path,
    best_xyz: Path,
    all_energies_csv: Optional[Path] = None,
    report_path: Optional[Path] = None,
    project_dir: Optional[Path] = None,
) -> List[Dict[str, object]]:
    project_dir = (project_dir or Path.cwd()).resolve()
    out_files = sorted(root.rglob("*.out"))
    rows: List[Dict[str, object]] = []

    for out_path in out_files:
        calc_dir = out_path.parent
        text = read_text(out_path)
        meta = parse_metadata(calc_dir)
        inferred = infer_from_name(calc_dir.name)

        energy = find_last_energy(text)
        freqs = parse_frequencies(text)
        has_imag, n_imag, lowest_freq = frequency_status(freqs)
        normal = normal_termination(text)
        conv = opt_converged(text)

        calc_name = str(meta.get("calculation_name") or calc_dir.name)
        opt_xyz_path = calc_dir / f"{calc_name}_optimized.xyz"
        atoms = parse_last_cartesian_block(text)
        if atoms:
            write_xyz(opt_xyz_path, atoms, f"Optimized geometry from {out_path.name}")
            xyz_file = portable_path(opt_xyz_path, project_dir)
        else:
            xyz_file = portable_path(meta.get("xyz_file") or "", project_dir)

        # true minimum требует наличия частотного расчёта. Для stage1 Opt без частот это False.
        has_freq = bool(freqs)
        is_true_min = bool(normal and conv and has_freq and has_imag == "False" and n_imag == "0")

        row: Dict[str, object] = {
            "calculation_name": calc_name,
            "geometry_type": meta.get("geometry_type", inferred.get("geometry_type", "")),
            "distance": meta.get("distance", inferred.get("distance", "")),
            "charge": meta.get("charge", inferred.get("charge", "")),
            "multiplicity": meta.get("multiplicity", inferred.get("multiplicity", "")),
            "method": meta.get("method", ""),
            "basis": meta.get("basis", ""),
            "total_energy_hartree": f"{energy:.12f}" if energy is not None else "",
            "relative_energy_ev": "",
            "normal_termination": str(normal),
            "optimization_converged": str(conv),
            "has_imaginary_frequencies": has_imag,
            "n_imaginary_frequencies": n_imag,
            "lowest_frequency_cm-1": lowest_freq,
            "is_true_minimum": str(is_true_min),
            "xyz_file": xyz_file,
            "output_file": portable_path(out_path, project_dir),
        }
        rows.append(row)

    # Сортировка: сначала строки с энергией, затем без энергии.
    def energy_key(row: Dict[str, object]) -> float:
        try:
            return float(row["total_energy_hartree"])
        except Exception:
            return math.inf

    rows.sort(key=energy_key)

    energies = [float(r["total_energy_hartree"]) for r in rows if str(r["total_energy_hartree"]).strip()]
    emin = min(energies) if energies else None
    if emin is not None:
        for r in rows:
            if str(r["total_energy_hartree"]).strip():
                rel = (float(r["total_energy_hartree"]) - emin) * HARTREE_TO_EV
                r["relative_energy_ev"] = f"{rel:.8f}"

    write_rows_csv(output_csv, rows, COLUMNS)
    if all_energies_csv is not None:
        write_rows_csv(all_energies_csv, rows, COLUMNS)

    # best_B6.xyz выбирается только среди true minima.
    best_xyz.parent.mkdir(parents=True, exist_ok=True)
    true_minima = [r for r in rows if r.get("is_true_minimum") == "True" and str(r.get("xyz_file", "")).strip()]
    if true_minima:
        src = resolve_project_path(str(true_minima[0]["xyz_file"]), project_dir)
        if src.exists():
            shutil.copyfile(src, best_xyz)
            print(f"Best true minimum saved to: {best_xyz}")
        else:
            print(f"WARNING: best xyz source does not exist: {src}")
    else:
        print("No true minimum found. best_B6.xyz was not created/updated.")

    print(f"Wrote CSV: {output_csv}")
    if all_energies_csv is not None:
        print(f"Wrote all energies CSV: {all_energies_csv}")
    if report_path is not None:
        write_final_report(report_path, rows, output_csv, best_xyz)
    print(f"Parsed output files: {len(out_files)}")
    return rows


def main() -> None:
    p = argparse.ArgumentParser(description="Collect ORCA results for B6 campaign.")
    p.add_argument("--root", default="calculations", help="Корень для поиска .out файлов.")
    p.add_argument("--csv", default="results/results.csv", help="CSV для записи результатов.")
    p.add_argument("--best-xyz", default="results/best_B6.xyz", help="Файл для лучшей структуры без мнимых частот.")
    p.add_argument("--all-energies-csv", default="", help="Дополнительная CSV-таблица всех энергий, если нужна.")
    p.add_argument("--report", default="", help="Итоговый текстовый отчет, если нужен.")
    p.add_argument("--project-dir", default=".", help="Корень проекта для записи переносимых относительных путей в CSV.")
    args = p.parse_args()

    project_dir = Path(args.project_dir).resolve()

    root = Path(args.root)
    if not root.is_absolute():
        root = project_dir / root
    output_csv = Path(args.csv)
    if not output_csv.is_absolute():
        output_csv = project_dir / output_csv
    best_xyz = Path(args.best_xyz)
    if not best_xyz.is_absolute():
        best_xyz = project_dir / best_xyz
    all_energies_csv = Path(args.all_energies_csv) if args.all_energies_csv else None
    if all_energies_csv is not None and not all_energies_csv.is_absolute():
        all_energies_csv = project_dir / all_energies_csv
    report_path = Path(args.report) if args.report else None
    if report_path is not None and not report_path.is_absolute():
        report_path = project_dir / report_path

    collect(
        root.resolve(),
        output_csv.resolve(),
        best_xyz.resolve(),
        all_energies_csv.resolve() if all_energies_csv is not None else None,
        report_path.resolve() if report_path is not None else None,
        project_dir,
    )


if __name__ == "__main__":
    main()
