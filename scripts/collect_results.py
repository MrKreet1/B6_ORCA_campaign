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
    """
    if not freqs:
        return "", "", ""
    negative = [f for f in freqs if f < -1.0]  # небольшой порог против численного шума около нуля
    lowest = min(freqs)
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
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(atoms)}\n")
        f.write(comment + "\n")
        for el, x, y, z in atoms:
            f.write(f"{el:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")


def collect(root: Path, output_csv: Path, best_xyz: Path) -> List[Dict[str, object]]:
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
            xyz_file = str(opt_xyz_path)
        else:
            xyz_file = str(meta.get("xyz_file") or "")

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
            "output_file": str(out_path),
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

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in COLUMNS})

    # best_B6.xyz выбирается только среди true minima.
    best_xyz.parent.mkdir(parents=True, exist_ok=True)
    true_minima = [r for r in rows if r.get("is_true_minimum") == "True" and str(r.get("xyz_file", "")).strip()]
    if true_minima:
        src = Path(str(true_minima[0]["xyz_file"]))
        if src.exists():
            shutil.copyfile(src, best_xyz)
            print(f"Best true minimum saved to: {best_xyz}")
        else:
            print(f"WARNING: best xyz source does not exist: {src}")
    else:
        print("No true minimum found. best_B6.xyz was not created/updated.")

    print(f"Wrote CSV: {output_csv}")
    print(f"Parsed output files: {len(out_files)}")
    return rows


def main() -> None:
    p = argparse.ArgumentParser(description="Collect ORCA results for B6 campaign.")
    p.add_argument("--root", default="calculations", help="Корень для поиска .out файлов.")
    p.add_argument("--csv", default="results/results.csv", help="CSV для записи результатов.")
    p.add_argument("--best-xyz", default="results/best_B6.xyz", help="Файл для лучшей структуры без мнимых частот.")
    args = p.parse_args()

    collect(Path(args.root).resolve(), Path(args.csv).resolve(), Path(args.best_xyz).resolve())


if __name__ == "__main__":
    main()
