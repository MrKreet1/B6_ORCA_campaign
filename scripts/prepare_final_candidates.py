#!/usr/bin/env python3
"""
prepare_final_candidates.py

Берёт лучшие сошедшиеся структуры из results.csv этапа 1 и создаёт отдельные
папки для финального Opt Freq расчёта на более точном уровне теории.

Важно: скрипт НЕ выбирает финальную структуру. Финальный выбор делает
collect_results.py после реальных частотных расчётов.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

Atom = Tuple[str, float, float, float]


def slug(text: str) -> str:
    text = text.strip().replace("/", "-")
    text = re.sub(r"[^A-Za-z0-9_.+-]+", "_", text)
    return text.strip("_")


def read_xyz(path: Path) -> List[Atom]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 8:
        raise ValueError(f"XYZ file too short: {path}")
    try:
        n = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"Invalid XYZ atom count in {path}") from exc
    atoms: List[Atom] = []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Bad XYZ line in {path}: {line}")
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    if len(atoms) != n:
        raise ValueError(f"Expected {n} atoms, got {len(atoms)} in {path}")
    return atoms


def resolve_existing_path(path_text: str, project: Path) -> Optional[Path]:
    if not path_text:
        return None

    raw = Path(path_text)
    candidates = [raw]
    if not raw.is_absolute():
        candidates.append(project / raw)

    marker = project.name
    parts = raw.parts
    if marker in parts:
        idx = parts.index(marker)
        rel = Path(*parts[idx + 1 :])
        candidates.append(project / rel)

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return None


def distance_fingerprint(atoms: Sequence[Atom]) -> List[float]:
    distances: List[float] = []
    for i in range(len(atoms)):
        _, xi, yi, zi = atoms[i]
        for j in range(i + 1, len(atoms)):
            _, xj, yj, zj = atoms[j]
            distances.append(math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2))
    return sorted(distances)


def same_geometry(a: Sequence[float], b: Sequence[float], tolerance: float) -> bool:
    if len(a) != len(b):
        return False
    return max(abs(x - y) for x, y in zip(a, b)) <= tolerance


def write_xyz(path: Path, atoms: List[Atom], comment: str) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(atoms)}\n")
        f.write(comment + "\n")
        for el, x, y, z in atoms:
            f.write(f"{el:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")


def atoms_to_block(atoms: List[Atom], charge: int, multiplicity: int) -> str:
    lines = [f"* xyz {charge} {multiplicity}"]
    for el, x, y, z in atoms:
        lines.append(f"{el:2s} {x:16.8f} {y:16.8f} {z:16.8f}")
    lines.append("*")
    return "\n".join(lines)


def make_input(
    atoms: List[Atom],
    charge: int,
    multiplicity: int,
    method: str,
    basis: str,
    extra_keywords: str,
    nprocs: int,
    maxcore: int,
) -> str:
    keywords = [method]
    if basis.strip():
        keywords.append(basis.strip())
    if extra_keywords.strip():
        keywords.extend(extra_keywords.split())
    keywords.extend(["Opt", "Freq"])

    return f"""! {' '.join(keywords)}

%pal
  nprocs {nprocs}
end

%maxcore {maxcore}

%scf
  MaxIter 500
end

%geom
  MaxIter 300
end

{atoms_to_block(atoms, charge, multiplicity)}
"""


def row_energy(row: Dict[str, str]) -> float:
    try:
        return float(row.get("total_energy_hartree", ""))
    except Exception:
        return float("inf")


def truthy(text: str) -> bool:
    return str(text).strip().lower() == "true"


def main() -> None:
    p = argparse.ArgumentParser(description="Prepare final B6 Opt Freq candidates.")
    p.add_argument("--project-dir", default=".")
    p.add_argument("--results-csv", default="results/results.csv")
    p.add_argument("--final-dir", default="calculations/final")
    p.add_argument("--n", type=int, default=10, help="Сколько лучших кандидатов взять.")
    p.add_argument("--best-structures-csv", default="results/best_structures.csv")
    p.add_argument("--dedupe-tolerance", type=float, default=0.05, help="Порог дедупликации по отсортированным B-B расстояниям, Angstrom.")
    p.add_argument("--no-dedupe", action="store_true", help="Отключить дедупликацию финальных кандидатов.")
    p.add_argument("--method", default="PBE0")
    p.add_argument("--basis", default="def2-TZVP")
    p.add_argument("--extra-keywords", default="D4 def2/J RIJCOSX TightSCF TightOpt")
    p.add_argument("--nprocs", type=int, default=8)
    p.add_argument("--maxcore", type=int, default=2500)
    p.add_argument("--allow-unconverged", action="store_true", help="Брать даже несошедшиеся stage1 структуры, если очень нужно.")
    args = p.parse_args()

    project = Path(args.project_dir).resolve()
    csv_path = Path(args.results_csv)
    if not csv_path.is_absolute():
        csv_path = project / csv_path
    final_root = Path(args.final_dir)
    if not final_root.is_absolute():
        final_root = project / final_root
    final_root.mkdir(parents=True, exist_ok=True)

    if not csv_path.exists():
        raise FileNotFoundError(f"results.csv not found: {csv_path}")

    with csv_path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))

    candidates: List[Dict[str, str]] = []
    for row in rows:
        if not row.get("total_energy_hartree"):
            continue
        if not row.get("xyz_file"):
            continue
        xyz = resolve_existing_path(row["xyz_file"], project)
        if xyz is None:
            continue
        row["_resolved_xyz_file"] = str(xyz)
        if not args.allow_unconverged:
            if not truthy(row.get("normal_termination", "")):
                continue
            if not truthy(row.get("optimization_converged", "")):
                continue
        candidates.append(row)

    candidates.sort(key=row_energy)
    selected: List[Tuple[Dict[str, str], List[Atom]]] = []
    fingerprints: List[List[float]] = []
    for row in candidates:
        atoms = read_xyz(Path(row["_resolved_xyz_file"]))
        fp = distance_fingerprint(atoms)
        if not args.no_dedupe and any(same_geometry(fp, known, args.dedupe_tolerance) for known in fingerprints):
            continue
        selected.append((row, atoms))
        fingerprints.append(fp)
        if len(selected) >= args.n:
            break

    if not selected:
        raise RuntimeError("No suitable candidates found. Check results.csv and xyz_file paths.")

    manifest: List[str] = []
    selected_rows: List[Dict[str, str]] = []
    for rank, (row, atoms) in enumerate(selected, start=1):
        src_xyz = Path(row["_resolved_xyz_file"])
        old_name = row.get("calculation_name") or src_xyz.stem
        calc_name = f"FINAL_rank{rank:02d}_{slug(old_name)}_{slug(args.method)}_{slug(args.basis)}_OptFreq"
        calc_dir = final_root / calc_name
        calc_dir.mkdir(parents=True, exist_ok=True)

        charge = int(row.get("charge") or 0)
        multiplicity = int(row.get("multiplicity") or 1)
        geom_type = row.get("geometry_type", "")
        distance = row.get("distance", "")

        xyz_path = calc_dir / f"{calc_name}.xyz"
        inp_path = calc_dir / f"{calc_name}.inp"
        meta_path = calc_dir / "metadata.json"

        write_xyz(xyz_path, atoms, f"Final candidate rank {rank}; source={old_name}")
        inp_path.write_text(
            make_input(atoms, charge, multiplicity, args.method, args.basis, args.extra_keywords, args.nprocs, args.maxcore),
            encoding="utf-8",
        )

        metadata = {
            "calculation_name": calc_name,
            "source_calculation_name": old_name,
            "source_xyz_file": str(src_xyz),
            "geometry_type": geom_type,
            "distance": distance,
            "charge": charge,
            "multiplicity": multiplicity,
            "method": args.method,
            "basis": args.basis,
            "task": "Opt Freq",
            "extra_keywords": args.extra_keywords,
            "nprocs": args.nprocs,
            "maxcore": args.maxcore,
            "xyz_file": str(xyz_path),
            "input_file": str(inp_path),
            "output_file": str(inp_path.with_suffix(".out")),
            "stage1_energy_hartree": row.get("total_energy_hartree", ""),
            "stage1_relative_energy_ev": row.get("relative_energy_ev", ""),
        }
        meta_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
        manifest.append(str(inp_path))
        selected_row = dict(row)
        selected_row.pop("_resolved_xyz_file", None)
        selected_row.update(
            {
                "rank": str(rank),
                "final_calculation_name": calc_name,
                "final_xyz_file": str(xyz_path),
                "final_input_file": str(inp_path),
            }
        )
        selected_rows.append(selected_row)
        print(f"Prepared rank {rank}: {inp_path}")

    (final_root / "final_inputs_manifest.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    best_structures_csv = Path(args.best_structures_csv)
    if not best_structures_csv.is_absolute():
        best_structures_csv = project / best_structures_csv
    best_structures_csv.parent.mkdir(parents=True, exist_ok=True)
    base_fields = [field for field in (list(rows[0].keys()) if rows else []) if not field.startswith("_")]
    fieldnames = ["rank", "final_calculation_name", "final_xyz_file", "final_input_file"] + [f for f in base_fields if f not in {"rank", "final_calculation_name", "final_xyz_file", "final_input_file"}]
    with best_structures_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected_rows)
    print(f"Wrote selected candidate table: {best_structures_csv}")
    print(f"Prepared {len(selected)} final Opt Freq jobs in {final_root}")


if __name__ == "__main__":
    main()
