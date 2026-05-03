#!/usr/bin/env python3
"""
generate_b6_inputs.py

Создаёт воспроизводимый набор входных файлов ORCA для поиска устойчивой
геометрии нейтрального кластера B6 при разных стартовых геометриях,
межатомных расстояниях и мультиплетностях.

Скрипт НЕ генерирует энергии и частоты. Энергии/частоты должны быть получены
только из реальных output-файлов ORCA и затем собраны collect_results.py.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import re
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

Atom = Tuple[str, float, float, float]
Coord = Tuple[float, float, float]


def slug(text: str) -> str:
    """Безопасное имя для папок/файлов."""
    text = text.strip().replace("/", "-")
    text = re.sub(r"[^A-Za-z0-9_.+-]+", "_", text)
    return text.strip("_")


def center(coords: Sequence[Coord]) -> List[Coord]:
    cx = sum(x for x, _, _ in coords) / len(coords)
    cy = sum(y for _, y, _ in coords) / len(coords)
    cz = sum(z for _, _, z in coords) / len(coords)
    return [(x - cx, y - cy, z - cz) for x, y, z in coords]


def atoms_from_coords(coords: Sequence[Coord]) -> List[Atom]:
    return [("B", x, y, z) for x, y, z in center(coords)]


def linear_chain(d: float) -> List[Atom]:
    return atoms_from_coords([((i - 2.5) * d, 0.0, 0.0) for i in range(6)])


def planar_ring(d: float) -> List[Atom]:
    # Для правильного шестиугольника сторона = радиус описанной окружности = d.
    coords: List[Coord] = []
    for i in range(6):
        a = 2.0 * math.pi * i / 6.0
        coords.append((d * math.cos(a), d * math.sin(a), 0.0))
    return atoms_from_coords(coords)


def compact_planar_triangle(d: float) -> List[Atom]:
    # Компактный фрагмент треугольной решётки 3+2+1.
    h = math.sqrt(3.0) * d / 2.0
    coords = [
        (-d, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (d, 0.0, 0.0),
        (-0.5 * d, h, 0.0),
        (0.5 * d, h, 0.0),
        (0.0, 2.0 * h, 0.0),
    ]
    return atoms_from_coords(coords)


def rhombic_planar(d: float) -> List[Atom]:
    # Параллелограмм 2x3 с углом 60°.
    h = math.sqrt(3.0) * d / 2.0
    coords: List[Coord] = []
    for row in range(2):
        for col in range(3):
            coords.append((col * d + 0.5 * row * d, row * h, 0.0))
    return atoms_from_coords(coords)


def rectangular_planar(d: float) -> List[Atom]:
    coords: List[Coord] = []
    for row in range(2):
        for col in range(3):
            coords.append(((col - 1.0) * d, (row - 0.5) * d, 0.0))
    return atoms_from_coords(coords)


def octahedral_3d(d: float) -> List[Atom]:
    # У октаэдра ребро = sqrt(2)*a, значит a = d/sqrt(2).
    a = d / math.sqrt(2.0)
    coords = [(a, 0.0, 0.0), (-a, 0.0, 0.0), (0.0, a, 0.0), (0.0, -a, 0.0), (0.0, 0.0, a), (0.0, 0.0, -a)]
    return atoms_from_coords(coords)


def trigonal_prism(d: float) -> List[Atom]:
    # Две равносторонние тройки, разделённые на d по оси z.
    r = d / math.sqrt(3.0)
    coords: List[Coord] = []
    for z in (-0.5 * d, 0.5 * d):
        for i in range(3):
            a = 2.0 * math.pi * i / 3.0 + math.pi / 6.0
            coords.append((r * math.cos(a), r * math.sin(a), z))
    return atoms_from_coords(coords)


def min_distance(coords: Sequence[Coord]) -> float:
    md = float("inf")
    for i in range(len(coords)):
        xi, yi, zi = coords[i]
        for j in range(i + 1, len(coords)):
            xj, yj, zj = coords[j]
            r = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)
            md = min(md, r)
    return md


def rotate_xyz(coords: Sequence[Coord], ax: float, ay: float, az: float) -> List[Coord]:
    """Простая 3D-ротация вокруг x/y/z."""
    sx, cx = math.sin(ax), math.cos(ax)
    sy, cy = math.sin(ay), math.cos(ay)
    sz, cz = math.sin(az), math.cos(az)
    out: List[Coord] = []
    for x, y, z in coords:
        y, z = y * cx - z * sx, y * sx + z * cx
        x, z = x * cy + z * sy, -x * sy + z * cy
        x, y = x * cz - y * sz, x * sz + y * cz
        out.append((x, y, z))
    return out


def random_3d(d: float, seed: int) -> List[Atom]:
    """Быстрая физически разумная случайная 3D-геометрия."""
    rng = random.Random(seed)
    if seed % 2 == 0:
        base = [(x, y, z) for _, x, y, z in trigonal_prism(d)]
    else:
        base = [(x, y, z) for _, x, y, z in octahedral_3d(d)]
    rotated = rotate_xyz(
        base,
        rng.uniform(0.0, 2.0 * math.pi),
        rng.uniform(0.0, 2.0 * math.pi),
        rng.uniform(0.0, 2.0 * math.pi),
    )
    jitter = min(0.18, 0.08 * d)
    coords = [
        (
            x + rng.uniform(-jitter, jitter),
            y + rng.uniform(-jitter, jitter),
            z + rng.uniform(-jitter, jitter),
        )
        for x, y, z in rotated
    ]
    return atoms_from_coords(coords)


def all_geometries(d: float, n_random: int, random_seed: int) -> List[Tuple[str, List[Atom]]]:
    geoms: List[Tuple[str, List[Atom]]] = [
        ("linear_chain", linear_chain(d)),
        ("planar_ring", planar_ring(d)),
        ("compact_planar_triangle", compact_planar_triangle(d)),
        ("rhombic_planar", rhombic_planar(d)),
        ("rectangular_planar", rectangular_planar(d)),
        ("octahedral_3d", octahedral_3d(d)),
        ("trigonal_prism", trigonal_prism(d)),
    ]
    for k in range(n_random):
        seed = random_seed + k
        geoms.append((f"random_3d_seed{seed}", random_3d(d, seed)))
    return geoms


def write_xyz(path: Path, atoms: Sequence[Atom], comment: str) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(atoms)}\n")
        f.write(comment + "\n")
        for el, x, y, z in atoms:
            f.write(f"{el:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")


def atoms_to_orca_block(atoms: Sequence[Atom], charge: int, multiplicity: int) -> str:
    lines = [f"* xyz {charge} {multiplicity}"]
    for el, x, y, z in atoms:
        lines.append(f"{el:2s} {x:16.8f} {y:16.8f} {z:16.8f}")
    lines.append("*")
    return "\n".join(lines)


def make_input(
    atoms: Sequence[Atom],
    charge: int,
    multiplicity: int,
    method: str,
    basis: str,
    task: str,
    extra_keywords: str,
    nprocs: int,
    maxcore: int,
) -> str:
    keywords = [method]
    if basis.strip():
        keywords.append(basis.strip())
    if extra_keywords.strip():
        keywords.extend(extra_keywords.strip().split())
    if task.strip():
        keywords.append(task.strip())

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

{atoms_to_orca_block(atoms, charge, multiplicity)}
"""


def parse_float_list(text: str) -> List[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def parse_int_list(text: str) -> List[int]:
    return [int(x.strip()) for x in text.split(",") if x.strip()]


def main() -> None:
    p = argparse.ArgumentParser(description="Generate B6 ORCA input campaign.")
    p.add_argument("--project-dir", default=".", help="Корень проекта. Обычно '.' если вы уже в B6_ORCA_campaign.")
    p.add_argument("--stage-dir", default="calculations/stage1", help="Папка расчётов внутри project-dir.")
    p.add_argument("--distances", default="3.5,3.0,2.5,2.2,2.0,1.8,1.6")
    p.add_argument("--multiplicities", default="1,3,5")
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--method", default="R2SCAN-3C")
    p.add_argument("--basis", default="")
    p.add_argument("--task", default="Opt")
    p.add_argument("--extra-keywords", default="TightSCF TightOpt")
    p.add_argument("--nprocs", type=int, default=8)
    p.add_argument("--maxcore", type=int, default=2500)
    p.add_argument("--n-random", type=int, default=5)
    p.add_argument("--random-seed", type=int, default=1000)
    args = p.parse_args()

    project = Path(args.project_dir).resolve()
    stage_root = project / args.stage_dir
    stage_root.mkdir(parents=True, exist_ok=True)

    distances = parse_float_list(args.distances)
    multiplicities = parse_int_list(args.multiplicities)

    jobs = []
    for d in distances:
        for geom_name, atoms in all_geometries(d, args.n_random, args.random_seed):
            for mult in multiplicities:
                calc_name = f"B6_{geom_name}_d{d:.2f}_q{args.charge}_m{mult}_{slug(args.method)}"
                calc_dir = stage_root / calc_name
                calc_dir.mkdir(parents=True, exist_ok=True)

                xyz_path = calc_dir / f"{calc_name}.xyz"
                inp_path = calc_dir / f"{calc_name}.inp"
                meta_path = calc_dir / "metadata.json"

                write_xyz(
                    xyz_path,
                    atoms,
                    f"B6 start; geometry={geom_name}; distance={d:.2f} A; charge={args.charge}; multiplicity={mult}",
                )
                inp_path.write_text(
                    make_input(atoms, args.charge, mult, args.method, args.basis, args.task, args.extra_keywords, args.nprocs, args.maxcore),
                    encoding="utf-8",
                )

                def rel(p: Path) -> str:
                    try:
                        return str(p.relative_to(project))
                    except ValueError:
                        return str(p)

                metadata = {
                    "calculation_name": calc_name,
                    "geometry_type": geom_name,
                    "distance": d,
                    "charge": args.charge,
                    "multiplicity": mult,
                    "method": args.method,
                    "basis": args.basis,
                    "task": args.task,
                    "extra_keywords": args.extra_keywords,
                    "nprocs": args.nprocs,
                    "maxcore": args.maxcore,
                    "xyz_file": rel(xyz_path),
                    "input_file": rel(inp_path),
                    "output_file": rel(inp_path.with_suffix(".out")),
                }
                meta_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
                jobs.append(str(inp_path))

    manifest = stage_root / "inputs_manifest.txt"
    manifest.write_text("\n".join(jobs) + "\n", encoding="utf-8")
    print(f"Created {len(jobs)} ORCA inputs in {stage_root}")
    print(f"Manifest: {manifest}")


if __name__ == "__main__":
    main()
