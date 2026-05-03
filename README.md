# B6_ORCA_campaign

Готовый пакет для автоматизированного поиска устойчивой геометрии кластера **B6** в **ORCA 6.1** на Linux/VPS. В архив уже включены заранее сгенерированные input-файлы этапа 1: **252 расчёта** = 7 расстояний × 12 стартовых геометрий × 3 мультиплетности.

Пакет делает:

- генерацию стартовых XYZ-геометрий B6;
- перебор геометрий, расстояний и мультиплетностей;
- создание отдельных папок расчётов;
- запуск ORCA через bash;
- повторный запуск неудачных расчётов;
- сбор энергий только из реальных `.out` файлов ORCA;
- подготовку финальных `Opt Freq` расчётов;
- проверку мнимых частот;
- выбор `best_B6.xyz` только среди структур без мнимых частот.

Численные энергии и частоты в этом пакете не заданы и не выдумываются.

---

## 1. Требования

- Linux VPS
- ORCA 6.1
- Python 3.8+
- CPU: до 8 cores
- RAM: 24 GB

Рекомендуемые настройки памяти:

```orca
%pal
  nprocs 8
end

%maxcore 2500
```

Пояснение: `2500 MB × 8 ≈ 20 GB`, остаётся запас для Linux и служебных процессов.

---

## 2. Структура проекта

```text
B6_ORCA_campaign/
├── README.md
├── setup_project.sh
├── scripts/
│   ├── generate_b6_inputs.py
│   ├── collect_results.py
│   ├── prepare_final_candidates.py
│   ├── run_all.sh
│   └── rerun_failed.sh
├── templates/
│   ├── stage1_opt_template.inp
│   └── final_opt_freq_template.inp
├── calculations/
│   ├── stage1/
│   └── final/
├── results/
└── logs/
```

---

## 3. Быстрый старт

```bash
cd ~/B6_ORCA_campaign
bash setup_project.sh

# Укажи путь к ORCA, если orca не находится через PATH:
export ORCA_CMD=/full/path/to/orca

# Проверка:
$ORCA_CMD --version
```

---

## 4. Этап 1 — генерация первичных расчётов

В архиве расчёты этапа 1 уже сгенерированы в `calculations/stage1`. Этот раздел нужен, если ты хочешь пересоздать их с другими параметрами.

По умолчанию будут созданы геометрии:

- linear_chain;
- planar_ring;
- compact_planar_triangle;
- rhombic_planar;
- rectangular_planar;
- octahedral_3d;
- trigonal_prism;
- несколько random_3d.

Расстояния:

```text
3.5, 3.0, 2.5, 2.2, 2.0, 1.8, 1.6 Å
```

Мультиплетности:

```text
1, 3, 5
```

Команда:

```bash
python3 scripts/generate_b6_inputs.py \
  --project-dir . \
  --stage-dir calculations/stage1 \
  --distances "3.5,3.0,2.5,2.2,2.0,1.8,1.6" \
  --multiplicities "1,3,5" \
  --charge 0 \
  --method R2SCAN-3C \
  --basis "" \
  --task Opt \
  --extra-keywords "TightSCF TightOpt" \
  --nprocs 8 \
  --maxcore 2500 \
  --n-random 5
```

Если `R2SCAN-3C` недоступен в твоей ORCA-сборке, попробуй:

```bash
--method B97-3C
```

или:

```bash
--method PBEH-3C
```

---

## 5. Запуск первичных расчётов

```bash
bash scripts/run_all.sh calculations/stage1
```

Если часть задач упала или не сошлась:

```bash
bash scripts/rerun_failed.sh calculations/stage1
```

---

## 6. Сбор результатов этапа 1

```bash
python3 scripts/collect_results.py \
  --root calculations/stage1 \
  --csv results/results.csv \
  --best-xyz results/best_B6.xyz
```

Просмотр таблицы:

```bash
column -s, -t < results/results.csv | less -S
```

На этапе 1 обычно нет частотного расчёта, поэтому `is_true_minimum` будет `False`. Это нормально.

---

## 7. Этап 2 — подготовка финальных Opt Freq кандидатов

Берём 5 лучших сошедшихся структур этапа 1:

```bash
python3 scripts/prepare_final_candidates.py \
  --project-dir . \
  --results-csv results/results.csv \
  --final-dir calculations/final \
  --n 5 \
  --method PBE0 \
  --basis def2-TZVP \
  --extra-keywords "D4 def2/J RIJCOSX TightSCF TightOpt" \
  --nprocs 8 \
  --maxcore 2500
```

---

## 8. Запуск финальных Opt Freq расчётов

```bash
bash scripts/run_all.sh calculations/final
```

Если нужно:

```bash
bash scripts/rerun_failed.sh calculations/final
```

---

## 9. Финальный сбор результатов

```bash
python3 scripts/collect_results.py \
  --root calculations/final \
  --csv results/final_results.csv \
  --best-xyz results/best_B6.xyz
```

Просмотр:

```bash
column -s, -t < results/final_results.csv | less -S
cat results/best_B6.xyz
```

---

## 10. Критерий настоящего минимума

Структура считается настоящим минимумом только если:

```text
normal_termination = True
optimization_converged = True
has_imaginary_frequencies = False
n_imaginary_frequencies = 0
is_true_minimum = True
```

Если самая низкая по энергии структура имеет мнимые частоты, она не выбирается как `best_B6.xyz`.
Будет выбрана следующая по энергии структура без мнимых частот.

---

## 11. Что делать при мнимой частоте

1. Открыть `.out` или `.hess` в визуализаторе частот.
2. Найти imaginary mode.
3. Сместить структуру вдоль этой моды в плюс и минус направлениях.
4. Сохранить две структуры.
5. Запустить для них новый `Opt Freq`.

Общий принцип:

```text
R_new = R_old ± scale × eigenvector_imaginary_mode
```

Типичный `scale`: `0.05–0.15 Å`.

---

## 12. Полный чек-лист

```bash
cd ~/B6_ORCA_campaign
bash setup_project.sh
export ORCA_CMD=/full/path/to/orca

python3 scripts/generate_b6_inputs.py --project-dir .
bash scripts/run_all.sh calculations/stage1
python3 scripts/collect_results.py --root calculations/stage1 --csv results/results.csv --best-xyz results/best_B6.xyz

python3 scripts/prepare_final_candidates.py --project-dir . --results-csv results/results.csv --n 5
bash scripts/run_all.sh calculations/final
python3 scripts/collect_results.py --root calculations/final --csv results/final_results.csv --best-xyz results/best_B6.xyz

column -s, -t < results/final_results.csv | less -S
cat results/best_B6.xyz
```
