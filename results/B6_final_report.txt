# Многостартовый DFT-поиск устойчивой геометрии нейтрального кластера B₆ методом ORCA 6.1

## Аннотация
В работе выполнен многостартовый DFT-поиск низкоэнергетической геометрии нейтрального кластера B₆. На первом этапе обработано `252` расчетов R2SCAN-3C Opt, из которых `243` завершились нормально и `243` дали сошедшуюся оптимизацию. На втором этапе выполнено `10` финальных PBE0-D4/def2-TZVP OptFreq расчетов. Лучшей найденной структурой в рамках данного набора является `FINAL_rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq` с мультиплетностью `3`, энергией `-148.631808591174` Hartree и минимальной ненулевой частотой `233.780000` cm⁻¹.

Ключевой результат: лучшие кандидаты после финальной оптимизации дают плоскую или практически плоскую структуру. Это согласуется с литературной тенденцией малых борных кластеров к 2D/квазиплоским мотивам, но вывод ограничен использованным набором стартов и выбранным уровнем DFT.

## Краткое содержание результата

| Параметр | Значение |
| --- | --- |
| Система | B₆, neutral |
| Заряд | 0 |
| Проверенные мультиплетности | 1, 3, 5 |
| Screening | R2SCAN-3C Opt |
| Финальный уровень | PBE0-D4/def2-TZVP OptFreq |
| Screening .out | 252 |
| Успешные screening | 243 |
| Финальные OptFreq .out | 10 |
| Истинные минимумы без мнимых частот | 10 |
| best_B6.xyz | results/best_B6.xyz |
| Лучшая энергия, Hartree | -148.631808591174 |
| Лучшая мультиплетность | 3 |
| Минимальная ненулевая частота, cm⁻¹ | 233.780000 |

## 1. Введение
Кластеры бора интересны из-за электронодефицитной природы атома B, делокализованного связывания и выраженной структурной конкуренции между плоскими, квазиплоскими и объемными мотивами. Для B₆ это означает, что результат нельзя надежно получить из одной заранее выбранной геометрии: разные стартовые структуры могут сходиться в разные локальные минимумы или, наоборот, показывать, что объемные старты переходят к плоской области поверхности потенциальной энергии.

## 2. Литературный контекст
Современные подходы к поиску минимумов атомных кластеров используют разнообразные начальные структуры, критерии уникальности и последующую квантово-химическую оптимизацию [1]. В работах по малым борным и борсодержащим кластерам типовой протокол включает DFT-оптимизацию, сравнение полных энергий, анализ электронных свойств и проверку устойчивости через частоты [2,3]. Для чистого B₆ особенно важны работы, где обсуждаются планарность, антиароматичность и химическое связывание B₆/B₆⁻ [4]. Более широкие обзоры по size-selected boron clusters также показывают тенденцию малых борных кластеров к плоским и квазиплоским структурам, связанную с делокализацией σ- и π-связей [5]. Поэтому в данной работе сравнивались не только плоские, но и 3D-старты, чтобы не навязывать геометрию заранее.

## 3. Цель и задачи работы
Цель работы: найти наиболее устойчивую геометрию нейтрального кластера B₆ в рамках заданного набора стартовых структур, мультиплетностей и уровней DFT.

Задачи: сгенерировать набор стартовых структур B₆; выполнить R2SCAN-3C Opt screening; сравнить энергии и сходимость; выбрать низкоэнергетические кандидаты; провести финальный PBE0-D4/def2-TZVP OptFreq; проверить отсутствие мнимых частот; сохранить `best_B6.xyz`.

## 4. Методика расчетов
- Программа: ORCA 6.1
- Кластер: B₆
- Заряд: 0
- Проверенные мультиплетности: 1, 3, 5
- Первичный метод: R2SCAN-3C Opt
- Финальный метод: PBE0-D4/def2-TZVP OptFreq
- Дисперсионная поправка: D4
- Число CPU: 8
- `%maxcore`: 2500 MB

Все энергии, частоты и координаты извлекались только из реальных ORCA `.out` файлов и производных CSV/XYZ файлов. Фиктивные или вручную придуманные значения не использовались.

### 4.1. Логика расчётного workflow
Расчётный проект был организован как последовательная кампания, в которой широкий набор стартовых структур сначала быстро оптимизируется на более дешёвом уровне, а затем низкоэнергетические кандидаты уточняются на более дорогом уровне с частотным анализом. Такой подход снижает риск того, что итог будет зависеть от одной произвольно выбранной геометрии.

1. Генерация стартовых XYZ и ORCA `.inp` файлов.
2. R2SCAN-3C Opt screening для всех стартов, расстояний и мультиплетностей.
3. Сбор `FINAL SINGLE POINT ENERGY`, признаков нормального завершения и сходимости оптимизации.
4. Отбор низкоэнергетических кандидатов с геометрической дедупликацией.
5. Финальный PBE0-D4/def2-TZVP OptFreq.
6. Частотный анализ и выбор `best_B6.xyz` только среди структур без мнимых частот.

### 4.2. ORCA-настройки screening
Ключевая строка screening-расчёта:

```orca
! R2SCAN-3C TightSCF TightOpt Opt

%pal
  nprocs 8
end

%maxcore 2500
```

Дополнительно использовались `MaxIter 500` для SCF и `MaxIter 300` для геометрической оптимизации. Полный шаблон хранится в `templates/stage1_opt_template.inp`.

### 4.3. ORCA-настройки финального этапа
Ключевая строка финального расчёта:

```orca
! PBE0 def2-TZVP D4 def2/J RIJCOSX TightSCF TightOpt Opt Freq
```

На финальном этапе выполнялись одновременно переоптимизация и расчёт частот (`Opt Freq`). Полный шаблон хранится в `templates/final_opt_freq_template.inp`.

### 4.4. Контроль качества парсинга
Для каждого `.out` файла проверялись строки `ORCA TERMINATED NORMALLY`, `THE OPTIMIZATION HAS CONVERGED`, `FINAL SINGLE POINT ENERGY` и, для финального этапа, блок `VIBRATIONAL FREQUENCIES`. В частотном анализе нулевые трансляционно-вращательные моды не интерпретировались как мнимые вибрационные частоты.

Рисунок 1: `results/figures/Figure_1_workflow.svg`.

## 5. Генерация стартовых геометрий B₆
Для уменьшения риска попадания в локальный минимум был использован многостартовый подход. Были сгенерированы плоские, квазиплоские и трехмерные стартовые структуры B₆ с различными начальными расстояниями B-B. Для каждой структуры были проверены мультиплетности 1, 3 и 5.

В расчетной кампании использовались следующие типы стартов: линейная цепочка; плоское кольцо; компактная плоская структура; ромбическая структура; прямоугольная структура; октаэдрическая 3D-структура; тригональная призма; случайные 3D-структуры. Дополнительно генератор поддерживает искаженное плоское кольцо, fused-triangle, квазиплоскую и пирамидальную 3D-структуру для расширенного набора.

В опубликованной расчетной кампании обработано `252` screening output-файла. Текущая расширенная версия генератора поддерживает набор до `384` screening-расчетов за счет дополнительных геометрий и расстояний; поэтому различие между числом `252` в обработанных результатах и `384` в README относится к разным состояниям расчетной кампании, а не к ошибке в таблицах.

### 5.1. Набор стартов и их назначение

| Стартовая геометрия | Тип | Зачем нужна в кампании |
| --- | --- | --- |
| linear_chain | 1D | Проверка вытянутого предела и возможной перестройки в компактную форму |
| planar_ring | 2D | Кольцевой плоский мотив B₆ |
| distorted_planar_ring | 2D | Проверка устойчивости кольца к нарушению симметрии |
| compact_planar_triangle | 2D | Компактный фрагмент треугольной борной сетки |
| rhombic_planar | 2D | Плоский ромбический мотив; важен для сравнения с плоскими минимумами |
| rectangular_planar | 2D | Альтернативный плоский мотив с иной топологией B-B контактов |
| fused_triangles_planar | 2D | Два соединённых треугольника как компактный борный мотив |
| quasi_planar | quasi-2D | Проверка слабого выхода атомов из плоскости |
| octahedral_3d | 3D | Высокосимметричный объёмный конкурент |
| trigonal_prism | 3D | Призматический объёмный конкурент |
| pentagonal_pyramid_3d | 3D | Пирамидальный объёмный старт |
| random_3d_seed* | 3D | Набор случайных, но физически разумных стартов |

В уже обработанном screening-наборе представлены `12` типов `geometry_type`: compact_planar_triangle, linear_chain, octahedral_3d, planar_ring, random_3d_seed1000, random_3d_seed1001, random_3d_seed1002, random_3d_seed1003, random_3d_seed1004, rectangular_planar, rhombic_planar, trigonal_prism.
Начальные расстояния B-B в обработанном наборе: 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 3.5 Å.

Рисунок 2: `results/figures/Figure_2_start_geometries.svg`.

## 6. Первичный screening: R2SCAN-3C Opt
На screening-этапе обработано `252` ORCA output-файлов. Нормально завершились `243` расчетов, сходимость оптимизации обнаружена у `243` расчетов. Полные энергии извлекались из строки `FINAL SINGLE POINT ENERGY`; расчеты без нормального завершения или без сходимости не рассматриваются как надежные финальные кандидаты.

### 6.1. Screening по мультиплетностям
Эта таблица показывает, как распределены расчёты по спиновым состояниям. Энергии разных мультиплетностей сравнивались только после успешной оптимизации на одном уровне теории.

| multiplicity | всего | normal | converged | лучшая E, Hartree | лучшая ΔE, eV | лучший расчёт |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 84 | 83 | 83 | -148.737100311600 | 0.16998512 | random_3d_seed1004_d1.60_q0_m1_R2SCAN-3C |
| 3 | 84 | 80 | 80 | -148.743347149551 | 0.00000000 | random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C |
| 5 | 84 | 80 | 80 | -148.696374246400 | 1.27819781 | random_3d_seed1002_d3.50_q0_m5_R2SCAN-3C |

### 6.2. Screening по типам стартовых геометрий
Таблица ниже нужна не для окончательного выбора минимума, а для контроля многостартового покрытия: она показывает, какие типы стартов давали низкоэнергетические структуры после R2SCAN-3C Opt.

| geometry_type | всего | normal | converged | лучшая E, Hartree | лучшая ΔE, eV | лучший расчёт |
| --- | --- | --- | --- | --- | --- | --- |
| compact_planar_triangle | 21 | 20 | 20 | -148.699661927371 | 1.18873545 | compact_planar_triangle_d2.50_q0_m3_R2SCAN-3C |
| linear_chain | 21 | 18 | 18 | -148.631492001145 | 3.04373365 | linear_chain_d3.50_q0_m1_R2SCAN-3C |
| octahedral_3d | 21 | 21 | 21 | -148.696303344346 | 1.28012715 | octahedral_3d_d2.00_q0_m3_R2SCAN-3C |
| planar_ring | 21 | 19 | 19 | -148.743337426637 | 0.00026457 | planar_ring_d2.50_q0_m3_R2SCAN-3C |
| random_3d_seed1000 | 21 | 21 | 21 | -148.743347149551 | 0.00000000 | random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C |
| random_3d_seed1001 | 21 | 20 | 20 | -148.696214213669 | 1.28255252 | random_3d_seed1001_d3.50_q0_m3_R2SCAN-3C |
| random_3d_seed1002 | 21 | 21 | 21 | -148.743339893109 | 0.00019746 | random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C |
| random_3d_seed1003 | 21 | 21 | 21 | -148.737099279024 | 0.17001322 | random_3d_seed1003_d2.50_q0_m1_R2SCAN-3C |
| random_3d_seed1004 | 21 | 20 | 20 | -148.737100311600 | 0.16998512 | random_3d_seed1004_d1.60_q0_m1_R2SCAN-3C |
| rectangular_planar | 21 | 21 | 21 | -148.728610602439 | 0.40100188 | rectangular_planar_d3.00_q0_m3_R2SCAN-3C |
| rhombic_planar | 21 | 21 | 21 | -148.743335807896 | 0.00030862 | rhombic_planar_d1.80_q0_m3_R2SCAN-3C |
| trigonal_prism | 21 | 20 | 20 | -148.737087566006 | 0.17033195 | trigonal_prism_d2.50_q0_m1_R2SCAN-3C |

### 6.3. Низкоэнергетическая область screening
Таблица 1. Топ-10 screening-результатов; полный набор приведен в `results/results.csv`.

| calculation_name | geometry_type | distance | multiplicity | method | total_energy_hartree | relative_energy_ev | normal_termination | optimization_converged |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C | random_3d_seed1000 | 2.5 | 3 | R2SCAN-3C | -148.743347149551 | 0.00000000 | True | True |
| random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C | random_3d_seed1000 | 3.0 | 3 | R2SCAN-3C | -148.743346682833 | 0.00001270 | True | True |
| random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C | random_3d_seed1000 | 2.2 | 3 | R2SCAN-3C | -148.743344753439 | 0.00006520 | True | True |
| random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C | random_3d_seed1002 | 3.5 | 3 | R2SCAN-3C | -148.743339893109 | 0.00019746 | True | True |
| planar_ring_d2.50_q0_m3_R2SCAN-3C | planar_ring | 2.5 | 3 | R2SCAN-3C | -148.743337426637 | 0.00026457 | True | True |
| rhombic_planar_d1.80_q0_m3_R2SCAN-3C | rhombic_planar | 1.8 | 3 | R2SCAN-3C | -148.743335807896 | 0.00030862 | True | True |
| rhombic_planar_d1.60_q0_m3_R2SCAN-3C | rhombic_planar | 1.6 | 3 | R2SCAN-3C | -148.743335803665 | 0.00030874 | True | True |
| rhombic_planar_d3.50_q0_m3_R2SCAN-3C | rhombic_planar | 3.5 | 3 | R2SCAN-3C | -148.743335803571 | 0.00030874 | True | True |
| rhombic_planar_d2.50_q0_m3_R2SCAN-3C | rhombic_planar | 2.5 | 3 | R2SCAN-3C | -148.743335793150 | 0.00030902 | True | True |
| rhombic_planar_d3.00_q0_m3_R2SCAN-3C | rhombic_planar | 3.0 | 3 | R2SCAN-3C | -148.743335774217 | 0.00030954 | True | True |

Рисунок 3: `results/figures/Figure_3_screening_top10.svg`.

## 7. Отбор финальных кандидатов
Финальные кандидаты выбирались из низкоэнергетических расчетов screening-этапа с нормальным завершением и сошедшейся оптимизацией. Для уменьшения дублирования структур используется сравнение отсортированных межатомных расстояний B-B; структуры с близкими distance fingerprints рассматриваются как геометрически повторяющиеся кандидаты. В текущем финальном наборе сохранены 10 OptFreq расчетов.

Практически это означает, что финальный этап не является повторением всех 252 screening-расчётов. Его задача - уточнить наиболее перспективную часть поверхности потенциальной энергии и проверить, являются ли структуры настоящими минимумами по частотам.

Распределение источников финальных кандидатов:

| Категория старта | число финальных расчётов |
| --- | --- |
| 3D/random | 4 |
| planar | 6 |

## 8. Финальный расчет: PBE0-D4/def2-TZVP OptFreq
Финальный этап включал `10` расчетов PBE0-D4/def2-TZVP OptFreq. Нормально завершились `10` расчетов, сходимость оптимизации обнаружена у `10` расчетов.

Таблица 2. Финальные расчеты PBE0-D4/def2-TZVP OptFreq.

| calculation_name | multiplicity | method | basis | total_energy_hartree | relative_energy_ev | lowest_frequency_cm-1 | n_imaginary_frequencies | is_true_minimum | xyz_file | output_file |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rank02_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631808591174 | 0.00000000 | 233.780000 | 0 | True | calculations/final/FINAL_rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations/final/FINAL_rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank03_random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631808572899 | 0.00000050 | 234.240000 | 0 | True | calculations/final/FINAL_rank03_B6_random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank03_B6_random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations/final/FINAL_rank03_B6_random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank03_B6_random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank08_rhombic_planar_d3.50_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631808437135 | 0.00000419 | 232.950000 | 0 | True | calculations/final/FINAL_rank08_B6_rhombic_planar_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank08_B6_rhombic_planar_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations/final/FINAL_rank08_B6_rhombic_planar_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank08_B6_rhombic_planar_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank07_rhombic_planar_d1.60_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631808352915 | 0.00000648 | 232.940000 | 0 | True | calculations/final/FINAL_rank07_B6_rhombic_planar_d1.60_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank07_B6_rhombic_planar_d1.60_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations/final/FINAL_rank07_B6_rhombic_planar_d1.60_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank07_B6_rhombic_planar_d1.60_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank01_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631808279344 | 0.00000849 | 234.040000 | 0 | True | calculations/final/FINAL_rank01_B6_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank01_B6_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations/final/FINAL_rank01_B6_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank01_B6_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank06_rhombic_planar_d1.80_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631807854760 | 0.00002004 | 232.880000 | 0 | True | calculations/final/FINAL_rank06_B6_rhombic_planar_d1.80_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank06_B6_rhombic_planar_d1.80_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations/final/FINAL_rank06_B6_rhombic_planar_d1.80_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank06_B6_rhombic_planar_d1.80_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank10_rhombic_planar_d3.00_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631807560673 | 0.00002804 | 232.840000 | 0 | True | calculations/final/FINAL_rank10_B6_rhombic_planar_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank10_B6_rhombic_planar_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations/final/FINAL_rank10_B6_rhombic_planar_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank10_B6_rhombic_planar_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank09_rhombic_planar_d2.50_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631807518628 | 0.00002919 | 232.850000 | 0 | True | calculations/final/FINAL_rank09_B6_rhombic_planar_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank09_B6_rhombic_planar_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations/final/FINAL_rank09_B6_rhombic_planar_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank09_B6_rhombic_planar_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank05_planar_ring_d2.50_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631801739118 | 0.00018645 | 233.060000 | 0 | True | calculations/final/FINAL_rank05_B6_planar_ring_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank05_B6_planar_ring_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations/final/FINAL_rank05_B6_planar_ring_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank05_B6_planar_ring_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank04_random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631801505515 | 0.00019281 | 233.660000 | 0 | True | calculations/final/FINAL_rank04_B6_random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank04_B6_random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations/final/FINAL_rank04_B6_random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq/FINAL_rank04_B6_random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |

### 8.1. Топ-5 уникальных геометрических групп
Финальная таблица содержит 10 строк, но строки не обязательно соответствуют 10 независимым минимумам. Для ориентировочной дедупликации ниже финальные XYZ сгруппированы по отсортированным расстояниям B-B с порогом 0.02 Å. Такая таблица помогает отделить физически разные мотивы от повторного попадания в один и тот же минимум.

| group | строк | представитель | m | min ΔE, eV | max ΔE, eV | исходные geometry_type |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 10 | rank02_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C | 3 | 0.00000000 | 0.00019281 | planar_ring, random_3d_seed1000, random_3d_seed1002, rhombic_planar |

### 8.2. Планарность финальных структур
Для оценки того, перешли ли 3D-старты в плоские или квазиплоские минимумы, для каждого финального XYZ была рассчитана лучшая плоскость по координатам атомов. В таблице приведены RMS-отклонение атомов от этой плоскости и максимальное абсолютное отклонение.

| calculation_name | source geometry | ΔE, eV | lowest ν, cm⁻¹ | RMS plane, Å | max plane, Å |
| --- | --- | --- | --- | --- | --- |
| rank02_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C | random_3d_seed1000 | 0.00000000 | 233.780000 | 0.00003 | 0.00004 |
| rank03_random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C | random_3d_seed1000 | 0.00000050 | 234.240000 | 0.00006 | 0.00009 |
| rank08_rhombic_planar_d3.50_q0_m3_R2SCAN-3C | rhombic_planar | 0.00000419 | 232.950000 | 0.00000 | 0.00000 |
| rank07_rhombic_planar_d1.60_q0_m3_R2SCAN-3C | rhombic_planar | 0.00000648 | 232.940000 | 0.00000 | 0.00000 |
| rank01_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C | random_3d_seed1000 | 0.00000849 | 234.040000 | 0.00005 | 0.00006 |
| rank06_rhombic_planar_d1.80_q0_m3_R2SCAN-3C | rhombic_planar | 0.00002004 | 232.880000 | 0.00000 | 0.00000 |
| rank10_rhombic_planar_d3.00_q0_m3_R2SCAN-3C | rhombic_planar | 0.00002804 | 232.840000 | 0.00000 | 0.00001 |
| rank09_rhombic_planar_d2.50_q0_m3_R2SCAN-3C | rhombic_planar | 0.00002919 | 232.850000 | 0.00000 | 0.00000 |
| rank05_planar_ring_d2.50_q0_m3_R2SCAN-3C | planar_ring | 0.00018645 | 233.060000 | 0.00001 | 0.00002 |
| rank04_random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C | random_3d_seed1002 | 0.00019281 | 233.660000 | 0.00003 | 0.00004 |

Рисунок 5: `results/figures/Figure_5_final_relative_energies.svg`.
Рисунок 6: `results/figures/Figure_6_min_energy_geometries.svg` - визуализация 10 финальных оптимизированных геометрий с минимальной энергией.

## 9. Частотный анализ
Структура считалась истинным минимумом только при одновременном выполнении трех условий: `ORCA TERMINATED NORMALLY`, сходимость оптимизации и отсутствие мнимых частот. Если структура имеет хотя бы одну мнимую частоту, она не считается финальным минимумом даже при низкой электронной энергии.

В ORCA output для финальных расчетов присутствуют шесть нулевых трансляционно-вращательных мод. Они не учитывались как мнимые вибрационные моды и не использовались при выборе `lowest_frequency_cm-1`; в таблицу записана минимальная ненулевая вибрационная частота после отсечения мод с |ν| <= 10 cm⁻¹.

В текущей финальной таблице структур с `n_imaginary_frequencies > 0` не найдено; истинных минимумов по указанному критерию: `10`. Для выбранной структуры `n_imaginary_frequencies = 0`, `lowest_frequency_cm-1 = 233.780000`.

### 9.1. Подробная сводка по частотам
Для B₆ всего 18 нормальных мод в декартовом представлении: 6 нулевых/трансляционно-вращательных и 12 ненулевых вибрационных. В таблице приведены первые ненулевые частоты, извлечённые из блоков `VIBRATIONAL FREQUENCIES` финальных ORCA output-файлов.

| calculation_name | нулевых мод | ненулевых мод | min ненулевая, cm⁻¹ | первые ненулевые частоты, cm⁻¹ | n_imag |
| --- | --- | --- | --- | --- | --- |
| rank02_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C | 6 | 12 | 233.78 | 233.78, 319.76, 329.76, 398.76, 516.97, 537.83 | 0 |
| rank03_random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C | 6 | 12 | 234.24 | 234.24, 319.95, 329.78, 398.74, 516.97, 538.02 | 0 |
| rank08_rhombic_planar_d3.50_q0_m3_R2SCAN-3C | 6 | 12 | 232.95 | 232.95, 319.35, 329.07, 398.28, 516.95, 537.72 | 0 |
| rank07_rhombic_planar_d1.60_q0_m3_R2SCAN-3C | 6 | 12 | 232.94 | 232.94, 319.35, 329.08, 398.29, 516.95, 537.72 | 0 |
| rank01_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C | 6 | 12 | 234.04 | 234.04, 319.88, 329.78, 398.77, 517.00, 537.98 | 0 |
| rank06_rhombic_planar_d1.80_q0_m3_R2SCAN-3C | 6 | 12 | 232.88 | 232.88, 319.37, 329.10, 398.26, 516.97, 537.75 | 0 |
| rank10_rhombic_planar_d3.00_q0_m3_R2SCAN-3C | 6 | 12 | 232.84 | 232.84, 319.37, 329.09, 398.24, 516.99, 537.75 | 0 |
| rank09_rhombic_planar_d2.50_q0_m3_R2SCAN-3C | 6 | 12 | 232.85 | 232.85, 319.37, 329.08, 398.24, 516.99, 537.75 | 0 |
| rank05_planar_ring_d2.50_q0_m3_R2SCAN-3C | 6 | 12 | 233.06 | 233.06, 319.72, 330.11, 398.45, 516.96, 537.80 | 0 |
| rank04_random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C | 6 | 12 | 233.66 | 233.66, 319.84, 329.65, 398.47, 517.11, 537.96 | 0 |

## 10. Обсуждение результатов
Самой устойчивой по финальной энергии оказалась структура `FINAL_rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq`. Она получена из старта `random_3d_seed1000`, имеет мультиплетность `3` и полную энергию `-148.631808591174` Hartree. Ее относительная энергия принята равной `0.00000000` eV.

Близкие по энергии кандидаты присутствуют: `10` финальных структур лежат в пределах 0.01 eV, а `10` структур - в пределах 0.05 eV от минимума. При этом `10` финальных строк отличаются от лучшей структуры менее чем на 0.001 eV, то есть намного меньше практической точности обычного DFT-сравнения изомеров. Несмотря на наличие 10 финальных строк, они, вероятно, представляют несколько очень близких или практически идентичных минимумов. Поэтому физически значимый вывод состоит не в различии между отдельными строками, а в устойчивом воспроизведении одной плоской низкоэнергетической структуры из разных стартовых геометрий.

3D-старты были конкурентоспособными как исходные кандидаты: `4` из `10` финальных расчетов происходят из 3D/random стартов. При этом анализ планарности оптимизированных финальных XYZ показывает, что `4` из них имеют RMS-отклонение от лучшей плоскости не больше 0.01 Å. Для выбранного `best_B6.xyz` RMS-отклонение от плоскости равно `0.0000` Å, максимальное отклонение `0.0000` Å. Поэтому итоговая структура является плоской или практически плоской, несмотря на то что лучший старт был random 3D.

Рисунок 4: `results/figures/Figure_4_best_B6.svg`.

Полученный результат согласуется с литературной тенденцией малых борных кластеров к плоским или квазиплоским структурам. Важно, что этот вывод сделан после проверки 3D-стартов, а не путем их исключения заранее.

### 10.1. Геометрические характеристики выбранной структуры
Для `best_B6.xyz` минимальное межатомное расстояние B-B равно `1.516993` Å, максимальное расстояние среди всех 15 пар атомов равно `4.061555` Å, среднее расстояние по всем парам равно `2.239707` Å. Если использовать простой геометрический cutoff 2.05 Å для близких B-B контактов, получается `9` коротких контактов.

Координаты выбранной структуры:

| atom | element | x, Å | y, Å | z, Å |
| --- | --- | --- | --- | --- |
| 1 | B | 0.96858400 | -0.20431000 | 0.99682400 |
| 2 | B | 1.40811300 | -1.44347300 | 0.24014500 |
| 3 | B | 0.20485900 | -0.76279700 | -0.54788700 |
| 4 | B | -0.20486100 | 0.76279800 | 0.54788500 |
| 5 | B | -1.40810400 | 1.44347500 | -0.24014500 |
| 6 | B | -0.96859100 | 0.20430800 | -0.99682200 |

Все попарные расстояния B-B в выбранной структуре:

| pair | distance, Å |
| --- | --- |
| B5-B6 | 1.516993 |
| B1-B2 | 1.516995 |
| B1-B4 | 1.585502 |
| B3-B6 | 1.585503 |
| B4-B5 | 1.591259 |
| B2-B3 | 1.591268 |
| B4-B6 | 1.811440 |
| B1-B3 | 1.811440 |
| B3-B4 | 1.922505 |
| B3-B5 | 2.750271 |
| B2-B4 | 2.750276 |
| B1-B6 | 2.809669 |
| B1-B5 | 3.145462 |
| B2-B6 | 3.145472 |
| B2-B5 | 4.061555 |

## 11. Ограничения расчёта
Следует учитывать, что полученный минимум является лучшим найденным минимумом в рамках использованного набора стартовых структур, мультиплетностей 1, 3 и 5 и выбранного уровня теории PBE0-D4/def2-TZVP. Для более строгого подтверждения глобального минимума можно расширить набор стартовых геометрий, выполнить дополнительную дедупликацию структур, проверить другие функционалы DFT и при необходимости провести более высокоуровневые single-point расчеты.

## 12. Вывод
В рамках выбранного набора стартовых геометрий, проверенных мультиплетностей и уровня теории PBE0-D4/def2-TZVP наиболее устойчивым найденным минимумом нейтрального кластера B₆ является структура `FINAL_rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq`, имеющая мультиплетность `3`, полную энергию `-148.631808591174` Hartree, относительную энергию `0.00000000` eV и не имеющая мнимых частот.

Этот результат не формулируется как доказательство абсолютного глобального минимума B₆; он является лучшим найденным минимумом в рамках выполненного многостартового DFT-набора и выбранного уровня теории.

## 13. Приложения
- `results/results.csv`: полный screening-набор.
- `results/final_results.csv`: финальные OptFreq энергии и частоты.
- `results/best_B6.xyz`: координаты выбранной структуры.
- `results/B6_final_report.txt`: текст отчета.
- `calculations/final/*/*.out`: ORCA output-файлы финальных расчетов.
- `results/figures/Figure_4_best_B6.svg`: изображение финальной структуры.
- `results/figures/Figure_6_min_energy_geometries.svg`: визуализация низкоэнергетических финальных геометрий.

### 13.1. Команды воспроизведения обработки данных
Ниже приведены команды, которые не запускают новые квантово-химические расчёты, а только пересобирают таблицы и отчёт из уже существующих ORCA output-файлов.

```bash
python3 scripts/collect_results.py \
  --root calculations/stage1 \
  --csv results/results.csv \
  --best-xyz results/best_B6.xyz \
  --all-energies-csv results/all_energies.csv

python3 scripts/collect_results.py \
  --root calculations/final \
  --csv results/final_results.csv \
  --best-xyz results/best_B6.xyz

python3 scripts/build_b6_report.py --project-dir .
```

### 13.2. Шаблоны ORCA input
Ниже приведены текущие шаблоны input-файлов, которые используются как документированные примеры настроек. Координаты в реальных `.inp` файлах генерируются отдельно для каждой стартовой структуры.

Screening-шаблон:

```orca
! R2SCAN-3C TightSCF TightOpt Opt

%pal
  nprocs 8
end

%maxcore 2500

%scf
  MaxIter 500
end

%geom
  MaxIter 300
end

* xyz 0 1
B   0.00000000   0.00000000   0.00000000
B   1.80000000   0.00000000   0.00000000
B   0.90000000   1.55884573   0.00000000
B   0.90000000   0.51961524   1.46969385
B  -0.90000000   0.51961524   0.48989795
B   2.70000000   0.51961524   0.48989795
*
```

Финальный OptFreq-шаблон:

```orca
! PBE0 def2-TZVP D4 def2/J RIJCOSX TightSCF TightOpt Opt Freq

%pal
  nprocs 8
end

%maxcore 2500

%scf
  MaxIter 500
end

%geom
  MaxIter 300
end

* xyz 0 1
B   0.00000000   0.00000000   0.00000000
B   1.80000000   0.00000000   0.00000000
B   0.90000000   1.55884573   0.00000000
B   0.90000000   0.51961524   1.46969385
B  -0.90000000   0.51961524   0.48989795
B   2.70000000   0.51961524   0.48989795
*
```

### 13.3. Контрольные файлы

| Файл | Назначение |
| --- | --- |
| results/results.csv | screening-таблица R2SCAN-3C Opt |
| results/all_energies.csv | дублирующая полная таблица энергий screening |
| results/final_results.csv | финальные PBE0-D4/def2-TZVP OptFreq результаты |
| results/best_B6.xyz | координаты выбранного минимума |
| results/B6_final_report.md | подробный отчёт в Markdown |
| results/B6_final_report.txt | текстовая копия отчёта |
| results/figures/*.svg | схемы workflow, геометрий и графики энергий |
| results/figures/Figure_6_min_energy_geometries.svg | топ финальных оптимизированных геометрий по энергии |

## Литература
[1] J. Burkhardt, Y. Jia, W.-L. Li. Structure Search with the Strategic Escape Algorithm. Journal of Chemical Theory and Computation, 2025, 21, 3765-3773. DOI: https://doi.org/10.1021/acs.jctc.4c01746
[2] Q.-S. Li, B. Song, L. Wen, L.-M. Yang, E. Ganz. Elucidation of Structures, Electronic Properties, and Chemical Bonding of Monophosphorus-Substituted Boron Clusters in Neutral, Negative, and Positively Charged PBn/PBn-/PBn+ (n = 4-8). Condensed Matter, 2022, 7, 66. https://www.mdpi.com/2410-3896/7/4/66
[3] Milon, D. Roy, F. Ahmed. A DFT study to investigate the physical, electrical, optical properties and thermodynamic functions of boron nanoclusters (MxB2n0; x=1,2, n=3,4,5). Heliyon, 2023, 9, e17886. DOI: https://doi.org/10.1016/j.heliyon.2023.e17886
[4] A. N. Alexandrova, A. I. Boldyrev, H.-J. Zhai, L.-S. Wang, E. Steiner, P. W. Fowler. Structure and Bonding in B6- and B6: Planarity and Antiaromaticity. Journal of Physical Chemistry A, 2003, 107, 1359-1369. DOI: https://doi.org/10.1021/jp0268866
[5] W.-L. Li, X. Chen, T. Jian, T.-T. Chen, J. Li, L.-S. Wang. From planar boron clusters to borophenes and metalloborophenes. Nature Reviews Chemistry, 2017, 1, 0071. DOI: https://doi.org/10.1038/s41570-017-0071
