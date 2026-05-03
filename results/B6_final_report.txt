# Многостартовый DFT-поиск устойчивой геометрии нейтрального кластера B₆ методом ORCA 6.1

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

Рисунок 1: `results\figures\Figure_1_workflow.svg`.

## 5. Генерация стартовых геометрий B₆
Для уменьшения риска попадания в локальный минимум был использован многостартовый подход. Были сгенерированы плоские, квазиплоские и трехмерные стартовые структуры B₆ с различными начальными расстояниями B-B. Для каждой структуры были проверены мультиплетности 1, 3 и 5.

В расчетной кампании использовались следующие типы стартов: линейная цепочка; плоское кольцо; компактная плоская структура; ромбическая структура; прямоугольная структура; октаэдрическая 3D-структура; тригональная призма; случайные 3D-структуры. Дополнительно генератор поддерживает искаженное плоское кольцо, fused-triangle, квазиплоскую и пирамидальную 3D-структуру для расширенного набора.

Рисунок 2: `results\figures\Figure_2_start_geometries.svg`.

## 6. Первичный screening: R2SCAN-3C Opt
На screening-этапе обработано `252` ORCA output-файлов. Нормально завершились `243` расчетов, сходимость оптимизации обнаружена у `243` расчетов. Полные энергии извлекались из строки `FINAL SINGLE POINT ENERGY`; расчеты без нормального завершения или без сходимости не рассматриваются как надежные финальные кандидаты.

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

Рисунок 3: `results\figures\Figure_3_screening_top10.svg`.

## 7. Отбор финальных кандидатов
Финальные кандидаты выбирались из низкоэнергетических расчетов screening-этапа с нормальным завершением и сошедшейся оптимизацией. Для уменьшения дублирования структур используется сравнение отсортированных межатомных расстояний B-B; структуры с близкими distance fingerprints рассматриваются как геометрически повторяющиеся кандидаты. В текущем финальном наборе сохранены 10 OptFreq расчетов.

## 8. Финальный расчет: PBE0-D4/def2-TZVP OptFreq
Финальный этап включал `10` расчетов PBE0-D4/def2-TZVP OptFreq. Нормально завершились `10` расчетов, сходимость оптимизации обнаружена у `10` расчетов.

Таблица 2. Финальные расчеты PBE0-D4/def2-TZVP OptFreq.

| calculation_name | multiplicity | method | basis | total_energy_hartree | relative_energy_ev | lowest_frequency_cm-1 | n_imaginary_frequencies | is_true_minimum | xyz_file | output_file |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rank02_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631808591174 | 0.00000000 | 233.780000 | 0 | True | calculations\final\FINAL_rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations\final\FINAL_rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank03_random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631808572899 | 0.00000050 | 234.240000 | 0 | True | calculations\final\FINAL_rank03_B6_random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank03_B6_random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations\final\FINAL_rank03_B6_random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank03_B6_random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank08_rhombic_planar_d3.50_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631808437135 | 0.00000419 | 232.950000 | 0 | True | calculations\final\FINAL_rank08_B6_rhombic_planar_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank08_B6_rhombic_planar_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations\final\FINAL_rank08_B6_rhombic_planar_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank08_B6_rhombic_planar_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank07_rhombic_planar_d1.60_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631808352915 | 0.00000648 | 232.940000 | 0 | True | calculations\final\FINAL_rank07_B6_rhombic_planar_d1.60_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank07_B6_rhombic_planar_d1.60_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations\final\FINAL_rank07_B6_rhombic_planar_d1.60_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank07_B6_rhombic_planar_d1.60_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank01_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631808279344 | 0.00000849 | 234.040000 | 0 | True | calculations\final\FINAL_rank01_B6_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank01_B6_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations\final\FINAL_rank01_B6_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank01_B6_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank06_rhombic_planar_d1.80_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631807854760 | 0.00002004 | 232.880000 | 0 | True | calculations\final\FINAL_rank06_B6_rhombic_planar_d1.80_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank06_B6_rhombic_planar_d1.80_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations\final\FINAL_rank06_B6_rhombic_planar_d1.80_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank06_B6_rhombic_planar_d1.80_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank10_rhombic_planar_d3.00_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631807560673 | 0.00002804 | 232.840000 | 0 | True | calculations\final\FINAL_rank10_B6_rhombic_planar_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank10_B6_rhombic_planar_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations\final\FINAL_rank10_B6_rhombic_planar_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank10_B6_rhombic_planar_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank09_rhombic_planar_d2.50_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631807518628 | 0.00002919 | 232.850000 | 0 | True | calculations\final\FINAL_rank09_B6_rhombic_planar_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank09_B6_rhombic_planar_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations\final\FINAL_rank09_B6_rhombic_planar_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank09_B6_rhombic_planar_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank05_planar_ring_d2.50_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631801739118 | 0.00018645 | 233.060000 | 0 | True | calculations\final\FINAL_rank05_B6_planar_ring_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank05_B6_planar_ring_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations\final\FINAL_rank05_B6_planar_ring_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank05_B6_planar_ring_d2.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |
| rank04_random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C | 3 | PBE0 | def2-TZVP | -148.631801505515 | 0.00019281 | 233.660000 | 0 | True | calculations\final\FINAL_rank04_B6_random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank04_B6_random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq_optimized.xyz | calculations\final\FINAL_rank04_B6_random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq\FINAL_rank04_B6_random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq.out |

Рисунок 5: `results\figures\Figure_5_final_relative_energies.svg`.

## 9. Частотный анализ
Структура считалась истинным минимумом только при одновременном выполнении трех условий: `ORCA TERMINATED NORMALLY`, сходимость оптимизации и отсутствие мнимых частот. Если структура имеет хотя бы одну мнимую частоту, она не считается финальным минимумом даже при низкой электронной энергии.

В ORCA output для финальных расчетов присутствуют шесть нулевых трансляционно-вращательных мод. Они не учитывались как мнимые вибрационные моды и не использовались при выборе `lowest_frequency_cm-1`; в таблицу записана минимальная ненулевая вибрационная частота после отсечения мод с |ν| <= 10 cm⁻¹.

В текущей финальной таблице структур с `n_imaginary_frequencies > 0` не найдено; истинных минимумов по указанному критерию: `10`. Для выбранной структуры `n_imaginary_frequencies = 0`, `lowest_frequency_cm-1 = 233.780000`.

## 10. Обсуждение результатов
Самой устойчивой по финальной энергии оказалась структура `FINAL_rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq`. Она получена из старта `random_3d_seed1000`, имеет мультиплетность `3` и полную энергию `-148.631808591174` Hartree. Ее относительная энергия принята равной `0.00000000` eV.

Близкие по энергии кандидаты присутствуют: `10` финальных структур лежат в пределах 0.01 eV, а `10` структур - в пределах 0.05 eV от минимума. Очень малые различия между несколькими финальными структурами указывают, что разные стартовые геометрии после оптимизации сходятся к одному и тому же или практически идентичному минимуму. Поэтому физически значимым является не различие между этими строками, а устойчивое воспроизведение одной низкоэнергетической плоской структуры из разных стартов.

3D-старты были конкурентоспособными как исходные кандидаты: `4` из `10` финальных расчетов происходят из 3D/random стартов. При этом анализ планарности оптимизированных финальных XYZ показывает, что `4` из них имеют RMS-отклонение от лучшей плоскости не больше 0.01 Å. Для выбранного `best_B6.xyz` RMS-отклонение от плоскости равно `0.0000` Å, максимальное отклонение `0.0000` Å. Поэтому итоговая структура является плоской или практически плоской, несмотря на то что лучший старт был random 3D.

Рисунок 4: `results\figures\Figure_4_best_B6.svg`.

Полученный результат согласуется с литературной тенденцией малых борных кластеров к плоским или квазиплоским структурам. Важно, что этот вывод сделан после проверки 3D-стартов, а не путем их исключения заранее.

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

## Литература
[1] J. Burkhardt, Y. Jia, W.-L. Li. Structure Search with the Strategic Escape Algorithm. Journal of Chemical Theory and Computation, 2025, 21, 3765-3773. DOI: https://doi.org/10.1021/acs.jctc.4c01746
[2] Q.-S. Li, B. Song, L. Wen, L.-M. Yang, E. Ganz. Elucidation of Structures, Electronic Properties, and Chemical Bonding of Monophosphorus-Substituted Boron Clusters in Neutral, Negative, and Positively Charged PBn/PBn-/PBn+ (n = 4-8). Condensed Matter, 2022, 7, 66. https://www.mdpi.com/2410-3896/7/4/66
[3] Milon, D. Roy, F. Ahmed. A DFT study to investigate the physical, electrical, optical properties and thermodynamic functions of boron nanoclusters (MxB2n0; x=1,2, n=3,4,5). Heliyon, 2023, 9, e17886. DOI: https://doi.org/10.1016/j.heliyon.2023.e17886
[4] A. N. Alexandrova, A. I. Boldyrev, H.-J. Zhai, L.-S. Wang, E. Steiner, P. W. Fowler. Structure and Bonding in B6- and B6: Planarity and Antiaromaticity. Journal of Physical Chemistry A, 2003, 107, 1359-1369. DOI: https://doi.org/10.1021/jp0268866
[5] W.-L. Li, X. Chen, T. Jian, T.-T. Chen, J. Li, L.-S. Wang. From planar boron clusters to borophenes and metalloborophenes. Nature Reviews Chemistry, 2017, 1, 0071. DOI: https://doi.org/10.1038/s41570-017-0071
