# B6_ORCA_campaign: отчет v2 без новых расчетов

Этот отчет версии 2 построен только из уже существующих файлов `calculations/stage1/**/*.out`, `calculations/final/**/*.out`, `results/*.csv` и `results/best_B6.xyz`. Новые ORCA-расчеты не запускались.

## 1. Что изменено относительно v1

- Из финальных `.out` извлечены `<S**2>`, ZPE, энтальпия и Gibbs free energy при 298.15 K.
- `results/final_results.csv` дополнен колонками `s2_actual`, `zpe_hartree`, `e_plus_zpe`, `gibbs_298` и производными относительными энергиями.
- Создан файл `results/best_B6_population.csv` с Mulliken/Loewdin зарядами и спиновыми плотностями выбранной структуры.
- Для всего сошедшегося screening выполнена дедупликация distance-fingerprint с порогом 0.02 A.
- Исправлена интерпретация финального этапа: 10 финальных расчетов сошлись в один минимум; это проверка воспроизводимости, а не ранжирование разных изомеров.

## 2. Версия ПО и среда

- ORCA: `6.1.1 RELEASE`.
- Host из ORCA output: `vmi3233575`.
- ОС запуска расчетов: Linux/VPS по структуре проекта и путям ORCA output (`/root/B6_ORCA_campaign/...`).
- Обработка отчета v2 выполнена локально Python-скриптом `scripts/extract_extras.py`.

## 3. Финальный этап: электронная энергия, ZPE и Gibbs

В таблице ниже `relative_e_plus_zpe_ev` пересчитана отдельно по `E + ZPE` и не заменяет исходную электронную `relative_energy_ev`.

| calculation | m | E, Eh | ZPE, Eh | E+ZPE, Eh | dE(E+ZPE), eV | G298, Eh | <S**2> |
| --- | --- | --- | --- | --- | --- | --- | --- |
| rank08_B6_rhombic_planar_d3.50_q0_m3_R2SCAN-3C | 3 | -148.631808437135 | 0.0199158400 | -148.611892597135 | 0.00000000 | -148.6397525700 | 2.207138 |
| rank07_B6_rhombic_planar_d1.60_q0_m3_R2SCAN-3C | 3 | -148.631808352915 | 0.0199159700 | -148.611892382915 | 0.00000583 | -148.6397523500 | 2.207133 |
| rank06_B6_rhombic_planar_d1.80_q0_m3_R2SCAN-3C | 3 | -148.631807854760 | 0.0199165100 | -148.611891344760 | 0.00003408 | -148.6397513600 | 2.207090 |
| rank10_B6_rhombic_planar_d3.00_q0_m3_R2SCAN-3C | 3 | -148.631807560673 | 0.0199164500 | -148.611891110673 | 0.00004045 | -148.6397512300 | 2.207090 |
| rank09_B6_rhombic_planar_d2.50_q0_m3_R2SCAN-3C | 3 | -148.631807518628 | 0.0199164600 | -148.611891058628 | 0.00004186 | -148.6397511700 | 2.207091 |
| rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C | 3 | -148.631808591174 | 0.0199227900 | -148.611885801174 | 0.00018493 | -148.6397421200 | 2.207087 |
| rank03_B6_random_3d_seed1000_d2.20_q0_m3_R2SCAN-3C | 3 | -148.631808572899 | 0.0199254600 | -148.611883112899 | 0.00025808 | -148.6403925500 | 2.207097 |
| rank01_B6_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C | 3 | -148.631808279344 | 0.0199252300 | -148.611883049344 | 0.00025981 | -148.6403929400 | 2.207069 |
| rank05_B6_planar_ring_d2.50_q0_m3_R2SCAN-3C | 3 | -148.631801739118 | 0.0199193900 | -148.611882349118 | 0.00027886 | -148.6397401800 | 2.207113 |
| rank04_B6_random_3d_seed1002_d3.50_q0_m3_R2SCAN-3C | 3 | -148.631801505515 | 0.0199224300 | -148.611879075515 | 0.00036794 | -148.6397358300 | 2.207093 |

Ключевое значение для выбранной структуры:

- calculation from `best_B6.xyz`: `FINAL_rank02_B6_random_3d_seed1000_d3.00_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq`
- point group from ORCA thermochemistry: `C2h`
- `<S**2>` actual/expected: `2.207087` / `2.000000`
- `E + ZPE`: `-148.611885801174` Eh
- `G(298.15 K)`: `-148.6397421200` Eh
- SOMO energies, eV: `-7.3816;-6.8593`
- HOMO/LUMO/gap, eV: `-6.8593` / `-4.1049` / `2.7544`
- minimum by `E + ZPE` within the same C2h basin: `FINAL_rank08_B6_rhombic_planar_d3.50_q0_m3_R2SCAN-3C_PBE0_def2-TZVP_OptFreq`, `dE(E+ZPE)=0.00000000` eV

## 4. Mulliken/Loewdin population для best_B6

Полная таблица сохранена в `results/best_B6_population.csv`.

| atom | Mulliken q | Mulliken spin | Loewdin q | Loewdin spin |
| --- | --- | --- | --- | --- |
| 1 | 0.048753 | 0.247332 | -0.011340 | 0.170679 |
| 2 | -0.013525 | 0.936048 | 0.092950 | 0.769103 |
| 3 | -0.035214 | -0.183369 | -0.081610 | 0.060226 |
| 4 | -0.035213 | -0.183367 | -0.081617 | 0.060228 |
| 5 | -0.013551 | 0.936034 | 0.092954 | 0.769089 |
| 6 | 0.048748 | 0.247322 | -0.011336 | 0.170675 |

## 5. Найденные уникальные минимумы screening

Дедупликация применена ко всем сошедшимся stage1-структурам: `243` попаданий. Порог distance-fingerprint: `0.02 A`. Сумма попаданий по группам равна `243`.

Полная таблица сохранена в `results/screening_unique_minima.csv`.

| group | representative | hits | m | dE, eV | geometry types |
| --- | --- | --- | --- | --- | --- |
| 1 | B6_random_3d_seed1000_d2.50_q0_m3_R2SCAN-3C | 10 | 3 | 0.00000000 | planar_ring:1;random_3d_seed1000:3;random_3d_seed1002:1;rhombic_planar:5 |
| 2 | B6_planar_ring_d2.00_q0_m3_R2SCAN-3C | 2 | 3 | 0.15130884 | planar_ring:2 |
| 3 | B6_random_3d_seed1004_d1.60_q0_m1_R2SCAN-3C | 13 | 1 | 0.16998512 | random_3d_seed1002:4;random_3d_seed1003:2;random_3d_seed1004:6;trigonal_prism:1 |
| 4 | B6_rectangular_planar_d3.00_q0_m3_R2SCAN-3C | 4 | 3 | 0.40100188 | random_3d_seed1000:1;rectangular_planar:3 |
| 5 | B6_rectangular_planar_d3.50_q0_m3_R2SCAN-3C | 1 | 3 | 0.41850499 | rectangular_planar:1 |
| 6 | B6_random_3d_seed1004_d3.50_q0_m1_R2SCAN-3C | 19 | 1 | 0.69699947 | planar_ring:2;random_3d_seed1000:7;random_3d_seed1002:3;random_3d_seed1004:1;rectangular_planar |
| 7 | B6_rhombic_planar_d2.00_q0_m3_R2SCAN-3C | 2 | 3 | 0.72019386 | rhombic_planar:2 |
| 8 | B6_planar_ring_d2.50_q0_m1_R2SCAN-3C | 1 | 1 | 0.87509235 | planar_ring:1 |
| 9 | B6_rhombic_planar_d2.50_q0_m1_R2SCAN-3C | 7 | 1 | 0.89025074 | planar_ring:1;rhombic_planar:6 |
| 10 | B6_random_3d_seed1004_d3.50_q0_m3_R2SCAN-3C | 3 | 3 | 0.91530824 | random_3d_seed1000:1;random_3d_seed1004:1;trigonal_prism:1 |
| 11 | B6_rectangular_planar_d1.60_q0_m3_R2SCAN-3C | 2 | 3 | 1.11032699 | rectangular_planar:2 |
| 12 | B6_random_3d_seed1002_d2.20_q0_m3_R2SCAN-3C | 13 | 3 | 1.16350933 | random_3d_seed1000:1;random_3d_seed1002:5;random_3d_seed1004:6;trigonal_prism:1 |
| 13 | B6_compact_planar_triangle_d2.50_q0_m3_R2SCAN... | 2 | 3 | 1.18873545 | compact_planar_triangle:2 |
| 14 | B6_rhombic_planar_d1.60_q0_m1_R2SCAN-3C | 1 | 1 | 1.24692411 | rhombic_planar:1 |
| 15 | B6_random_3d_seed1002_d3.50_q0_m5_R2SCAN-3C | 8 | 5 | 1.27819781 | random_3d_seed1000:2;random_3d_seed1002:1;rhombic_planar:5 |
| 16 | B6_octahedral_3d_d2.00_q0_m3_R2SCAN-3C | 23 | 3 | 1.28012715 | octahedral_3d:7;random_3d_seed1000:1;random_3d_seed1001:7;random_3d_seed1002:1;random_3d_seed10 |
| 17 | B6_rhombic_planar_d2.50_q0_m5_R2SCAN-3C | 1 | 5 | 1.28138081 | rhombic_planar:1 |
| 18 | B6_rhombic_planar_d3.50_q0_m5_R2SCAN-3C | 4 | 5 | 1.30912375 | rectangular_planar:3;rhombic_planar:1 |
| 19 | B6_planar_ring_d3.50_q0_m3_R2SCAN-3C | 2 | 3 | 1.42484460 | compact_planar_triangle:1;planar_ring:1 |
| 20 | B6_trigonal_prism_d1.80_q0_m5_R2SCAN-3C | 1 | 5 | 1.53209339 | trigonal_prism:1 |
| 21 | B6_rectangular_planar_d1.60_q0_m1_R2SCAN-3C | 2 | 1 | 1.61932424 | rectangular_planar:2 |
| 22 | B6_random_3d_seed1000_d3.50_q0_m5_R2SCAN-3C | 1 | 5 | 1.65794035 | random_3d_seed1000:1 |
| 23 | B6_compact_planar_triangle_d3.50_q0_m5_R2SCAN... | 3 | 5 | 1.68286421 | compact_planar_triangle:3 |
| 24 | B6_compact_planar_triangle_d1.60_q0_m5_R2SCAN... | 2 | 5 | 1.70446028 | compact_planar_triangle:2 |
| 25 | B6_planar_ring_d1.60_q0_m1_R2SCAN-3C | 2 | 1 | 1.78185909 | planar_ring:2 |
| 26 | B6_random_3d_seed1003_d3.00_q0_m5_R2SCAN-3C | 17 | 5 | 2.00185961 | random_3d_seed1000:3;random_3d_seed1002:5;random_3d_seed1003:2;random_3d_seed1004:6;trigonal_pr |
| 27 | B6_planar_ring_d1.80_q0_m3_R2SCAN-3C | 2 | 3 | 2.10142461 | planar_ring:2 |
| 28 | B6_compact_planar_triangle_d2.00_q0_m3_R2SCAN... | 3 | 3 | 2.23491136 | compact_planar_triangle:3 |
| 29 | B6_random_3d_seed1003_d2.00_q0_m1_R2SCAN-3C | 11 | 1 | 2.27182826 | octahedral_3d:3;random_3d_seed1001:4;random_3d_seed1003:4 |
| 30 | B6_random_3d_seed1002_d1.60_q0_m5_R2SCAN-3C | 2 | 5 | 2.30182379 | random_3d_seed1000:1;random_3d_seed1002:1 |
| 31 | B6_rectangular_planar_d2.00_q0_m3_R2SCAN-3C | 1 | 3 | 2.36979432 | rectangular_planar:1 |
| 32 | B6_compact_planar_triangle_d2.20_q0_m5_R2SCAN... | 2 | 5 | 2.52032598 | compact_planar_triangle:2 |
| 33 | B6_planar_ring_d3.50_q0_m5_R2SCAN-3C | 5 | 5 | 2.55140800 | planar_ring:5 |
| 34 | B6_random_3d_seed1001_d2.50_q0_m5_R2SCAN-3C | 18 | 5 | 2.66309362 | octahedral_3d:7;random_3d_seed1001:6;random_3d_seed1003:5 |
| 35 | B6_planar_ring_d1.60_q0_m5_R2SCAN-3C | 1 | 5 | 2.89840162 | planar_ring:1 |
| 36 | B6_planar_ring_d3.00_q0_m5_R2SCAN-3C | 1 | 5 | 2.95118095 | planar_ring:1 |
| 37 | B6_octahedral_3d_d3.00_q0_m1_R2SCAN-3C | 7 | 1 | 2.99474971 | octahedral_3d:3;random_3d_seed1001:3;random_3d_seed1003:1 |
| 38 | B6_trigonal_prism_d3.50_q0_m1_R2SCAN-3C | 2 | 1 | 3.04367285 | linear_chain:1;trigonal_prism:1 |
| 39 | B6_compact_planar_triangle_d2.20_q0_m1_R2SCAN... | 7 | 1 | 3.06886387 | compact_planar_triangle:7 |
| 40 | B6_rectangular_planar_d1.60_q0_m5_R2SCAN-3C | 1 | 5 | 3.20540705 | rectangular_planar:1 |
| 41 | B6_trigonal_prism_d1.60_q0_m3_R2SCAN-3C | 1 | 3 | 3.37281124 | trigonal_prism:1 |
| 42 | B6_rectangular_planar_d2.20_q0_m5_R2SCAN-3C | 1 | 5 | 3.38224990 | rectangular_planar:1 |
| 43 | B6_trigonal_prism_d1.60_q0_m1_R2SCAN-3C | 1 | 1 | 3.54340925 | trigonal_prism:1 |
| 44 | B6_linear_chain_d1.60_q0_m5_R2SCAN-3C | 8 | 5 | 3.54635517 | linear_chain:8 |
| 45 | B6_trigonal_prism_d2.20_q0_m3_R2SCAN-3C | 4 | 3 | 3.75652156 | trigonal_prism:4 |
| 46 | B6_rectangular_planar_d1.80_q0_m5_R2SCAN-3C | 2 | 5 | 3.83475595 | rectangular_planar:2 |
| 47 | B6_trigonal_prism_d2.20_q0_m1_R2SCAN-3C | 3 | 1 | 4.05967959 | trigonal_prism:3 |
| 48 | B6_linear_chain_d1.60_q0_m3_R2SCAN-3C | 3 | 3 | 4.19813741 | linear_chain:3 |
| 49 | B6_trigonal_prism_d1.60_q0_m5_R2SCAN-3C | 1 | 5 | 4.65383674 | trigonal_prism:1 |
| 50 | B6_trigonal_prism_d2.00_q0_m5_R2SCAN-3C | 3 | 5 | 4.96901762 | trigonal_prism:3 |
| 51 | B6_linear_chain_d3.00_q0_m1_R2SCAN-3C | 1 | 1 | 6.25860710 | linear_chain:1 |
| 52 | B6_linear_chain_d2.50_q0_m1_R2SCAN-3C | 4 | 1 | 6.56973428 | linear_chain:4 |
| 53 | B6_octahedral_3d_d2.20_q0_m1_R2SCAN-3C | 1 | 1 | 7.16771133 | octahedral_3d:1 |
| 54 | B6_linear_chain_d1.60_q0_m1_R2SCAN-3C | 1 | 1 | 8.28174665 | linear_chain:1 |

## 6. Интерпретация финального этапа

10 финальных PBE0-D4/def2-TZVP OptFreq расчетов не являются набором независимых финальных изомеров. После оптимизации они приходят к одному C2h-минимуму с различиями полной энергии порядка микровольт-электронвольт. Поэтому финальный этап следует трактовать как проверку воспроизводимости минимума из разных стартов.

## 7. Точечная группа и геометрия best_B6

ORCA thermochemistry определяет точечную группу выбранной структуры как `C2h` при собственном анализе симметрии. Координаты `best_B6.xyz` также имеют центр инверсии с парными атомами B1/B6, B2/B5 и B3/B4.

Характерные расстояния B-B в `best_B6.xyz`:

| pairs | distance, A | count |
| --- | --- | --- |
| B1-B2, B5-B6 | 1.52 | 2 |
| B1-B4, B2-B3, B3-B6, B4-B5 | 1.59 | 4 |
| B1-B3, B4-B6 | 1.81 | 2 |
| B3-B4 | 1.92 | 1 |
| B2-B4, B3-B5 | 2.75 | 2 |
| B1-B6 | 2.81 | 1 |
| B1-B5, B2-B6 | 3.15 | 2 |
| B2-B5 | 4.06 | 1 |

Все попарные расстояния:

| pair | distance, A |
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

## 8. Подраздел 10.2: сравнение с Alexandrova et al., JPC A 2003

Сравнение выполнено с работой A. N. Alexandrova, A. I. Boldyrev, H.-J. Zhai, L.-S. Wang, E. Steiner, P. W. Fowler, *J. Phys. Chem. A* 2003, 107, 1359-1369, DOI: `10.1021/jp0268866`.
Числа литературы ниже приведены как округленные характерные длины C2h-мотива; для строгой публикационной версии их стоит сверить с печатной таблицей/рисунком статьи.

| quantity | this work | Alexandrova et al. 2003 | comment |
| --- | --- | --- | --- |
| Term/state | triplet, m=3 UKS | triplet neutral B6 | Term label is not assigned by this ORCA output. |
| Point group | C2h | C2h | Matches the C2h motif discussed by Alexandrova et al. |
| Short B-B bonds, A | 1.52 x2 | about 1.52 | B1-B2 and B5-B6 in this work. |
| Side B-B bonds, A | 1.59 x4 | about 1.59 | Four equivalent/near-equivalent side bonds after rounding. |
| Long internal B-B bonds, A | 1.81 x2 | about 1.81 | B1-B3 and B4-B6 in this work. |
| Central B-B distance, A | 1.92 x1 | about 1.92 | B3-B4 in this work. |

## 9. Ограничения

- Спиновый порядок установлен только на уровне screening R2SCAN-3C и затем проверен для выбранных triplet-кандидатов на финальном уровне.
- Мультиреференсность не оценивалась; значения `<S**2>` для UKS приведены как диагностические, но не являются полноценной проверкой характера волновой функции.
- Финальный набор не содержит альтернативных изомеров после оптимизации: он содержит разные старты, пришедшие в один C2h-минимум.
- Литературное сравнение по геометрии сделано по характерным длинам B-B; это не заменяет полный benchmark на одинаковом уровне теории.

## 10. Файлы v2

- `results/B6_final_report_v2.md`: этот отчет Markdown.
- `results/B6_final_report_v2.txt`: текстовая копия отчета.
- `results/final_results.csv`: расширенная финальная таблица.
- `results/best_B6_population.csv`: Mulliken/Loewdin population для best_B6.
- `results/screening_unique_minima.csv`: уникальные stage1-минимумы после дедупликации.
