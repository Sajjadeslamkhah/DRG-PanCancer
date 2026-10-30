# DRG pan-cancer — v3 revision analyses

Generated 2026-10-30 00:31 | runtime 199.5 min

## R1a Structure null

| test                                         |   observed |   null_median |   null_q025 |   null_q975 |   empirical_p_lower |   n_null |   empirical_p_upper |
|:---------------------------------------------|-----------:|--------------:|------------:|------------:|--------------------:|---------:|--------------------:|
| core-vs-cytoskeleton median rho              |    -0.0357 |        0.0859 |     -0.1566 |      0.3262 |              0.1189 |     1000 |            nan      |
| core coherence (median pairwise rho)         |     0.0427 |        0.001  |     -0.0255 |      0.0646 |            nan      |      200 |              0.0597 |
| cytoskeleton coherence (median pairwise rho) |     0.1223 |        0.0173 |     -0.0118 |      0.0641 |            nan      |      200 |              0.005  |

## R1b Prognostic null (OS)

| score           |   cancers |   obs_n_FDR_sig |   null_median_n_FDR_sig |   null_q95_n_FDR_sig |   p_n_FDR_sig |   obs_pooled_HR |   null_median_pooled_HR |   p_pooled_logHR_abs |   obs_mean_abs_logHR |   null_median_mean_abs_logHR |   p_mean_abs_logHR |
|:----------------|----------:|----------------:|------------------------:|---------------------:|--------------:|----------------:|------------------------:|---------------------:|---------------------:|-----------------------------:|-------------------:|
| DRG_24          |        28 |               9 |                       3 |                    8 |        0.041  |          1.1489 |                  1.0719 |               0.0929 |               0.2243 |                       0.1691 |             0.0509 |
| core_10         |        28 |               0 |                       3 |                    7 |        1      |          1.0688 |                  1.035  |               0.3107 |               0.1863 |                       0.1594 |             0.1858 |
| cytoskeleton_14 |        28 |               4 |                       3 |                    7 |        0.4226 |          1.1143 |                  1.0737 |               0.2178 |               0.1921 |                       0.1688 |             0.2318 |

## R2 Metabolic specificity (median partial rho)

| feature             | pathway        | adjustment           |   median_rho |   n_pos |   n_neg |   k |
|:--------------------|:---------------|:---------------------|-------------:|--------:|--------:|----:|
| DRG_score           | OXPHOS         | purity               |      -0.0606 |       3 |      10 |  33 |
| DRG_score           | OXPHOS         | purity+NRF2          |      -0.1302 |       2 |      13 |  33 |
| DRG_score           | OXPHOS         | purity+NRF2+microenv |      -0.0641 |       2 |      11 |  33 |
| DRG_score           | PPP_NADPH_full | purity               |       0.2948 |      26 |       0 |  33 |
| DRG_score           | PPP_NADPH_full | purity+NRF2          |       0.2642 |      26 |       0 |  33 |
| DRG_score           | PPP_NADPH_full | purity+NRF2+microenv |       0.2625 |      24 |       0 |  33 |
| DRG_score           | oxPPP_strict   | purity               |       0.0557 |       8 |       3 |  33 |
| DRG_score           | oxPPP_strict   | purity+NRF2          |       0.025  |       4 |       4 |  33 |
| DRG_score           | oxPPP_strict   | purity+NRF2+microenv |       0.0087 |       3 |       2 |  33 |
| LRPPRC              | OXPHOS         | purity               |       0.0344 |       7 |       7 |  33 |
| LRPPRC              | OXPHOS         | purity+NRF2          |       0.0566 |       7 |       6 |  33 |
| LRPPRC              | OXPHOS         | purity+NRF2+microenv |       0.0206 |       5 |       6 |  33 |
| LRPPRC              | PPP_NADPH_full | purity               |       0.3022 |      23 |       0 |  33 |
| LRPPRC              | PPP_NADPH_full | purity+NRF2          |       0.3228 |      24 |       0 |  33 |
| LRPPRC              | PPP_NADPH_full | purity+NRF2+microenv |       0.3288 |      27 |       0 |  33 |
| LRPPRC              | oxPPP_strict   | purity               |       0.0717 |      11 |       8 |  33 |
| LRPPRC              | oxPPP_strict   | purity+NRF2          |       0.0531 |      11 |       7 |  33 |
| LRPPRC              | oxPPP_strict   | purity+NRF2+microenv |       0.0178 |       9 |       5 |  33 |
| NDUFS1              | OXPHOS         | purity               |      -0.0106 |       7 |       6 |  33 |
| NDUFS1              | OXPHOS         | purity+NRF2          |      -0.0307 |       6 |       6 |  33 |
| NDUFS1              | OXPHOS         | purity+NRF2+microenv |      -0.0593 |       6 |       9 |  33 |
| NDUFS1              | PPP_NADPH_full | purity               |       0.3013 |      26 |       0 |  33 |
| NDUFS1              | PPP_NADPH_full | purity+NRF2          |       0.2781 |      24 |       0 |  33 |
| NDUFS1              | PPP_NADPH_full | purity+NRF2+microenv |       0.2608 |      26 |       0 |  33 |
| NDUFS1              | oxPPP_strict   | purity               |      -0.0349 |       3 |      10 |  33 |
| NDUFS1              | oxPPP_strict   | purity+NRF2          |      -0.0861 |       1 |      12 |  33 |
| NDUFS1              | oxPPP_strict   | purity+NRF2+microenv |      -0.0809 |       0 |      11 |  33 |
| SLC7A11             | OXPHOS         | purity               |      -0.0428 |       4 |       4 |  33 |
| SLC7A11             | OXPHOS         | purity+NRF2          |      -0.091  |       1 |      12 |  33 |
| SLC7A11             | OXPHOS         | purity+NRF2+microenv |      -0.0818 |       1 |      10 |  33 |
| SLC7A11             | PPP_NADPH_full | purity               |       0.3205 |      26 |       1 |  33 |
| SLC7A11             | PPP_NADPH_full | purity+NRF2          |       0.2151 |      25 |       1 |  33 |
| SLC7A11             | PPP_NADPH_full | purity+NRF2+microenv |       0.2075 |      24 |       1 |  33 |
| SLC7A11             | oxPPP_strict   | purity               |       0.0771 |      11 |       4 |  33 |
| SLC7A11             | oxPPP_strict   | purity+NRF2          |       0.028  |       3 |       5 |  33 |
| SLC7A11             | oxPPP_strict   | purity+NRF2+microenv |      -0.0103 |       2 |       6 |  33 |
| core_subset         | OXPHOS         | purity               |       0.2387 |      22 |       0 |  33 |
| core_subset         | OXPHOS         | purity+NRF2          |       0.2155 |      21 |       0 |  33 |
| core_subset         | OXPHOS         | purity+NRF2+microenv |       0.2222 |      23 |       0 |  33 |
| core_subset         | PPP_NADPH_full | purity               |       0.3983 |      29 |       0 |  33 |
| core_subset         | PPP_NADPH_full | purity+NRF2          |       0.3759 |      28 |       0 |  33 |
| core_subset         | PPP_NADPH_full | purity+NRF2+microenv |       0.3661 |      30 |       0 |  33 |
| core_subset         | oxPPP_strict   | purity               |       0.2136 |      21 |       1 |  33 |
| core_subset         | oxPPP_strict   | purity+NRF2          |       0.1711 |      20 |       2 |  33 |
| core_subset         | oxPPP_strict   | purity+NRF2+microenv |       0.1349 |      15 |       1 |  33 |
| cytoskeleton_subset | OXPHOS         | purity               |      -0.3945 |       1 |      27 |  33 |
| cytoskeleton_subset | OXPHOS         | purity+NRF2          |      -0.3682 |       1 |      29 |  33 |
| cytoskeleton_subset | OXPHOS         | purity+NRF2+microenv |      -0.3117 |       1 |      28 |  33 |
| cytoskeleton_subset | PPP_NADPH_full | purity               |      -0.0378 |       9 |       7 |  33 |
| cytoskeleton_subset | PPP_NADPH_full | purity+NRF2          |       0.0425 |       8 |       7 |  33 |
| cytoskeleton_subset | PPP_NADPH_full | purity+NRF2+microenv |       0.031  |       4 |       4 |  33 |
| cytoskeleton_subset | oxPPP_strict   | purity               |      -0.2157 |       2 |      18 |  33 |
| cytoskeleton_subset | oxPPP_strict   | purity+NRF2          |      -0.2258 |       2 |      18 |  33 |
| cytoskeleton_subset | oxPPP_strict   | purity+NRF2+microenv |      -0.1359 |       0 |      14 |  33 |

## R2c SLC7A11 gene-level null

| gene    | pathway      | adjustment   |   observed_median_rho |   null_median |   null_q95 |   empirical_p |   n_null |
|:--------|:-------------|:-------------|----------------------:|--------------:|-----------:|--------------:|---------:|
| SLC7A11 | oxPPP_strict | purity+NRF2  |                 0.028 |       -0.0518 |      0.229 |        0.3027 |     1000 |

## R3a Subset Cox adjusted for purity + microenvironment (FDR<0.1)

| cancer   | model                 | term   |   n |   events |     HR |     lo |     hi |      p |    FDR |
|:---------|:----------------------|:-------|----:|---------:|-------:|-------:|-------:|-------:|-------:|
| BLCA     | core+purity+microenv  | core   | 395 |      173 | 1.2437 | 1.0663 | 1.4507 | 0.0055 | 0.037  |
| BLCA     | joint+purity+microenv | core   | 395 |      173 | 1.243  | 1.0667 | 1.4486 | 0.0053 | 0.036  |
| HNSC     | core+purity+microenv  | core   | 506 |      213 | 1.2231 | 1.0438 | 1.4333 | 0.0128 | 0.069  |
| KIRC     | core+purity+microenv  | core   | 496 |      168 | 0.7632 | 0.6402 | 0.9099 | 0.0026 | 0.0322 |
| KIRC     | joint+purity+microenv | core   | 496 |      168 | 0.7692 | 0.6457 | 0.9163 | 0.0033 | 0.0345 |
| KIRP     | core+purity+microenv  | core   | 282 |       42 | 1.4951 | 1.0676 | 2.0939 | 0.0193 | 0.0867 |
| LAML     | cyto+purity+microenv  | cyto   |  90 |       60 | 1.5065 | 1.0874 | 2.0873 | 0.0138 | 0.0929 |
| LAML     | core+purity+microenv  | core   |  90 |       60 | 1.5059 | 1.1518 | 1.9689 | 0.0028 | 0.0322 |
| LAML     | joint+purity+microenv | core   |  90 |       60 | 1.534  | 1.1627 | 2.0238 | 0.0025 | 0.0345 |
| LGG      | cyto+purity+microenv  | cyto   | 504 |      125 | 1.7067 | 1.2336 | 2.3613 | 0.0012 | 0.0169 |
| LGG      | joint+purity+microenv | cyto   | 504 |      125 | 1.7034 | 1.2314 | 2.3565 | 0.0013 | 0.0175 |
| LIHC     | core+purity+microenv  | core   | 353 |      126 | 1.3035 | 1.0906 | 1.5579 | 0.0036 | 0.0322 |
| LIHC     | joint+purity+microenv | core   | 353 |      126 | 1.3017 | 1.0886 | 1.5564 | 0.0038 | 0.0345 |
| MESO     | cyto+purity+microenv  | cyto   |  80 |       68 | 2.0686 | 1.4051 | 3.0454 | 0.0002 | 0.0062 |
| MESO     | joint+purity+microenv | cyto   |  80 |       68 | 2.0996 | 1.4377 | 3.0664 | 0.0001 | 0.0033 |
| UVM      | cyto+purity+microenv  | cyto   |  80 |       23 | 2.1496 | 1.2111 | 3.8155 | 0.0089 | 0.0805 |

## R3b Correlations adjusted for purity + microenvironment (median)

| index                  |   median |
|:-----------------------|---------:|
| n                      | 293      |
| rho_core_cyto_adj      |   0.0281 |
| p_core_cyto_adj        |   0.151  |
| rho_cyto_OXPHOS_adj    |  -0.3065 |
| p_cyto_OXPHOS_adj      |   0      |
| rho_cyto_purity        |  -0.2452 |
| rho_cyto_Stromal_CAF   |   0.3974 |
| rho_cyto_Endothelial   |   0.2326 |
| rho_cyto_Leukocyte     |   0.1441 |
| rho_cyto_Smooth_muscle |   0.3184 |
| FDR_core_cyto          |   0.2912 |
| FDR_cyto_OXPHOS        |   0      |

## R4a Glioma molecular adjustment

| cancer   | score   | adjustment                     | covariates_used                       | source              |   n |   events |     HR |     lo |     hi |      p |
|:---------|:--------|:-------------------------------|:--------------------------------------|:--------------------|----:|---------:|-------:|-------:|-------:|-------:|
| LGG      | DRG     | age+sex+IDH                    | age10+male+IDHmut                     | PanCanAtlas subtype | 508 |      123 | 1.2866 | 1.059  | 1.5631 | 0.0112 |
| LGG      | DRG     | age+sex+IDH+codel              | age10+male+IDHmut+codel               | PanCanAtlas subtype | 508 |      123 | 1.2048 | 0.986  | 1.4721 | 0.0684 |
| LGG      | DRG     | age+sex+IDH+codel+grade        | age10+male+IDHmut+codel+grade3        | PanCanAtlas subtype | 507 |      123 | 1.1936 | 0.9703 | 1.4684 | 0.0941 |
| LGG      | DRG     | age+sex+IDH+codel+grade+purity | age10+male+IDHmut+codel+grade3+purity | PanCanAtlas subtype | 500 |      123 | 1.2398 | 1.0039 | 1.5313 | 0.046  |
| LGG      | cyto    | age+sex+IDH                    | age10+male+IDHmut                     | PanCanAtlas subtype | 508 |      123 | 1.4655 | 1.191  | 1.8033 | 0.0003 |
| LGG      | cyto    | age+sex+IDH+codel              | age10+male+IDHmut+codel               | PanCanAtlas subtype | 508 |      123 | 1.354  | 1.0847 | 1.6901 | 0.0074 |
| LGG      | cyto    | age+sex+IDH+codel+grade        | age10+male+IDHmut+codel+grade3        | PanCanAtlas subtype | 507 |      123 | 1.2633 | 1.0109 | 1.5787 | 0.0398 |
| LGG      | cyto    | age+sex+IDH+codel+grade+purity | age10+male+IDHmut+codel+grade3+purity | PanCanAtlas subtype | 500 |      123 | 1.2907 | 1.031  | 1.6158 | 0.026  |
| LGG      | core    | age+sex+IDH                    | age10+male+IDHmut                     | PanCanAtlas subtype | 508 |      123 | 0.9673 | 0.8304 | 1.1268 | 0.6697 |
| LGG      | core    | age+sex+IDH+codel              | age10+male+IDHmut+codel               | PanCanAtlas subtype | 508 |      123 | 0.9744 | 0.8364 | 1.1351 | 0.7389 |
| LGG      | core    | age+sex+IDH+codel+grade        | age10+male+IDHmut+codel+grade3        | PanCanAtlas subtype | 507 |      123 | 1.0092 | 0.8563 | 1.1894 | 0.913  |
| LGG      | core    | age+sex+IDH+codel+grade+purity | age10+male+IDHmut+codel+grade3+purity | PanCanAtlas subtype | 500 |      123 | 1.0342 | 0.8756 | 1.2214 | 0.6923 |
| GBM      | DRG     | age+sex+IDH                    | age10+male+IDHmut                     | PanCanAtlas subtype | 117 |       88 | 1.0677 | 0.8412 | 1.3551 | 0.5902 |
| GBM      | DRG     | age+sex+IDH+codel              | age10+male+IDHmut                     | PanCanAtlas subtype | 117 |       88 | 1.0677 | 0.8412 | 1.3551 | 0.5902 |
| GBM      | DRG     | age+sex+IDH+codel+grade        | age10+male+IDHmut                     | PanCanAtlas subtype | 117 |       88 | 1.0677 | 0.8412 | 1.3551 | 0.5902 |
| GBM      | DRG     | age+sex+IDH+codel+grade+purity | age10+male+IDHmut+purity              | PanCanAtlas subtype | 112 |       84 | 0.9854 | 0.7703 | 1.2606 | 0.9069 |
| GBM      | cyto    | age+sex+IDH                    | age10+male+IDHmut                     | PanCanAtlas subtype | 117 |       88 | 1.1768 | 0.9248 | 1.4976 | 0.1854 |
| GBM      | cyto    | age+sex+IDH+codel              | age10+male+IDHmut                     | PanCanAtlas subtype | 117 |       88 | 1.1768 | 0.9248 | 1.4976 | 0.1854 |
| GBM      | cyto    | age+sex+IDH+codel+grade        | age10+male+IDHmut                     | PanCanAtlas subtype | 117 |       88 | 1.1768 | 0.9248 | 1.4976 | 0.1854 |
| GBM      | cyto    | age+sex+IDH+codel+grade+purity | age10+male+IDHmut+purity              | PanCanAtlas subtype | 112 |       84 | 1.13   | 0.8817 | 1.4481 | 0.3344 |
| GBM      | core    | age+sex+IDH                    | age10+male+IDHmut                     | PanCanAtlas subtype | 117 |       88 | 0.8482 | 0.6851 | 1.05   | 0.1306 |
| GBM      | core    | age+sex+IDH+codel              | age10+male+IDHmut                     | PanCanAtlas subtype | 117 |       88 | 0.8482 | 0.6851 | 1.05   | 0.1306 |
| GBM      | core    | age+sex+IDH+codel+grade        | age10+male+IDHmut                     | PanCanAtlas subtype | 117 |       88 | 0.8482 | 0.6851 | 1.05   | 0.1306 |
| GBM      | core    | age+sex+IDH+codel+grade+purity | age10+male+IDHmut+purity              | PanCanAtlas subtype | 112 |       84 | 0.8039 | 0.6517 | 0.9917 | 0.0415 |

## R4b IDH-stratified

| cancer   | score   | stratum      |   n |   events |     HR |     lo |     hi |      p |
|:---------|:--------|:-------------|----:|---------:|-------:|-------:|-------:|-------:|
| LGG      | DRG     | IDH-wildtype |  94 |       50 | 1.108  | 0.8195 | 1.498  | 0.5052 |
| LGG      | DRG     | IDH-mutant   | 414 |       73 | 1.2788 | 0.9809 | 1.6672 | 0.0691 |
| LGG      | cyto    | IDH-wildtype |  94 |       50 | 1.1852 | 0.8507 | 1.6513 | 0.3153 |
| LGG      | cyto    | IDH-mutant   | 414 |       73 | 1.6369 | 1.2694 | 2.1108 | 0.0001 |
| LGG      | core    | IDH-wildtype |  94 |       50 | 0.9804 | 0.7971 | 1.2057 | 0.851  |
| LGG      | core    | IDH-mutant   | 414 |       73 | 0.6707 | 0.5066 | 0.888  | 0.0053 |
| GBM      | DRG     | IDH-wildtype | 110 |       84 | 1.0839 | 0.8538 | 1.3761 | 0.5081 |
| GBM      | cyto    | IDH-wildtype | 110 |       84 | 1.1699 | 0.9173 | 1.492  | 0.2061 |
| GBM      | core    | IDH-wildtype | 110 |       84 | 0.8866 | 0.7185 | 1.094  | 0.2617 |

## R5a PH tests (q<0.1 shown)

| cancer   | model           | term   |   PH_p |   PH_p_min_any_covariate |   PH_q |
|:---------|:----------------|:-------|-------:|-------------------------:|-------:|
| LGG      | joint core+cyto | cyto   | 0.0029 |                   0.0029 | 0.0803 |

## R5b Piecewise HRs

| cancer   | model               | term   |   cut_years |   HR_early |   lo_early |   hi_early |   p_early |   HR_late |   lo_late |   hi_late |   p_late |
|:---------|:--------------------|:-------|------------:|-----------:|-----------:|-----------:|----------:|----------:|----------:|----------:|---------:|
| CESC     | DRG+clinical+purity | DRG    |      1.6537 |     1.6418 |     1.1408 |     2.3629 |    0.0076 |    1.0493 |    0.7585 |    1.4517 |   0.7713 |
| CESC     | joint core+cyto     | core   |      1.6619 |     1.5802 |     1.1519 |     2.1676 |    0.0046 |    0.9058 |    0.6448 |    1.2726 |   0.5685 |
| LGG      | joint core+cyto     | cyto   |      2.2286 |     2.1826 |     1.74   |     2.7378 |    0      |    1.2647 |    0.9789 |    1.6341 |   0.0724 |
| LIHC     | DRG+clinical+purity | DRG    |      1.0322 |     1.5829 |     1.2017 |     2.0851 |    0.0011 |    0.8333 |    0.6099 |    1.1387 |   0.2524 |
| LIHC     | joint core+cyto     | core   |      1.1376 |     1.5757 |     1.2759 |     1.9458 |    0      |    0.9816 |    0.7591 |    1.2694 |   0.8875 |
| LUAD     | DRG+clinical+purity | DRG    |      1.6715 |     1.5197 |     1.2239 |     1.887  |    0.0002 |    0.9943 |    0.8009 |    1.2344 |   0.9586 |
| LUAD     | joint core+cyto     | cyto   |      1.6619 |     1.3676 |     1.0673 |     1.7525 |    0.0133 |    1.0971 |    0.8982 |    1.34   |   0.3637 |
| READ     | joint core+cyto     | cyto   |      1.6701 |     1.2787 |     0.6898 |     2.3705 |    0.435  |    0.3531 |    0.1575 |    0.7918 |   0.0115 |

## R6b LASSO benchmark

| cancer   |    n |   events | clinical_covariates   |   C_DRG |   C_clinical |   C_clinical_plus_DRG |   AUC_IPCW_1y |   AUC_IPCW_3y |   AUC_IPCW_5y |   deltaC_median |   deltaC_q025 |   deltaC_q975 |   pct_splits_DRG_adds | robust_v2   |
|:---------|-----:|---------:|:----------------------|--------:|-------------:|----------------------:|--------------:|--------------:|--------------:|----------------:|--------------:|--------------:|----------------------:|:------------|
| BLCA     |  404 |      177 | age10+male+stage_num  |  0.5687 |       0.6567 |                0.6515 |        0.5915 |        0.5795 |        0.5525 |          0.0035 |       -0.0412 |        0.0367 |                    54 | False       |
| BRCA     | 1057 |      140 | age10+male+stage_num  |  0.5081 |       0.7394 |                0.7124 |        0.4894 |        0.514  |        0.52   |         -0.0212 |       -0.0843 |        0.0105 |                    12 | False       |
| CESC     |  285 |       71 | age10+stage_num       |  0.561  |       0.6241 |                0.6209 |        0.5527 |        0.5574 |        0.5961 |          0.011  |       -0.1383 |        0.0865 |                    60 | False       |
| COAD     |  421 |       92 | age10+male+stage_num  |  0.5046 |       0.7529 |                0.7034 |        0.5192 |        0.4989 |        0.5    |         -0.0251 |       -0.1037 |        0.0167 |                    16 | False       |
| ESCA     |  178 |       75 | age10+male+stage_num  |  0.478  |       0.6464 |                0.5692 |        0.4491 |        0.4645 |        0.2788 |         -0.0594 |       -0.2569 |        0.007  |                     4 | False       |
| GBM      |  153 |      121 | age10+male            |  0.5017 |       0.6049 |                0.5489 |        0.5107 |        0.5    |      nan      |         -0.0403 |       -0.1695 |        0.0507 |                    18 | False       |
| HNSC     |  519 |      220 | age10+male+stage_num  |  0.5642 |       0.6038 |                0.6268 |        0.5554 |        0.5795 |        0.518  |          0.0272 |       -0.0148 |        0.0684 |                    80 | False       |
| KIRC     |  528 |      174 | age10+male+stage_num  |  0.6448 |       0.763  |                0.7635 |        0.6587 |        0.6503 |        0.6856 |          0.003  |       -0.0388 |        0.0263 |                    58 | True        |
| KIRP     |  267 |       39 | age10+male+stage_num  |  0.7709 |       0.8313 |                0.8114 |        0.8354 |        0.7374 |        0.6218 |         -0.0073 |       -0.1357 |        0.0593 |                    42 | True        |
| LAML     |  149 |       92 | age10+male            |  0.5848 |       0.6616 |                0.6552 |        0.6528 |        0.5438 |        0.7122 |         -0.0117 |       -0.1019 |        0.0811 |                    32 | False       |
| LGG      |  511 |      125 | age10+male            |  0.7538 |       0.7495 |                0.8134 |        0.8495 |        0.7577 |        0.6518 |          0.0656 |       -0.0011 |        0.1137 |                    96 | True        |
| LIHC     |  341 |      116 | age10+male+stage_num  |  0.6079 |       0.615  |                0.6339 |        0.6386 |        0.6015 |        0.5722 |          0.0221 |       -0.0484 |        0.0983 |                    66 | False       |
| LUAD     |  484 |      177 | age10+male+stage_num  |  0.5581 |       0.6567 |                0.6413 |        0.5789 |        0.5743 |        0.5089 |         -0.0177 |       -0.0564 |        0.0394 |                    30 | False       |
| LUSC     |  485 |      210 | age10+male+stage_num  |  0.5014 |       0.5671 |                0.5542 |        0.5033 |        0.5179 |        0.5026 |         -0.0141 |       -0.0649 |        0.0524 |                    26 | False       |
| MESO     |   85 |       73 | age10+male+stage_num  |  0.6525 |       0.5116 |                0.6229 |        0.6488 |        0.8211 |        0.8026 |          0.1082 |       -0.0357 |        0.2748 |                    94 | True        |
| OV       |  301 |      182 | age10+stage_num       |  0.5094 |       0.6186 |                0.5654 |        0.4904 |        0.503  |        0.5043 |         -0.0443 |       -0.1269 |        0.0066 |                     4 | False       |
| PAAD     |  174 |       93 | age10+male+stage_num  |  0.617  |       0.5423 |                0.6154 |        0.6401 |        0.6615 |        0.7342 |          0.0668 |       -0.0287 |        0.1945 |                    94 | False       |
| SARC     |  259 |       98 | age10+male            |  0.531  |       0.5837 |                0.5574 |        0.5738 |        0.5302 |        0.4986 |         -0.0111 |       -0.0854 |        0.0598 |                    38 | False       |
| SKCM     |  407 |      189 | age10+male+stage_num  |  0.5968 |       0.6208 |                0.6586 |        0.6506 |        0.6257 |        0.6271 |          0.0356 |       -0.0113 |        0.0764 |                    80 | False       |
| STAD     |  371 |      146 | age10+male+stage_num  |  0.5    |       0.6372 |                0.6027 |        0.5    |        0.5    |        0.521  |         -0.0239 |       -0.1022 |        0.0166 |                     6 | False       |
| UCEC     |  528 |       87 | age10+stage_num       |  0.6249 |       0.7538 |                0.7355 |        0.6401 |        0.6416 |        0.6255 |         -0.0195 |       -0.0775 |        0.0163 |                    18 | False       |

## R6c Random-gene-set LASSO null

| cancer   |   splits |   observed_median_test_C |   null_median |   null_q95 |   empirical_p |   n_null |
|:---------|---------:|-------------------------:|--------------:|-----------:|--------------:|---------:|
| KIRC     |       10 |                   0.645  |        0.6538 |     0.6883 |        0.7129 |      100 |
| KIRP     |       10 |                   0.7709 |        0.7576 |     0.822  |        0.396  |      100 |
| LGG      |       10 |                   0.7544 |        0.7741 |     0.8081 |        0.802  |      100 |
| MESO     |       10 |                   0.62   |        0.6341 |     0.7146 |        0.6535 |      100 |
