# Data

`nhanes_dr1tot.csv` holds six daily nutrient intakes from the NHANES
2015–2016 Dietary Interview, Total Nutrient Intakes, First Day, file
[`DR1TOT_I.xpt`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/DR1TOT_I.xpt).
The data are published by the U.S. National Center for Health Statistics and
are in the public domain.

| Column | Nutrient                  | NHANES variable |
|--------|---------------------------|-----------------|
| RET    | Retinol (mcg)             | `DR1TRET`  |
| AC     | Alpha-carotene (mcg)      | `DR1TACAR` |
| LZ     | Lutein + zeaxanthin (mcg) | `DR1TLZ`   |
| VA     | Vitamin A, RAE (mcg)      | `DR1TVARA` |
| BC     | Beta-carotene (mcg)       | `DR1TBCAR` |
| VK     | Vitamin K (mcg)           | `DR1TVK`   |

The file has one row per respondent in `DR1TOT_I` (9,544 rows). For 1,217
respondents without a valid first-day dietary recall, all six values are
missing. `scripts/07_data_application.R` drops these rows and uses the
8,327 complete observations.
