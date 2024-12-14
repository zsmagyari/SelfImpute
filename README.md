
This repository contain the full dataset associated to a proposed scientific publication entitled
# Experimental Comparative Study on Self-Imputation Methods and their Quality Assessment for Monthly River Flow Data with Gaps: case Study to Mures River
authored by Zsolt Magyari-Sáska, Ionel Haidu and Attila Magyari-Sáska

## Content of folders
* **Originals** - original data with no gaps
* **Masked** - 1000 samples for each station with different [05,10,15,20] gap percentage category each, used in forward imputation
* **Imputed** - Ratio method, Kalman Filter, Random Forest, Gradient Bosst, Extrem Gradient Boost, CatBooster imputed and reimputed datasets for all 4 percentage category for each station
* **Remasked** - remasked datasets created after forward imputation for each imputation method and for each gap percentage category, used in backward imputation

## function.R
- contains 4 functions used in gap creation and imputation performance assessment
- description of usage inside
- the functions are linked to the presented folder structure


