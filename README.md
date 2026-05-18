# Wave Heights Long-Term Trend Analysis

![wave-height-trends](https://github.com/user-attachments/assets/9f8322a2-1a37-42e2-bc74-b4fc6b400a78)

## Purpose

This repository contains a standalone C++ command-line program for analysing long-term trends in significant wave height (`swh`) from an `input.csv` time series. The program reads the CSV file, extracts the date and significant wave height, removes the mean seasonal cycle, performs trend and variability tests, ranks decades and months, and writes a complete plain-text report to `report.txt`.

The workflow is intended for long metocean records where the objective is to assess whether significant wave heights show systematic long-term changes after accounting for the dominant monthly seasonal signal.

## Main Capabilities

| Capability | Description |
|---|---|
| CSV header-based reading | Reads `input.csv` using the column names in the header, not fixed column positions. |
| Required variables | Uses `datetime` and `swh` for the full statistical analysis. |
| Additional variables | Accepts the full operational CSV structure with wave, wind and direction columns. Variables not required for this analysis are kept in the file but ignored by the trend engine. |
| Data validation | Skips empty, malformed or non-numeric records and reports the number of skipped rows. |
| Chronological sorting | Sorts all valid records by year and month before processing. |
| Seasonal adjustment | Removes the average monthly seasonal cycle from the SWH series. |
| Decadal grouping | Groups deseasonalized values by calendar decade and excludes the final incomplete decade. |
| Modified Mann–Kendall test | Computes the Mann–Kendall trend statistic and applies a lag-1 autocorrelation variance adjustment when required. |
| Sen’s slope | Estimates the robust long-term slope of the deseasonalized SWH series. |
| Seasonal Mann–Kendall test | Splits the deseasonalized series into monthly sub-series and combines monthly Mann–Kendall statistics. |
| ANOVA analysis | Tests whether mean deseasonalized SWH differs between decade groups and between months. |
| Tukey HSD screening | Performs approximate pairwise decade comparisons when the decadal ANOVA F-statistic exceeds the internal threshold. |
| Text report | Writes a detailed human-readable analysis report to `report.txt`. |

## Input File

The executable expects an input file named exactly:

```text
input.csv
```

The file must be placed in the same directory from which the executable is run.

### Expected CSV Structure

The recommended production header is:

```csv
datetime,swh,mwp,mwd,wind,dwi,u10,v10
```

| Column | Required by this program | Meaning | Used in calculations |
|---|---:|---|---:|
| `datetime` | Yes | Date or date-time string. The program reads the year from characters 1–4 and the month from characters 6–7. | Yes |
| `swh` | Yes | Significant wave height. | Yes |
| `mwp` | No | Mean wave period. | No |
| `mwd` | No | Mean wave direction. | No |
| `wind` | No | Wind speed. | No |
| `dwi` | No | Wind direction. | No |
| `u10` | No | 10 m eastward wind component. | No |
| `v10` | No | 10 m northward wind component. | No |

Only `datetime` and `swh` are mandatory. The other columns may remain in the CSV because they are part of the operational dataset format, but they are not used by this specific wave-height trend analysis.

### Minimal Valid Example

```csv
datetime,swh,mwp,mwd,wind,dwi,u10,v10
1979-01-01 00:00:00,2.15,7.8,285,10.4,292,-4.1,9.6
1979-01-01 01:00:00,2.08,7.6,287,10.1,290,-4.0,9.3
1979-01-01 02:00:00,2.21,7.9,284,10.7,293,-4.2,9.8
```

### Input Parsing Rules

| Rule | Program behaviour |
|---|---|
| Header matching | Column names are matched case-insensitively after trimming whitespace. |
| UTF-8 BOM | A UTF-8 byte-order mark at the start of the header is removed automatically. |
| CSV fields | Quoted fields and escaped quotes are supported. |
| Date parsing | The program expects the date string to begin with `YYYY-MM`, for example `1979-01-01 00:00:00`. |
| Numeric parsing | `swh` must be a finite numeric value. |
| Invalid records | Rows with missing required fields, invalid dates or invalid `swh` values are skipped. |
| Empty records | Empty lines are ignored. |

## Output File

The program writes one output file:

```text
report.txt
```

| Output section | Contents |
|---|---|
| Report header | Analysis period based on the processed calendar-decade range. |
| Introduction | Plain-language explanation of the purpose and statistical approach. |
| Basic decadal statistics | Count, mean and standard deviation of deseasonalized SWH for each decade group. |
| Decade ranking | Full decade groups ranked by average deseasonalized SWH. |
| Modified Mann–Kendall test | Total pairs, S statistic, adjusted variance and Z value. |
| Sen’s slope | Robust slope estimate for the deseasonalized time series. |
| Decadal ANOVA | Between-decade comparison of mean deseasonalized SWH. |
| Tukey HSD screening | Approximate pairwise decade comparison when applicable. |
| Seasonal Mann–Kendall test | Month-by-month Mann–Kendall statistics and combined seasonal statistic. |
| Monthly ranking | Months ranked by average deseasonalized SWH. |
| Monthly ANOVA | Between-month comparison of deseasonalized SWH. |
| Monthly decadal analysis | For each month, decade-level count, mean, standard deviation and ANOVA diagnostics. |
| Final conclusions | Plain-language interpretation and a note that statistical trends do not prove causation. |

## Statistical Workflow

| Step | Operation | Result |
|---:|---|---|
| 1 | Read `input.csv` and locate `datetime` and `swh` columns from the header. | Valid SWH observations with year and month. |
| 2 | Remove invalid rows and sort the remaining records chronologically. | Clean ordered time series. |
| 3 | Compute the mean SWH for each calendar month using the full valid record. | Monthly climatological baseline. |
| 4 | Subtract each monthly mean from the corresponding raw SWH value. | Deseasonalized SWH anomaly series. |
| 5 | Group deseasonalized values by calendar decade. | Decadal samples for ranking and ANOVA. |
| 6 | Apply the modified Mann–Kendall trend test to the full deseasonalized series. | Overall monotonic trend statistic. |
| 7 | Estimate Sen’s slope using pairwise slopes or a random pair sample for very large datasets. | Robust slope in SWH units per year. |
| 8 | Split the deseasonalized series by month and apply the Seasonal Mann–Kendall method. | Combined monthly trend statistic. |
| 9 | Perform decadal and monthly ANOVA calculations. | Difference tests between groups. |
| 10 | Write all statistics, rankings and conclusions to `report.txt`. | Complete text report. |

## Methods Implemented

### Deseasonalization

The program first calculates the average SWH for each calendar month over the complete valid dataset. Each observation is then converted into a deseasonalized anomaly:

```text
swh_deseasonalized = swh_observed - mean_swh_for_same_calendar_month
```

This removes the dominant mean seasonal cycle and allows the long-term tests to focus on persistent changes rather than regular month-to-month variability.

### Modified Mann–Kendall Trend Test

The Mann–Kendall test is a non-parametric test for monotonic trend. The program computes the S statistic using an inversion-counting algorithm based on merge sort, which is efficient for long time series.

| Quantity | Meaning |
|---|---|
| `S` | Mann–Kendall trend statistic. Positive values indicate an increasing tendency; negative values indicate a decreasing tendency. |
| `varS` | Variance of `S`. |
| `r` | Lag-1 autocorrelation coefficient of the deseasonalized series. |
| `Z` | Standardized trend statistic. The report uses `|Z| > 1.96` as the 5% significance screening criterion. |

When the lag-1 autocorrelation is positive, the program inflates the Mann–Kendall variance using:

```text
variance_factor = (1 + r) / (1 - r)
```

This reduces the risk of overstating trend significance in positively autocorrelated series.

### Sen’s Slope Estimator

Sen’s slope is computed as the median of pairwise slopes:

```text
slope(i,j) = (swh_deseasonalized[j] - swh_deseasonalized[i]) / (time[j] - time[i])
```

| Dataset size | Program behaviour |
|---|---|
| Up to 1,000,000 possible pairs | Uses all pairwise slopes. |
| More than 1,000,000 possible pairs | Uses a random sample of 1,000,000 valid pairs. |

The output slope is expressed in SWH units per year.

### Seasonal Mann–Kendall Test

The deseasonalized series is split into 12 monthly sub-series. The Mann–Kendall statistic and variance are computed for each month, then combined into a single seasonal test statistic.

| Monthly component | Description |
|---|---|
| Monthly count | Number of valid deseasonalized records for that month. |
| Monthly `S` | Mann–Kendall statistic for that month. |
| Monthly variance | Variance of the monthly statistic. |
| Monthly `Z` | Standardized monthly statistic. |
| Combined `S` | Sum of valid monthly S statistics. |
| Combined variance | Sum of valid monthly variances. |
| Overall `Z` | Combined seasonal Mann–Kendall statistic. |

### ANOVA and Tukey HSD Screening

The program performs one-way ANOVA calculations for grouped deseasonalized SWH data.

| Analysis | Groups compared | Internal interpretation rule |
|---|---|---|
| Decadal ANOVA | Calendar-decade groups | `F > 2` suggests meaningful differences between decades. |
| Monthly ANOVA | Calendar months | `F > 2` suggests meaningful differences between months. |
| Monthly decadal ANOVA | Decade groups within each month | `F > 2` suggests meaningful differences between decades for that month. |
| Tukey HSD screening | Pairwise decade groups | Applied when the decadal ANOVA F-statistic exceeds 2. |

The Tukey HSD calculation uses an approximate critical value for practical screening. It is useful for identifying which decade pairs are most different, but it should be interpreted as an engineering/statistical screening result rather than a replacement for a full distribution-specific post-hoc analysis.

## Compilation

The program is written in standard C++17 and has no third-party library dependencies. OpenMP support is enabled at compilation.

### Recommended Windows / MSYS2 MinGW-w64 Command

```bash
g++ -std=c++17 -O3 -fopenmp -Wall -Wextra -pedantic wave_height_trends.cpp -o wave_height_trends.exe -static -static-libgcc -static-libstdc++
```

### General Linux Command

```bash
g++ -std=c++17 -O3 -fopenmp -Wall -Wextra -pedantic wave_height_trends.cpp -o wave_height_trends
```

### Compiler Flags

| Flag | Purpose |
|---|---|
| `-std=c++17` | Uses the C++17 language standard. |
| `-O3` | Enables high-level compiler optimisation. |
| `-fopenmp` | Enables OpenMP support. |
| `-Wall` | Enables common compiler warnings. |
| `-Wextra` | Enables additional compiler warnings. |
| `-pedantic` | Enforces stricter standard-compliance diagnostics. |
| `-static` | Produces a more portable static binary where supported. |
| `-static-libgcc` | Statically links GCC runtime support. |
| `-static-libstdc++` | Statically links the C++ standard library. |

## Running the Program

Place these files in the same folder:

| File | Role |
|---|---|
| `wave_height_trends.exe` or `wave_height_trends` | Compiled executable. |
| `input.csv` | Input time series file. |

Then run:

### Windows

```bash
./wave_height_trends.exe
```

### Linux

```bash
./wave_height_trends
```

On successful completion, the console prints:

```text
Analysis complete. Please see report.txt for the detailed report and final conclusions.
```

## Console Messages

| Message | Meaning | Action |
|---|---|---|
| `Error opening input.csv` | The input file was not found in the working directory. | Place `input.csv` beside the executable or run the executable from the correct folder. |
| `CSV file is empty.` | The CSV file has no readable header line. | Check that the file is not empty. |
| `CSV header must contain at least these columns: datetime,swh` | Required columns are missing. | Ensure that the header contains `datetime` and `swh`. |
| `No valid data parsed from CSV. Required usable columns: datetime and swh.` | No row contains both a valid date and a valid SWH value. | Check date format and numeric SWH values. |
| `Warning: skipped ... invalid row(s)` | Some rows were rejected but valid rows remain. | Review the skipped rows if the count is unexpectedly high. |

## Interpretation Guide

| Report value | Main interpretation |
|---|---|
| Positive Mann–Kendall `S` | Later values tend to be higher than earlier values. |
| Negative Mann–Kendall `S` | Later values tend to be lower than earlier values. |
| `|Z| > 1.96` | Trend is significant at the 5% screening level used in the report. |
| Positive Sen’s slope | Increasing deseasonalized SWH trend in units per year. |
| Negative Sen’s slope | Decreasing deseasonalized SWH trend in units per year. |
| High decadal mean | The decade has above-average deseasonalized SWH relative to other decades. |
| ANOVA `F > 2` | The program treats group differences as meaningful enough for further interpretation. |
| Tukey difference above critical difference | The corresponding pair of decades is flagged as significantly different by the internal approximate screening. |

## Practical Notes

| Topic | Recommendation |
|---|---|
| Data length | Longer records provide more reliable trend and decade comparisons. |
| Missing data | The program skips invalid rows but does not infill missing time steps. |
| Time resolution | Hourly, 3-hourly, daily or other regular records may be used, provided the date begins with `YYYY-MM`. |
| Units | SWH units are preserved. If `swh` is in metres, Sen’s slope is reported in metres per year. |
| Direction and wind columns | Directional and wind variables can remain in the CSV but are not part of this SWH-only trend analysis. |
| Causation | Statistical trend detection does not establish the physical cause of the trend. |

## Repository Structure

```text
.
├── wave_height_trends.cpp   # C++17 source code
├── input.csv                # Input data file, provided by the user
├── report.txt               # Output report, generated after execution
└── README.md                # Project documentation
```

## Technical References

| Reference | Relevance |
|---|---|
| Hirsch, R. M., Slack, J. R., & Smith, R. A. (1982). Techniques of trend analysis for monthly water quality data. *Water Resources Research*. | Seasonal Mann–Kendall trend analysis. |
| Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association*. | Sen’s slope estimator. |
| Tukey, J. W. (1949). Comparing individual means in the analysis of variance. *Biometrics*. | Tukey post-hoc comparison concept. |
| Box, G. E. P., & Jenkins, G. M. (1976). *Time Series Analysis: Forecasting and Control*. | Time-series autocorrelation context. |
| Cochrane, D., & Orcutt, G. H. (1949). Application of least squares regression to relationships containing auto-correlated error terms. | Autocorrelation adjustment background. |
