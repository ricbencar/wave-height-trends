# Wave Heights Long-Term Trend Analysis

## Overall Description:

This program analyzes long-term trends in significant wave heights (SWH) from a CSV file ("input.csv"). It applies a series of statistical techniques to reduce methodological biases and produces a detailed report ("report.txt").
![wave-height-trends](https://github.com/user-attachments/assets/9f8322a2-1a37-42e2-bc74-b4fc6b400a78)
The key steps include:

1. **CSV Data Reading and Parsing:**
   - Reads the CSV file containing time series data of wave heights.
   - Extracts the date/time and corresponding SWH values.
   - Discards invalid or improperly formatted rows.

2. **Data Sorting and Chronological Ordering:**
   - Sorts the parsed data in ascending order by year and month.

3. **Seasonal Effect Removal (Deseasonalization):**
   - Calculates the average SWH for each calendar month over the entire record.
   - Subtracts the monthly average from each measurement to remove seasonal cycles.

4. **Grouping into Full Decades:**
   - Groups the deseasonalized data into complete decades (e.g., 1940s, 1950s, etc.), including only decades with data for all 10 years.

5. **Advanced Statistical Trend Analysis:**
   - **Modified Mann–Kendall Test:**
     - Computes the trend statistic (S) using inversion counts.
     - Adjusts the variance for lag-1 autocorrelation to avoid inflated significance.
     - Outputs a standardized Z value to indicate trend significance.

   - **Sen’s Slope Estimator:**
     - Computes the median of all pairwise slopes as a robust estimate of the trend.
     - To prevent excessive memory use for large datasets, if the number of pairs is huge, a fixed number of pairs are sampled randomly.

   - **Seasonal Mann–Kendall Test (Framework Outline):**
     - Outlines how to analyze trends for each month separately to further reduce seasonal biases.

   - **One-Way ANOVA with Tukey HSD Post-hoc Test:**
     - Tests whether the mean deseasonalized SWH differs among full decades.
     - If the F-statistic exceeds a rough threshold (F > 2), Tukey HSD identifies which decades differ.

6. **Ranking of Full Decades:**
   - Computes the average deseasonalized SWH for each full decade.
   - Ranks decades from highest to lowest average SWH.

7. **Report Generation:**
   - Compiles all processing details, statistical test results, and conclusions into a very detailed report ("report.txt").

In summary, by carefully preprocessing the data, removing seasonal effects, adjusting for autocorrelation, and grouping into complete decades, the program minimizes biases and produces robust, evidence-based conclusions about long-term changes in wave heights.

## Compile with:

```bash
g++ -O3 -fopenmp -Wall wave_height_trends.cpp -o wave_height_trends -static -static-libgcc -static-libstdc++
```

## Compilation Details (brief):
- `-O3`: High-level optimizations.
- `-fopenmp`: Enable parallel processing with OpenMP.
- `-Wall`: Enable all warnings.
- `-static, -static-libgcc, -static-libstdc++`: Produce a fully statically linked executable.

## Technical References

   - Hirsch, R.M., Slack, J.R., & Smith, R.A. (1982). Techniques of trend analysis for monthly water quality data. Water Resources Research.

   - Sen, P.K. (1968). Estimates of the regression coefficient based on Kendall's tau. Journal of the American Statistical Association.

   - Tukey, J.W. (1949). Comparing Individual Means in the Analysis of Variance. Biometrics.

   - Box, G.E.P. & Jenkins, G.M. (1976). Time Series Analysis: Forecasting and Control.

   - Cochrane, D., & Orcutt, G.H. (1949). Application of Least Squares Regression to Relationships Containing Auto-Correlated Error Terms.
